#!/usr/bin/env perl
use strict;
use diagnostics;
use Data::Dumper;
use Getopt::Long;
use 5.012;

my $usage = <<_EOUSAGE_;

#########################################################################################
# duplexesFinder.pl --list <filename> --min_len <num> --max_len <num> --overhang <num> --count <num> [--perfect]
# --list	Sam files to be found duplexes in
# --min_len	The minimum of read length
# --max_len	The maxmum of read length
# --overhang	The 5'-distance of read pairs (can be 0)
# --count	The least number of read counts
# --perfect	Only reserve perfect match reads
# Example:	./duplexesFinder.pl --list sam.list --min_len 21 --max_len 24 --overhang 2 --count 2 --perfect
###########################################################################################

_EOUSAGE_

	;

# 声明全局变量
our $list;
our $min_len;
our $max_len;
our $overhang;
our $count;
our $perfect;

# 呈递参数	
&GetOptions( 'list=s' => \$list,
			 'min_len=s' => \$min_len,
			 'max_len=s' => \$max_len,
			 'overhang=i' => \$overhang,
			 'count=s' => \$count,
			 'perfect!' => \$perfect,
			 );

# 检查传入的参数
die $usage unless ($list && $min_len && $max_len && $count);
			 
open LIST, $list or die "Can't read the list file!\n";

my $displayDistance = $overhang +1;

while (<LIST>) {
	chomp;
	my $samFile = $_.".sam";
	say $samFile;
	open SAM, $samFile or die "Can't open the sam file!\n";
	my $outFile = $_."_duplexes";
	open OUT, ">", $outFile or die "Can't output!\n";

	# 添加两个hash用于存储读段信息
	my %posSeq;
	my %negSeq;

	# 逐行读入sam文件
	while(<SAM>){
		# 去换行符
		chomp;
		# 切片并提取以下4列的信息
		my ($flag, $loc, $cigar, $seq, $nm) = (split(/\t/, $_))[1, 3, 5, 9, 13];
		# 直接滤掉不匹配的读段
		next if ($flag & 0x4);
		# 将NM tag中表示编辑距离的整数提取出来
		(my $editD = $nm) =~ s/NM:i:(\d+)/$1/;
		# 如果是perfect模式，把有错配和带有gap的都去掉
		if ($perfect) {
			next unless ($editD == 0);
		}
		# 提取序列长度
		my $len = length($seq);
		# 跳过那些不关心的长度
		next if (($len < $min_len) or ($len > $max_len));
		# 将长度和起始位点作为hash键
		my $len_loc = (join "\t", ($len,$loc));
		# 10为反向，其余为正向序列
		if ($flag & 0x10) {
			# 使用nested hash储存序列并计数，比如：
=pod
			$VAR1 = {
					  '22	213' => {
									'TTGTGTGATTAAATTAATTCAA' => 12
								  },
					  '21	4007' => {
									 'CTTGAACCTGATGAAGATTTC' => 59
								   },
					}
=cut
			$negSeq{$len_loc}{$seq} ++;
		}
		else {
			$posSeq{$len_loc}{$seq} ++;
		}
	}
	
	# 观察复杂数据结构用
	# print Dumper(\%posSeq);
	# print Dumper(\%negSeq);
	
	
	# 遍历正向外层hash
	for my $pMark (sort keys %posSeq) {
		# 切片，提取正向读段长度和起始位点
		my ($pLen, $pLoc) = split /\t/,$pMark;
			# 使用外层hash的键解内层hash的引用并提取内层hash的键即序列（虽然hash只有一对键值也要用下标转换成标量）
			my $pSeq =  (keys %{ $posSeq{$pMark} })[0];
			# 根据上一行取出的键从内层hash取值即个数
			my $pCount = $posSeq{$pMark}->{$pSeq};
			# 如果正向深度太低，就不考虑反向了
			if ($pCount > $count) {
				# 根据出头长度计算反向序列位点
				my $nLoc = $pLoc - $overhang;
				# 生成反向外层hash的键
				my $nMark = (join "\t", ($pLen,$nLoc));
				# 判断外层键在反向hash中是否存在
				if(exists $negSeq{$nMark}) {
					my $nSeq =  (keys %{ $negSeq{$nMark} })[0];
					my $nCount = $negSeq{$nMark}->{$nSeq};
					# 反向深度过滤
					if ($nCount > $count) {
						# 为便于观察，位点差一位的用0补齐
						$nLoc = "0".$nLoc unless (length($nLoc) == length($pLoc));
						say OUT qq{\#$pLen}.qq{}.qq{ }x($displayDistance).qq{$pLoc-$pSeq($pCount)};
						say OUT qq{\#$pLen $nLoc-$nSeq($nCount)};
						print OUT "\n";
					}
				}	
			}
	}
	close SAM;
	close OUT;
}

close LIST;