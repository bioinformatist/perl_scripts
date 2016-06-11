#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $usage = <<_EOUSAGE_;

####################################################################################
# perl duplexesCover.pl --file_list <filename> 
# --file_list	a list of input file names without any suffix
#
####################################################################################

_EOUSAGE_
	;
	

my $file_list = "";

# 呈递参数
&GetOptions("file_list=s" => \$file_list);
# 检查传入的参数
die $usage unless (-s $file_list);

my %gn;

open LIST,'<',$file_list or die "cannot open list file:$!\n"; 
while (<LIST>){
     chomp; 
	 # 对每个独立的文件分别计数
	 undef %gn;
	 # 打印样本名称
	 print "$_\t";
     my $filename = $_."_duplexes"; 
	 my $outname = $_.".cover"; 
     open DUPLEXES,'<',$filename or die "cannot open duplexes file:$!\n";
	 open COVER,'>',$outname or die "cannot create output file:$!\n";
	 # 跳过空文件
	 if(-z $filename){
	     print "$filename is NULL!\n";
		 next;
	} 
	 # 每次循环处理一个duplexes文件
	 while (<DUPLEXES>){
	     chomp; 
		 # 跳过空行
	     next unless /^#/; 
		 # 解析正向读段
	     if ($.%3 == 1){
		     # 捕获读段长度、读段起始位点、读段个数三个信息
	         if (/^#(\d+)\s+(\d+)-.+\((\d+)\)/){
			     # 以位置为键创建一个匿名数组的哈希，匿名数组第一个元素存储该读段在该位置贡献的depth
				 # 遍历该读段覆盖的全部位点进行上述操作
		         $gn{$_}[0] += $3 for ($2..($1 + $2 -1))
		    }
	    }
		 # 解析反向读段
		 elsif ($.%3 == 2){
	         if (/^#(\d+)\s+(\d+)-.+\((\d+)\)/){
			     # 以位置为键创建一个匿名数组的哈希，匿名数组第二个元素存储该读段在该位置贡献的depth
				 # 遍历该读段覆盖的全部位点进行上述操作
		         $gn{$_}[1] += $3 for ($2..($1 + $2 - 1))
		    }
	    }
	} 

	# 通过声明标量上下文语境直接输出hash键的个数，即duplexes（不论正反）在基因组上的覆盖长度
	print scalar keys %gn;
	print "\n";
	
# 遍历哈希元素，对位置按照数值大小排序
	for (sort { $a <=> $b } keys %gn){
		 # 从未出现过的位置赋值为0
		 $gn{$_}[0] =0 unless defined $gn{$_}[0]; 
		 $gn{$_}[1] =0 unless defined $gn{$_}[1]; 
		 # 输出位点、正向depth和反向depth三个信息
		 print COVER "$_\t$gn{$_}->[0]\t$gn{$_}->[1]\n";
	}
	
	close DUPLEXES;
	close COVER;

}

close LIST;

