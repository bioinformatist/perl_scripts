#!/usr/bin/perl
use warnings;
use 5.012;

my $usage = <<_EOUSAGE_;

############################################################################################
#
#  perl $0 <sam_file> <CCS_in_line> <kept_sam> <abandaned_sam> [distance]
#  example:
#         perl $0 2048_2064.sam merged_ccs.line kept aban 10
#
############################################################################################

_EOUSAGE_
	;

###############################
##   全局变量（含默认设置）  ##
###############################
my $in_sam = shift; # 输入的sam文件
my $in_ccs = shift; # 输入的CCS文件
my $primer = "AAGCAGTGGTATCAACGCAGAGTAC"; # 3' primer
my $polyA = "A" x 25;
my $kept_sam = shift; # 输出保留的SAM
my $aban_sam = shift; # 输出干掉的SAM
my $distance = shift; #最大编辑距离

$distance //= 5; # 默认的窗口序列与引物比对时最大编辑距离
unless (-s $in_sam && -s $in_ccs) {#至少需要两个参数
	die $usage;
}

#################
##  主程序开始 ##
#################
main: {
	# 先读一波CCS文件，把序列id和序列读入hash
	open CCS,'<',$in_ccs or die "can't open the input CCS file $in_ccs:$!\n";

	my %ccs;
	my $id;

	while (<CCS>){
	  $_ =~ s/[\r\n]+//g;
		if(/^>(.+ccs)/){
			$id = $1;
		}else{
			my $query = $_;
			$ccs{$id} = $query;
		}
	}


	close CCS;

	open SAM, $in_sam or die "Can't open 2048_2064 file!\n";
	open KEPT, ">", $kept_sam or die "Can't output result!\n";
	open ABAN, ">", $aban_sam or die "Can't output result!\n";

	my $count;
	my $kept;

	while (<SAM>){
		$count ++;
		$_ =~ s/[\r\n]+//g;
		my @F = split /\t/, $_;
		(my $qname = $F[0]) =~ s/^(m.+\/)\d+_\d+_CCS/$1ccs/;
		# 只要当前sam对应的CCS序列中有嵌合或整合错误的迹象（除去开头结尾60bp，仍可见疑似接头或polyA就马上干掉）
		if ((primer_match($ccs{$qname}, $primer)) || (primer_match($ccs{$qname}, rev_com($primer))) || (primer_match($ccs{$qname}, $polyA)) || (primer_match($ccs{$qname}, rev_com($polyA))) || primer_find($ccs{$qname})) {
			say ABAN $_;
		}else{
			$kept++;
			say KEPT $_;
		}
	}
	
	close SAM;
	close KEPT;
	close ABAN;
	
	# 输出统计信息
	my $ratio = $kept / $count;
	say "All records: $count";
	say "Kept ratio: $ratio";
}

##############
##  子程序  ##
##############
# 计算窗口长度的ccs序列与判定序列的汉明距离
sub primer_match{
    # 每次只能考虑read第i个位置开始的固定长度的子序列，与primer匹配
	my $read_seq = shift;   # 得到第1个参数，read序列
	my $primer_seq = shift; # 得到第2个参数，判定序列
	my $read_len = length($read_seq);
	my $window_size = 25; # 滑动窗口长度
	my $offset = 60; # 不考虑前后60nt

    # 采用hamming距离，滑动窗口提取固定长度片段与primer比较
	for(my $i = $offset; $i < ($read_len - $offset - $window_size + 1); $i ++){  # 不考虑前面60nt，依次按照窗口25长向右滑动
		my $read_substr = substr($read_seq, $i + 1, $window_size);# 截取window_size长度片段
		my $editDistance = hamming($read_substr,$primer_seq);# 求截取的片段与引物的汉明距离
		if ($editDistance <= $distance){ # 如果汉明距离小于指定阈值
			return 1;
			last;          # 找到了，就不继续找了，跳出循环
		}
	}
}

# 子程序：求2个字符串之间的hamming距离
sub hamming($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }

# 子程序：反向互补序列
sub rev_com {
	my $seq = shift;
	$seq =~ s/[\r\n]+//;
	my @seq = split qq{},$seq;
	@seq = reverse @seq;
	$seq = join q{},@seq;
	$seq =~ tr/ATCGatcg/TAGCtagc/;
	return $seq;
}

sub primer_find{
	my $read_seq = shift;
	my $read_substrL = substr($read_seq, 0, 25);
	(my $read_substrR) = map{/.+[AGCT]{25}$/} $read_seq;
	my $edL5 = hamming($read_substrL, rev_com($primer));
	my $edR5 = hamming($read_substrR, rev_com($primer));
	my $edL3 = hamming($read_substrL, $primer);
	my $edR3 = hamming($read_substrR, $primer);
	unless ((($edL5 <= $distance) && ($edR3 <= $distance)) || (($edR5 <= $distance) && ($edL3 <= $distance))){
		return 1;
	}
}
