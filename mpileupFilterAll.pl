#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
=pod
my $usage = <<_EOUSAGE_;
###########################################################################################
# perl $0 --inputfile <FILE> --refRetained [0|1] --isFiltered [0|1] --minNum [INT] --keptNum [INT] --out [FILE]
# Required(1):
#  --inputfile         The name of a input mpileup file
# Options(2):
#  --refRetained  Set 0 for retaining reference base, 1 for removing reference base [1]
#  --isFiltered   Set 0 for no further filteration , 1 for further filteration [1]			
#  --minNum       Minimum of a mutation depth [5]  
#  --keptNum      Number of retaining mutation types [2]  
#  --out          The name of a output file [$inputfile.snp]    
#
#########################################################################################
_EOUSAGE_
	;
=cut
##############################
##  全局变量（含默认设置）  ##
##############################
our $inputfile = shift;
our $minNum1 = 2; # 最小覆盖数目，默认为5
our $minNum2 = 5;
our $keptNum = 2; # 各个SNP类型中保留的最大数目，默认为2
our $refRetained = 0; # 默认不保留参考碱基
our $isFiltered = 0; # 默认进一步过滤变异类型，只保留前$keptNum个
our $out; # 指定输出文件,若不指定，默认输出文件名为$inputfile.snp
$out = $inputfile.".snp" unless $out;
our @pil;
our $refBase;
=pod
####################
##  输入参数处理  ##
####################
&GetOptions( 'inputfile=s' => \$inputfile,
		'refRetained=i' => \$refRetained,
		'isFiltered=i' => \$isFiltered,
		'minNum=i' => \$minNum,
		'keptNum=i' => \$keptNum,
		'out=s' => \$out
		    );

die $usage unless (-s $inputfile);
=cut

##################
##  主程序开始  ##
##################
main: {
    open IN,'<',$inputfile or die "Can't open mpileup file\n";
    open OUT,'>',$out or die "Can't create outSNP file\n"; 
    while (<IN>){
	    chomp;
	    # 获取染色体名称、参考碱基的位置、参考碱基、正常样本覆盖碱基、异常样本覆盖碱基
	    @pil = split /\t/,$_;
		# 跳过depth低于阈值的行
		if (($pil[3] < 5) || ($pil[6] < 5)){
		    next;
		}
	    # 将参考碱基转化为大写
	    $refBase = uc($pil[2]);
	    #my $sample1Base = uc($pil[4]);
	    #my $sample2Base = uc($pil[7]);
        # 统计正常样本
	    my $normal = &genetyping($pil[3],$pil[4],$pil[5], $minNum1);
	    # 统计异常样本
	    my $abnormal = &genetyping($pil[6],$pil[7],$pil[8], $minNum2);
	    # 比较正常样本和异常样本
	    &comp($normal,$abnormal);
    }
    close IN;
    close OUT;
}
##############
##  子程序  ##
##############
# 解析覆盖该位点的每条reads与该位点的匹配方式，提取序列
sub genetyping {
    my $depth = shift; # 某位点reads覆盖的数目，仅仅包括点、逗、单碱基替换以及删除的覆盖数
    my $read_bases1 = shift; # 覆盖该位点的每条reads与该位点的匹配方式
	my $qual = shift; # 覆盖该位点碱基的碱基质量，也是仅仅包括点、逗、单碱基替换以及删除的情形，不包括插入、缺失、删除、^、$等
	my $thold = shift;
	
	$read_bases1 =~ s/\.|\,/$refBase/g; # 将点逗转化为参考碱基
	my $read_bases = uc($read_bases1);
	my $mapStr_len = length($read_bases);
    my @bases = split //,$read_bases;
	my ($Inbase,$Delbase) = ("","");
	my $min_qual = 20; # # 记录覆盖位点对应的碱基质量低于阈值的不考虑
	my $qualsnums = &convertQualAsciiToNumsPhred33($qual);
	my $str_type = "";
	my @qual = (); # 定义数组，存储覆盖该位点碱基的碱基质量
	my %type = (); # 定义哈希，存储类型与次数
	my $j = -1; # 记录覆盖位点对应的碱基质量的位置，实际中从0计数，方便索引
	foreach(@$qualsnums) {
        push @qual,$_;
	}
	for (my $i=0;$i < $mapStr_len;$i++){ # 每次读取一个字母
	    my $letter = $bases[$i];
		if (($letter eq q{^})){ # 如果出现“^”，后面紧跟的ASCII码字符减去33即为mapping质量，不考虑
			$i = $i+1; # 切换至ASCII码字符位置，注意特殊情况^+. ^C.等等
		}elsif ($letter eq q{$}){ # 如果出现“$”
			;
		}elsif ($letter eq q{+}){ # 分析插入子串，例如+20TTTTTTTTTTTT....TG。只要是插入，全部算作一种类型
			($Inbase,$i) = &parseInDelStr($read_bases,$i);
			$str_type = "+".$Inbase;
			$type{$str_type}++;
		}elsif ($letter eq q{-}){ # 分析缺失子串，例如-20TTTTTTTTTTTT....TG。只要是缺失，全部算作一种类型
			($Delbase,$i) = &parseInDelStr($read_bases,$i);
			$str_type = "-".$Delbase;
			$type{$str_type}++;
		}else{ # 如果出现“*”、参考碱基、单碱基替换
		    $j++; # 碱基质量计数加1
			if ($qual[$j] >= $min_qual){
			    $type{$letter}++;
			}
		}
	}
    for (keys %type){ 
	    if ($type{$_} < $thold){
		    delete $type{$_}; # 删除低于阈值的键值对
		}
	}
	# 判断是否输出参考碱基
	if ($refRetained){ 
	    delete $type{$refBase}; # 不考虑参考碱基
	}
	my $len = keys %type; # 哈希中键的个数
	my @key_sorted = sort{$type{$b} <=> $type{$a}} keys %type;# 对哈希值从大到小排序
	# 判断是否进一步筛选变异类型
	if ($isFiltered){
	    if ($len > $keptNum){ 
	        delete $type{$_} for @key_sorted[$keptNum..$#key_sorted]; # 默认只保留前两名
	    }else{
	        ;
	    }
	}
	#print Dumper(\%type);
	return \%type;
}

# 这个子程序用于解析插入或缺失子串的情形
sub parseInDelStr {
    my $reads = shift; # 输入read与该位点匹配方式字符串
    my $loc = shift; # 输入“+”或“-”的位置（0-based）
	my @ta = split //,$reads;
	$loc++; # 指到"+"或"-"下一个字符，即数字的第一位
	my $strInDel = ""; # 定义一个字符串，用于存放"+"后面的数字
    while ($ta[$loc] =~ m/^(\d)$/){ # 如果是数字
	    $strInDel .= $1;
		$loc++; # 直到遍历完“+”或“-”后的数字
	}
	# 循环结束后，$loc指向数字后面所有碱基的第1个
	my $baseInDel = substr($reads,$loc,$strInDel); # 截取“+数字”后面部分
	my $end = $loc + $strInDel - 1; # 插入/缺失子串的最后一位
	return ($strInDel.$baseInDel,$end); # 返回插入/缺失子串以及插入子串的最后一位坐标
}

# 这个子程序用于计算碱基质量，默认为Phred33质量的测序数据
sub convertQualAsciiToNumsPhred33 {
    my $quality = shift;
    my @nums;
    my @ascii = split //,$quality;
    foreach(@ascii) {
	    push(@nums,(ord($_) - 33)); # ord()用于返回ASCII数值
    }
    return \@nums;
}

# 这个子程序用于比较mpileup文件两个样本genetyping的后者较前者而言的不同之处
sub comp {
    my ($type1,$type2) = @_;
	# 提取第一个样本返回的类型
	my @sampe1_keys = keys %$type1;
	my $sample1 = join("/",@sampe1_keys);
	# 提取第二个样本返回的类型
	my @sampe2_keys = keys %$type2;
	my $sample2 = join("/",@sampe2_keys);
	foreach my $key2 (keys %$type2){ # 遍历第二个样本所有键
		unless (exists $type1->{$key2}){ # 如果在第一个样本中不存在
		    if ($sample1){ # 样本1要满足最小depth
			    print OUT "$pil[0]\t$pil[1]\t$refBase\t/$sample1\t/$sample2\n"; # 输出该位点的突变类型信息
				last; # 上面一旦打印即打印所有，防止重复输出
			}
		}
    } 
}
