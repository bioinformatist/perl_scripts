#!/usr/bin/perl
use 5.10.0;
use Data::Dumper;

=pod
这个程序有3种思路：
1. 最简单：几乎不使用RAM，开启句柄后每查找完一条，使用seek定位到文件头，既慢又伤硬盘（读写频繁）。缺点：由于不是在RAM内操作，速度极慢
2. 稍微复杂一点：用一层hash，每次遍历hash的键，检查是不是和当前读入的CCS一致，一致则输出。缺点：遍历hash的键，效率较低，相当于遍历一个array，没有发挥出hash的威力
3. 最复杂：$0使用的逻辑。采用%hash{key}->{key}->value的结构将CCS和subreads紧密联系，准确定位，只读一次硬盘，内存效率也最高
=cut

my $ccs_file = shift;
my $sub_file = shift;
my $output = shift;

open CCS, $ccs_file or die "No CCS file:$!\n";
open SUB, $sub_file or die "No subreads file:$!\n";
open OUT, ">", $output or die "Cannot output:$!\n";


my %hash;
my $id;
my $id_header;

while (<SUB>) {
	$_ =~ s/[\r\n]+//;
	if (/^>/) {
		$id = $_;
		($id_header) = map {/(^>.+\/\d+\/).+$/} $_;
	} else {
		$hash{$id_header}{$id} .= $_;
	}
}

# 调试可用如下语句，将STDOUT重定向至文件即可
# print Dumper (\%hash);
close SUB;

while (<CCS>) {
	$_ =~ s/[\r\n]+//;
	if (/^>/) {
		$id = $_;
		($id_header) = map {/(^>.+\/\d+\/)ccs$/} $_;
	} else {
		say OUT $id;
		say OUT $_;
		# 此处注意，$hash{key}->{key}此类unblessed hash reference是不允许在列表上下文中使用的，要使用显式的解引用写法
		for (keys %{ $hash{$id_header} }) {
			say OUT $_;
			say OUT $hash{$id_header}->{$_};	
		}
		say OUT "";
	}
}

close CCS;
close OUT;

__END__