#!/usr/bin/perl

use 5.012;

my $list = shift;
my $out_file = shift;


my %ccs;

for (glob "*.fa") {
	open CCS, $_ or die "Cannot open CCS file $_:$!\n";
	my $id;
	while (<CCS>) {
		$_ =~ s/[\r\n]+//;
		if (/^>/) {
			($id) = map {m*^>(.+/\d+/)ccs*} $_;
		}else{
			$ccs{$id} .= $_;
		}
	}
}


close CCS;

open LIST, $list or die "Cannot open list:$!\n";
open OUT, ">", $out_file or die "Cannot output results:$!\n";

while (<LIST>) {
	$_ =~ s/[\r\n]+//;
	my ($title) = map {m%^(.+/\d+/).*%} $_;
	say OUT ">".$title."ccs";
	say OUT $ccs{$title};
}

close LIST;
close OUT;

__END__