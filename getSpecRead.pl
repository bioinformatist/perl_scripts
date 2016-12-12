#!/usr/bin/perl

use 5.012;


my $list = shift;
my $out_file = shift;
my $base_per_line = shift;


my %read;


for (glob "*.fa") {
	open READ, $_ or die "Cannot open reads file $_:$!\n";
	my $id;
	while (<READ>) {
		$_ =~ s/[\r\n]+//;
		if (/^>/) {
			($id) = map {m*^>(.+/\d+)/.+$*} $_;
			$read{$id} = $_."\n";
		}else{
			$read{$id} .= $_;
		}
	}
}


close READ;


open LIST, $list or die "Cannot open list:$!\n";
open OUT, ">", $out_file or die "Cannot output results:$!\n";


while (<LIST>) {
	# Skip blank line
	next if /^\s*$/;
	$_ =~ s/[\r\n]+//;
	my ($title) = map {m%^(.+/\d+)/.*%} $_;
	if (exists $read{$title}) {
		my @single_read = split /\n/, $read{$title};
		say OUT $single_read[0];
		$single_read[1] =~ s/(.{$base_per_line})/$1\n/gs;
		say OUT $single_read[1];
	}
}


close LIST;
close OUT;

__END__