#!/usr/bin/perl -w
use strict;

### check cmd argument
if(scalar(@ARGV) == 0){
	print STDERR "PLEASE SPECIFY INPUT FILE PATH\n";
	exit(-1);
}

my $progenesisFile = $ARGV[0];
my $filteredProgenesisiFile = $progenesisFile;

### check input file
unless($filteredProgenesisiFile =~ /\.csv$/){
	print STDERR "MAKE SURE INPUT FILE IS A .CSV FILE\n";
	exit(-1);
}
### add suffix to out file
$filteredProgenesisiFile =~ s/\.csv$/_FILTERED.csv/;


my $c=0;
open( IN, "<$progenesisFile" ) || die("Couldn't open file $progenesisFile\n");
open(OUT, ">$filteredProgenesisiFile") || die("Couldn't open file $filteredProgenesisiFile\n");

my $peptideColNb =0;
while ( my $line = <IN> ) {

	chomp($line);
	my @fields = split ",", $line;
	
	if($c == 2){

		print OUT $line."\n";
		
		### get peptide sequnece column number
		while( !($fields[$peptideColNb] =~ "Sequence") ){
			$peptideColNb++;
		}
					
	}elsif($c > 2){
		print OUT $line."\n" if (defined $fields[$peptideColNb] && length($fields[$peptideColNb]) > 2) ;
		
	}else{
		print OUT $line."\n";
	}

	print "CHECKING LINE: $c\r";
	$c++;
	
}

print "CREATED FILE $filteredProgenesisiFile\n";

close(IN);
close(OUT);
exit(0);