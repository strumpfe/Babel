#!/usr/bin/perl
# FastaComposition.pl

# Usage: FastaComposition.pl <filename>

use strict;
use Getopt::Long;
use carp;

my $file;
my $seq = "";
my $seq2 = "";
my $name;
my $total;
my $total_g;
my $total_a;
my $total_t;
my $total_c;
my $total_n;
my $seqlen;
my $seqlen_g;
my $seqlen_a;
my $seqlen_t;
my $seqlen_c;
my $seqlen_n;
my @sequences;
my $per_gc;

my $totallength;
my $n50;
my $verbose;
my $help;

&GetOptions(
    'file:s'      => \$file,
    'total'       => \$totallength,
    'n50'         => \$n50,
    'verbose'     => \$verbose,
    'help'        => \$help,
    );

$/=">";

open (FILE, "<$file") || die "Failed to open file: '$file'\n";
while (<FILE>) {
    chop $_;

    # don't process the first (false) entry
    next if ($_ eq "");

    # ensure all chars are upper-case
    tr/[a-z]/[A-Z]/;

    # parse sequence name, and any trailing characters of header
    ($name) = (/^(\S+)(.+)\n/);

    # make string of just the sequence lines (still contains newline chars)
    $seq = substr($_,length($1));
    if (length($2) > 0) { $seq2=substr($seq,length($2)); $seq = $seq2;}

   
    # calculate stats
    $seqlen_g = $seq =~ tr/G/G/;
    $seqlen_a = $seq =~ tr/A/A/;
    $seqlen_t = $seq =~ tr/T/T/;
    $seqlen_c = $seq =~ tr/C/C/;
    $seqlen_n = $seq =~ tr/N/N/;

    # exclude padding N's from [gc] calculation
    $seqlen = ($seqlen_g + $seqlen_a + $seqlen_t + $seqlen_c);
    $per_gc = (int ( (($seqlen_g + $seqlen_c) / $seqlen) * 100)) / 100;
    $seqlen = $seqlen + $seqlen_n;

    # add sequence totals to overall file totals
    $total = $total + $seqlen;
    $total_g = $total_g + $seqlen_g;
    $total_a = $total_a + $seqlen_a;
    $total_t = $total_t + $seqlen_t;
    $total_c = $total_c + $seqlen_c;
    $total_n = $total_n + $seqlen_n;
    
    # push sequence length to array for N50 calculation (if requested)
    push (@sequences, $seqlen);

    # Output formats
    unless ($totallength) {
	print "$name\t $seqlen\t$per_gc\tG: $seqlen_g\tA: $seqlen_a\tT: $seqlen_t\tC: $seqlen_c\tN: $seqlen_n\n";
    }
}
close SHOTGUN;
$/="\n";


# Total length output 
if ( $totallength ) {
    print "Total length of seqeunces in $file is $total\n";
    $seqlen = ($total_g + $total_a + $total_t + $total_c);
    $per_gc = (int ( (($total_g + $total_c) / $seqlen) * 100)) / 100;
    print "$file\t $total\t$per_gc\tG: $total_g\tA: $total_a\tT: $total_t\tC: $total_c\tN: $total_n\n";
}

# Contig N50 value
if ( $n50 ) {
    my @sort = sort {$b <=> $a} @sequences;    
    my $seqn50;
    my $totn50;
    foreach my $val(@sort){
	$totn50+=$val;
	$seqn50++;
	if($totn50 >= $total/2){
	    print "N50 length is $totn50 and N50 value is: $val from $seqn50 sequences\n";
	    last; 
	}
    }
}

exit(0);
