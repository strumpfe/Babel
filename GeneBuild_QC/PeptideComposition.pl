#!/usr/bin/perl
# PeptideComposition.pl

# Usage: PeptideComposition.pl <filename>

use strict;
use Getopt::Long;
use carp;

my $file;
my $seq = "";
my $seq2 = "";
my $name;
my $total;
my $seqlen;
my @sequences;

my %lengths;
my $outsize;

my $no_x_containing;
my $no_x_start;
my $no_x_end;
my $no_x_long;

my $no_short_predictions_30;
my $no_short_predictions_60;
my $no_short_predictions_100;

my %comp;
my $residue;
my $total_length;
my $no_sequences;
my %sizebin;
my $binsize = 50;
my $whichbin;
my $maxbin = 20;
my $binreport;

my $minlen = 250;
my $full;
my $stat;
my $totallength;
my $report;
my $verbose;
my $help;
my $bin;

&GetOptions(
    'file:s'      => \$file,
    'minlen:i'    => \$minlen,
    'total'       => \$total,
    'full'        => \$full,
    'report'      => \$report,
    'stat'        => \$stat,
    'bin'         => \$bin,
    'binreport'   => \$binreport,
    'verbose'     => \$verbose,
    'help'        => \$help,
    );

$/=">";

open (FILE, "<$file") || die "Failed to open file: '$file'\n";
while (<FILE>) {
    chomp $_;

    # don't process the first (false) entry
    next if ($_ eq "");

    # ensure all chars are upper-case
    tr/[a-z]/[A-Z]/;

    # parse sequence name, and any trailing characters of header
#    ($name) = (/^(\S+)(.+)\n/);
    ($name) = (/^(\S+)\n/);

    # make string of just the sequence lines (still contains newline chars)
    $seq = substr($_,length($1));
    if (length($2) > 0) { $seq2=substr($seq,length($2)); $seq = $seq2;}

    # calculate stats
    $seqlen = length($seq);
    
    # bin for size
    $whichbin = int($seqlen/$binsize);
#    print " $name $seqlen $whichbin\n";
    $sizebin{$whichbin}++;
    if ($whichbin >= $maxbin) {$sizebin{$maxbin}++};

    if (($whichbin == 1) && ($binreport)) {
	print ">$name$seq";
    }


    for ($residue = 0; $residue <= $seqlen; $residue++) {
#	print substr($seq,$residue,1) . "\n" if ($verbose);
	$comp{substr($seq,$residue,1)}++;
    }
   
    # Housekeeping
    $totallength = $totallength + $seqlen;   # Increase overall length
    $no_sequences++;                         # Increment number of sequences
    $lengths{$name}=$seqlen;                 # Store sequence length in hash

    # Output formats
    if ( $full ) {
	if ( $verbose ) {
	    print "$name\t$seqlen\tAla: $comp{A} Arg: $comp{R} Asn: $comp{N} Asp: $comp{D} Cys: $comp{C} Gln: $comp{Q} Glu: $comp{E} Gly: $comp{G} His: $comp{H} Ile: $comp{I} Leu: $comp{L} Lys: $comp{K} Met: $comp{M} Phe: $comp{F} Pro: $comp{P} Ser: $comp{S} Thr: $comp{T} Trp: $comp{W} Tyr: $comp{Y} Val: $comp{V} Unkn: $comp{X}\n";
	}
	else {
	    print "$name\t$seqlen\tA: $comp{A} R: $comp{R} N: $comp{N} D: $comp{D} C: $comp{C} Q: $comp{Q} E: $comp{E} G: $comp{G} H: $comp{H} I: $comp{I} L: $comp{L} K: $comp{K} M: $comp{M} F: $comp{F} P: $comp{P} S: $comp{S} T: $comp{T} W: $comp{W} Y: $comp{Y} V: $comp{V} Unkn: $comp{X}\n";
	}
    }
    elsif ($total) {
	printf "%.50s\t$seqlen\n", $name;
    }
    elsif ($stat && $verbose) {
	printf "%.50s\t$seqlen\n", $name;
    }

    # Report 
    if ( $report ) {



    # composition based reporting
	if ( $comp{X} > 0 ) {
	    print "$name\t$comp{X} UNKN residues from $seqlen residue prediction\n";
	    $no_x_containing++;
#	    print substr($seq,1,1) . "<< >>" . substr($seq,-2,1) . "\n";
	    if ( substr($seq,1,1) eq "X" ) {
		    print "$name\tStarts with an UNKN residue, partial prediction or potential 5'-UTR issue?\n";# if ($verbose);
		    $no_x_start++;
	    }
	    if ( substr($seq,-2,1) eq "X" ) {
		    print "$name\tEnds with an UNKN residue, partial prediction?\n";# if ($verbose);
		    $no_x_end++;
	    }
        if ( $comp{X} > 5 ) {
            print "$name\tContains more than 5 aa ambiguous residues\n";# if ($verbose);
            $no_x_long++;
        }
    }

    # Length based reporting
    if ( $seqlen <= 30 ) {
        print "$name\tLength <= 30 aa\n" if ($verbose);
        $no_short_predictions_30++;
    }
    if ( $seqlen <= 60 ) {
        print "$name\tLength <= 60 aa\n" if ($verbose);
        $no_short_predictions_60++;
    }
    if ( $seqlen <= 100 ) {
        print "$name\tLength <= 100 aa\n" if ($verbose);
        $no_short_predictions_100++;
    }

	

    }

    # Clean up after yourself
    undef (%comp);


}
close SHOTGUN;
$/="\n";

if ($report) {
    print "// File: $file\n";
    print "// No. transcripts: $no_sequences\n";
    print "//\n";
    print "// Bad character reports\n";
    print "// No. entries with UNKN X's        : $no_x_containing\n";
    print "// No. entries with start UNKN X's  : $no_x_start\n";
    print "// No. entries with end UNKN X's    : $no_x_end\n";
    print "// No. entries with more than 5 X's : $no_x_long\n";
    print "//\n// Length reports\n";
    print "// No. entries less than 30 aa      : $no_short_predictions_30\n";
    print "// No. entries less than 60 aa      : $no_short_predictions_60\n";
    print "// No. entries less than 100 aa     : $no_short_predictions_100\n";
    print "\n";
}


if ($total) {
    print "// Number of entries: $no_sequences\n";
    print "// Total length is $totallength residues\n";   
    print "// Mean length is ".(int($totallength/$no_sequences)). " from $no_sequences sequences\n";
#    print "// Potential number of loci is ".(int($totallength/$minlen)). " based on mean length of $minlen aa\n";
}

if ( ($bin) &! ($binreport) ) {
    my $i;
    print "// Number of entries: $no_sequences\n";
    for ($i=0; $i<20; $i++) {
	if ($sizebin{$i} > 0) {
	    print "Bin " . ($i * $binsize) . "-" . (($i * $binsize)+$binsize) . "\t" . $sizebin{$i} . "\n";
	}
	else {
	    print "Bin " . ($i * $binsize) . "-" . (($i * $binsize)+$binsize) . "\t0\n";	    
	}
    }
    
    print "Bin    >" . ($maxbin * $binsize) .  "\t" . $sizebin{$maxbin} . "\n";	    

}



if ($stat) {

    my $meanlen = int( ($totallength/$no_sequences) + 0.5);
    my $meandiff;
    foreach my $i (sort keys %lengths) {
        $meandiff = $meandiff + ( ($lengths{$i} - $meanlen) * ($lengths{$i} - $meanlen) );
    }
    my $standdev = sqrt($meandiff / $no_sequences);
    my $lowerlimit_2s = int( ($meanlen - (2 * $standdev)) + 0.5);
    my $lowerlimit_1s = int( ($meanlen - $standdev) + 0.5);
    my $upperlimit_1s = int( ($meanlen + $standdev) + 0.5);
    my $upperlimit_2s = int( ($meanlen + (2 * $standdev)) + 0.5);


    print "// Number of entries: $no_sequences\n";
    print "// Total length is $totallength residues\n";   
    print "// Mean length is $meanlen\n";
    print "// Standard deviation is $standdev\n";
    print "//  sd size range is $lowerlimit_1s - $upperlimit_1s\n";
    print "// 2sd size range is $lowerlimit_2s - $upperlimit_2s\n";

    my $outsize_1s = 0;
    foreach my $i (sort keys %lengths) {
        # σ
        if ( $lengths{$i} < $lowerlimit_1s) { $outsize_1s++; printf "// sd small %.60s \t$lengths{$i} \t[file: $file mean: $meanlen]\n", $i; }
        if ( $lengths{$i} > $upperlimit_1s) { $outsize_1s++; printf "// sd big   %.60s \t$lengths{$i} \t[file: $file mean: $meanlen]\n", $i; }
    }
    print "// No of entries outsize  sd: $outsize_1s\n";
   
    my $outsize_2s = 0;
    foreach my $i (sort keys %lengths) {
        # 2σ
        if ( $lengths{$i} < $lowerlimit_2s) { $outsize_2s++; printf "// 2sd small %.60s \t$lengths{$i} \t[file: $file mean: $meanlen]\n", $i; }
        if ( $lengths{$i} > $upperlimit_2s) { $outsize_2s++; printf "// 2sd big   %.60s \t$lengths{$i} \t[file: $file mean: $meanlen]\n", $i; }
    }
    print "// No of entries outsize 2sd: $outsize_2s\n";

}



exit(0);
