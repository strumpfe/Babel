#!/usr/bin/perl
# Simple_genebuild_stats.pl

# Usage: Simple_genebuild_stats.pl <gff3_filename>

use strict;
use Getopt::Long;
use carp;

my $file;                # GFF3 filename
my @f;
my %features;            # hash of feture type counts
my %featurelength;       # cumulative length of all features
my %trans_exon_length;   # Length of transcript represented as exons
my %trans_cds_length;    # Length of transcript represented as CDS
my %trans_intron_length; # Length of transcript represented as introns

my $featurespan;

my $id;
my $parent;
my %noexon;
my %nocds;
my %nointron;

my $lastparent;
my $lastscaffold;
my $laststart;
my $laststop;
my $laststrand;
my $gene;
my $intronID;
my $intronsize;

my $bin;
my $binsize = 50;    # Size of bins to count length frequencies
my %bincount;

my $stats;
my $exons;
my $singleexons;
my $introns;
my $report;
my $verbose;
my $help;

#---------------------------------------------------------#
&GetOptions(
    'file:s'      => \$file,
    'report'      => \$report,
    'stats'       => \$stats,
    'exons'       => \$exons,
    'singleexons' => \$singleexons,
    'introns'     => \$introns,
    'verbose'     => \$verbose,
    'help'        => \$help,
    );
#---------------------------------------------------------#


#
# Open GFF3 file
# 

open (FILE, "<$file") || die "Failed to open file: '$file'\n";
while (<FILE>) {
    chop $_;

    # don't process the comment lines
    next if (/^\#/);

    # Split the file into tab-delimited fields
    @f = split/\t/;

    # Housekeeping
    $features{$f[2]}++;

    # Length statistics
    if ( ($f[2] eq 'exon') or ($f[2] eq 'CDS') or ($f[2] eq "five_prime_utr") or ($f[2] eq "three_prime_utr")) {
            $featurespan = $f[4] - $f[3] + 1;
            $featurelength{$f[2]} = $featurelength{$f[2]} + $featurespan;
    }

    # exon - transcript relationship
    if ($f[2] eq 'exon') {
            ($parent) = $f[8] =~ (/Parent=(\S+);/);
            $noexon{$parent}++;
            $trans_exon_length{$parent} = $trans_exon_length{$parent} + $f[4] - $f[3] + 1;
            print "// $parent\t$noexon{$parent} $trans_exon_length{$parent}\n" if ($verbose);
    }

    # CDS - transcript relationship
    if ($f[2] eq 'CDS') {
            ($parent) = $f[8] =~ (/Parent=(\S+);/);
            $nocds{$parent}++;
            $trans_cds_length{$parent} = $trans_cds_length{$parent} + $f[4] - $f[3] + 1;
    }


    # Handling introns...
    if ($f[2] eq "exon") {
        ($parent) = $f[8] =~ (/Parent=(\S+);/);
        ($gene) = $parent =~ (/(\S+)\-/);
        if ( $parent eq $lastparent ) {
            print "// Matching parent for exon from $parent => intron\n" if ($verbose);
            print "// End of previous exon : $lastscaffold ($laststrand) $laststop\n" if ($verbose);
            print "// Current exon         : $f[0] ($f[6]) $f[3]\n" if ($verbose);
            print "// Intron coordinates   : " . ($laststop + 1) . " - " . ($f[3] - 1) . "\n" if ($verbose);
            print "// Intron length        : " . ($f[3] - $laststop + 1)  . "\n" if ($verbose);
            
            $nointron{$parent}++;
            $intronID = $parent."-".$nointron{$parent};
            if ( $f[6] eq '+') { $intronsize = $f[3] - $laststop - 1; }
            if ( $f[6] eq '-') { $intronsize = $laststart - $f[4] - 1; }
            if ($intronsize < 0) {
                    print "// Reverse coordinates from  " . ($f[3] - $laststop + 1) . " to " . ($laststop - $f[3]+ 1) . "\n" if ($verbose);
                    $intronsize = $laststop - $f[3] + 1;  
            }
            
            # Intron length bins
            $trans_intron_length{$intronID} = $intronsize;
            $bin = int($trans_intron_length{$intronID}/$binsize);
            $bincount{$bin}++;

            # print the GFF3 intron line
            if ( $f[6] eq '+' ) {
                print "$f[0]\t$f[1]\tintron\t".($laststop + 1)."\t".($f[3] - 1)."\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$parent\"; length \"$intronsize\"\n";
            }
            elsif ( $f[6] eq '-' ) {
                print "$f[0]\t$f[1]\tintron\t".($f[4] + 1)."\t".($laststart - 1)."\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$parent\"; length \"$intronsize\"\n";
            }

            $features{intron}++;

            print "\n" if ($verbose);
        } 

        # reset pointers based on current exon
        $lastparent   = $parent;
        $lastscaffold = $f[0];
        $laststrand   = $f[6];
        $laststart    = $f[3];
        $laststop     = $f[4];
    }



}
close (FILE);

#exit(0);

&feature_report if ($stats);
&exon_count     if ($exons);
&single_exons   if ($singleexons);
&intron_count   if ($introns);

exit(0);

#----------------------------------------------------------------------------------------------------------------------------------------------------------#

##
## Subroutines
##

sub feature_report {

    print  "+------------+--------+\n";
    print  "| Feature    |  Count |\n";
    print  "+------------+--------+\n";
    printf "| Gene       | %6s |\n", $features{gene};
    printf "| mRNA       | %6s |\n", $features{mRNA};
    printf "| Pseudogene | %6s |\n", $features{pseudogene};
    print  "+------------+--------+----------+\n";
    printf "| exon       | %6s | %8s |\n", $features{exon},  $featurelength{exon};
    printf "| CDS        | %6s | %8s |\n", $features{CDS},  $featurelength{CDS};
    print  "+------------+--------+----------+\n";
    printf "| UTR 5'     | %6s | %8s |\n", $features{five_prime_utr},  $featurelength{five_prime_utr};
    printf "| UTR 3'     | %6s | %8s |\n", $features{three_prime_utr}, $featurelength{three_prime_utr};
    print  "+------------+--------+----------+\n";

    my $total_ncRNA = $features{tRNA} + $features{rRNA} + $features{miRNA} + $features{snRNA} + $features{snoRNA} + $features{misc_RNA} ;

    printf  "| ncRNA      | %6s |\n", $total_ncRNA;
    print  "+------------+--------+\n";
    printf "| tRNA       | %6s |\n", $features{tRNA};
    printf "| rRNA       | %6s |\n", $features{rRNA};
    printf "| miRNA      | %6s |\n", $features{miRNA};
    printf "| snRNA      | %6s |\n", $features{snRNA};
    printf "| snoRNA     | %6s |\n", $features{snoRNA};
    printf "| misc_RNA   | %6s |\n", $features{misc_RNA};
    print  "+------------+--------+\n";

}


sub exon_count {

    my %exonfreq;
    my $exonmax = 1;

    foreach my $i (sort keys %noexon) {
        print "$i\t$noexon{$i}\n" if ($verbose);
        $exonfreq{$noexon{$i}}++;
        if ($noexon{$i} > $exonmax) { $exonmax = $noexon{$i}; }
    }

    print "+-------+-------+--------+\n";
    print "| Exons | Count | % mRNA |\n";
    print "+-------+-------+--------+\n";
    for (my $i = 1; $i <= $exonmax; $i++ ) {
        if ( $exonfreq{$i} eq "" ) { $exonfreq{$i} = 0; }
        printf "| %5s | %5s | %6.2f |\n", $i, $exonfreq{$i}, &percent($exonfreq{$i},$features{mRNA});
    }
    print "+-------+-------+--------+\n";

    print "No. of exons:  $features{exon}\n";
    print "No. of mRNAs:  $features{mRNA}\n";
    printf "Mean no exons: %.2f\n", ($features{exon}/$features{mRNA}); 

}


sub single_exons {

    my $i;

    print "+---------------+--------+\n";
    print "| Locus         | Length |\n";
    print "+---------------+--------+\n";
     foreach $i (sort { $trans_exon_length{$b} <=> $trans_exon_length{$a} } keys %trans_exon_length) {
        if ($noexon{$i} == 1) {
           printf "| $i | %6s |\n", $trans_exon_length{$i};
        }
    }
    print "+---------------+--------+\n";
}


sub intron_count {

    my %intronfreq;
    my $intronmax = 1;

    foreach my $i (sort keys %nointron) {
        print "$i\t$nointron{$i}\n" if ($verbose);
        $intronfreq{$nointron{$i}}++;
        if ($nointron{$i} > $intronmax) { $intronmax = $nointron{$i}; }
    }

    print "+---------+-------+--------+\n";
    print "| Introns | Count | % mRNA |\n";
    print "+---------+-------+--------+\n";
    for (my $i = 1; $i <= $intronmax; $i++ ) {
        if ( $intronfreq{$i} eq "" ) { $intronfreq{$i} = 0; }
        printf "| %7s | %5s | %6.2f |\n", $i, $intronfreq{$i}, &percent($intronfreq{$i},$features{mRNA});
    }
    print "+---------+-------+--------+\n";

    print "No. of introns:  $features{intron}\n";
    print "No. of mRNAs:  $features{mRNA}\n";
    printf "Mean no introns: %.2f\n", ($features{intron}/$features{mRNA}); 

    # Intron length frequency distribution
    print "+---------------+-----------+\n";
    print "| Intron size   | Frequency |\n";
    print "+---------------+-----------+\n";
    
    my $binrange;
    foreach my $i (sort { $a <=> $b } keys %bincount) {
        $binrange = (($i*$binsize)-$binsize) .'-'. ($i*$binsize); 
  #    printf "| %13s | %6s |\n", $binrange,$bincount{$i};
        printf "%12s |", $bincount{$i};
    }
    print "+---------------+-----------+\n";


#    print "+------------------+--------+\n";
#    print "| Intron           | Length |\n";
#    print "+------------------+--------+\n";
#     foreach my $i (sort { $trans_intron_length{$a} <=> $trans_intron_length{$b} } keys %trans_intron_length) {
     #   if ($nointron{$i} == 1) {
#           printf "| %16s | %6s |\n", $i, $trans_intron_length{$i};
     #   }
#    }
#    print "+------------------+--------+\n";

}


sub percent {
    my ($i,$tot) = @_;
    my $p = int ((($i / $tot)*1000) + 0.5 )/10;
    return $p;
}


