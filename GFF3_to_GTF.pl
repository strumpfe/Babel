#!/usr/bin/perl
# GFF3_to_GTF.pl

# Usage: GFF3_to_GTF.pl -options <filename>

use strict;
use Getopt::Long;
use carp;

my $file;           # File name of input (GFF3 format)
my $file_root;      # Root of file name (i.e. name minus .gff3)
my $outfile;        # Output file name (root.gtf)
my $errfile;        # Error file name (root.err)

my $gene;
my $transcript;
my $utr;
my $coding_start;
my $coding_end;
my $ATG_start;
my $ATG_end;
my $STOP_start;
my $STOP_end;
my $parsed_CDS_exon;
my %exon_no;

my $verbose;
my $verboser;
my $help;


#---------------------------------------------------------#
&GetOptions(
    'file:s'    =>  \$file,
    'verbose'   =>  \$verbose,
    'verboser'  =>  \$verboser,
    'help'      =>  \$help,
    );
#---------------------------------------------------------#



##
## Generate outfile name (based on input GFF3)
##

print "// Reading GFF3 file: $file\n" if ($verbose);
$file_root = $file;
$file_root =~ s/.gff3//g;
$outfile = $file_root . ".gtf";
$errfile = $file_root . ".err";
print "// Writing GTF  file: $outfile\n" if ($verbose);
print "// Writing errors to: $errfile.$$\n" if ($verbose);


##
## Calculate the number of exons and coding exons per transcript
##

my ($no_exons,$no_cds_exons) = &get_data($file);
my %no_exons                 = %$no_exons;
my %no_cds_exons             = %$no_cds_exons;

#exit(0);

##
## Parse gff file
##

open (ERROR,  "> $errfile.$$") || die "Failed to open error file: '$errfile.$$'\n";
open (OUTPUT, "> $outfile")    || die "Failed to open output file: '$outfile'\n";
open (FILE,   "< $file")       || die "Failed to open input file: '$file'\n";
while (<FILE>) {

    my @f = split /\t/;

    # gene
    if ( $f[2] eq "gene" ) {

        # Test whether the GFF3 has biotype declaration on the gene object
        unless (/biotype/) {
            print "# $0 " . gmtime( time()) ."\n\n";
            print "GFF3 gene line has no biotype declaration. Aborting.\n";
            print "Biotypes where added to the GFF3 for the VB-2013-12 release.\n";
            print "Re-dump GFF3 from the core database and try again.\n\n";
            exit(0);
        }

        ($gene) = (/ID=(\S+);biotype/);
        print OUTPUT "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\";\n";
    }
    # mRNA
    elsif ( grep /^$f[2]$/i, ("mRNA","tRNA","miRNA","snRNA","snoRNA","misc_RNA","rRNA","pseudogene") ) {

        ($transcript) = (/ID=(\S+);Parent/);

        print "// New transcript $transcript [CDS exons $no_cds_exons{$transcript} of $no_exons{$transcript}]\n" if ($verbose);
	    print ERROR "$gene\t$transcript\t$no_exons{$transcript}\t$no_cds_exons{$transcript}\n";
        print OUTPUT "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$transcript\";\n";

        undef $coding_start;
        undef $coding_end;
        undef $utr;
        $parsed_CDS_exon = 0;
    }
    # CDS
    elsif ( $f[2] eq "CDS" ) {

        ##############
        # +ve strand #
        ##############
        if ( $f[6] eq "+" ) {
            $parsed_CDS_exon++;

            # First exon
            if ( !$coding_start ) {
                $ATG_start = $f[3];
                $ATG_end   = $f[3] + 2;
                print OUTPUT "$f[0]\t$f[1]\tstart_codon\t$ATG_start\t$ATG_end\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$transcript\";\n";
                print OUTPUT "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$transcript\"; exon_number \"$parsed_CDS_exon of $no_exons{$transcript}\";\n";
		        $coding_start = 1;
            } #_ end of first exon block
            # Final exon
            elsif ( !$coding_end ) {
                if ( $parsed_CDS_exon == $no_cds_exons{$transcript} ) {
                    $STOP_start = $f[4] - 2;
                    $STOP_end   = $f[4];
                    print OUTPUT "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$transcript\"; exon_number \"$parsed_CDS_exon of $no_exons{$transcript}\";\n";
                    print OUTPUT "$f[0]\t$f[1]\tstop_codon\t$STOP_start\t$STOP_end\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$transcript\";\n";
                    $coding_end = 1;
                }
                else {
                    print OUTPUT "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$transcript\"; exon_number \"$parsed_CDS_exon of $no_exons{$transcript}\";\n";
                }
            } #_ end of final exon block

        } #_ end of +ve strand parsing

        ################
        ## -ve strand ##
        ################

        if ( $f[6] eq  "-" ) {
            $parsed_CDS_exon++;

            # First exon
            if ( !$coding_start ) {
                $ATG_start = $f[4] - 2;
                $ATG_end   = $f[4];
                print OUTPUT "$f[0]\t$f[1]\tstart_codon\t$ATG_start\t$ATG_end\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$transcript\";\n";
                print OUTPUT "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$transcript\"; exon_number \"$parsed_CDS_exon of $no_exons{$transcript}\";\n";
                $coding_start = 1;
            } #_ end of first exon block
            # Final exon
            elsif  ( !$coding_end ) {
                if ( $parsed_CDS_exon == $no_cds_exons{$transcript} ) {
                    $STOP_start = $f[3];
                    $STOP_end   = $f[3] + 2;
                    print OUTPUT "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$transcript\"; exon_number \"$parsed_CDS_exon of $no_exons{$transcript}\";\n";
                    print OUTPUT "$f[0]\t$f[1]\tstop_codon\t$STOP_start\t$STOP_end\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$transcript\";\n";
                    $coding_end = 1;
                }
                else {
                    print OUTPUT "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$transcript\"; exon_number \"$parsed_CDS_exon of $no_exons{$transcript}\";\n";
                }
            } #_ end of final exon block
        } #_ end of -ve strand parsing
    } #_ end of CDS parsing

    ##########
    ## Exon ##
    ##########
    elsif ( $f[2] eq "exon" ) {
        $exon_no{$transcript}++;
        print OUTPUT "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$transcript\"; exon_number \"$exon_no{$transcript} of $no_exons{$transcript}\";\n";
    } #_ end of exon parsing

    # five_prime_utr
    elsif ( $f[2] eq "five_prime_utr" ) {
        print OUTPUT "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$transcript\";\n";
        if ( $utr > 0 ) { $parsed_CDS_exon++; }
        $utr = 1;
    }
    # three_prime_utr  (really redundant with above)
    elsif ( $f[2] eq "three_prime_utr" ) {
        print OUTPUT "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tgene_id \"$gene\"; transcript_id \"$transcript\";\n";
    }
    else {
        print;
    }

}
close (FILE);
close (ERROR);
close (OUTPUT);


exit(0);

#----------------------------------------------------------------------------------------------------------------------------------------------------------#

##
## Subroutines
##


##
## Get exon & CDS count per transcript from GFF file
##

sub get_data {

    my $infile = shift;             # string  : path of the input GFF file
    my @f;                          # array   : GFF fields
    my %no_exons;                   # hash    : number of exons keyed by transcript name
    my %no_cds_exons;               # hash    : number of coding exons keyed by transcript name

    open (GFF, "<$infile");
    while (<GFF>) {

        chomp;
    	@f = split (/\t/);

        # parse exon information
        if ( ($f[2] eq "exon") && ($f[8] =~ /Parent=(\S+.+)\;/) ) {
	    $no_exons{$1}++;
        $no_cds_exons{$1} = 0;
	    print "// Increment exon count for $1 to $no_exons{$1}\n" if ($verboser);
	    }

        # parse CDS information
        if ( ($f[2] eq "CDS") && ($f[8] =~ /Parent=(\S+.+)\;/) ) {
	    $no_cds_exons{$1}++;
	    print "// Set CDS count for $1 to $no_exons{$1} (i.e. this is the last exon number)\n" if ($verboser);
        }
    }
    close GFF;

    return (\%no_exons,\%no_cds_exons);

} #_ end of get_data sub


#-----------------------------------------------------------------------------------------------#


__DATA__
gene
mRNA
exon
CDS
five_prime_utr
three_prime_utr
__END__



=pod
 
=head1 NAME - GFF3_to_GTF.pl

=head2 USAGE 

This script writes GTF format files based on GFF3 input for the VectorBase project.

=head2 ARGUMENTS

B<GFF3_to_GTF.pl> arguments:

=head3 Datafiles

=over 4

=item -file, GFF3 file to convert to GTF format 

=back

=head3 File types

=over 4

=item -bam, Select BAM (bam) file configuration for ini file

=item -bigwig, Select bigwig (bw) file configuration for ini file
 
=back

=head3 Selecting tracks

=over 4

=item -species, Binomial name of organism to use as filter, e.g. Aedes aegypti

=back

=head3 Misc

=over 4

=item -verbose, Verbose mode toggle on extra command line output
 
=item -verboser, (Even more) Verbose mode toggle to see full debug infromation
 
=item -help, these help pages

=back

=head1 AUTHOR

Dan Lawson (lawson@ebi.ac.uk)

=cut




