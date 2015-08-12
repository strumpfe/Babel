#!/usr/bin/perl
use strict;
use Getopt::Long;
use IO::Handle;


my $verbose;
my $verboser;
my $help;
my $summary;
my $organism;        # Code (4 letters) to designate which species to process
my @organism;
# other
my $write;           # write output to file
my $terminal;        # write output to terminal
# Summary
my %summary_data;
my %species_count;
# get organism and spreadsheet data
my %organism2code;
my %code2organism;
my %orgkeys;
&getdata;

#---------------------------------------------------------#
GetOptions (
    # Misc
    "verbose"       => \$verbose,
    "help"          => \$help,
    # Species
    "organism=s{2}" => \@organism,
    # run mode
    "summary"       => \$summary,
    # Output
    "write"         => \$write,
    "terminal"      => \$terminal,
    );
#---------------------------------------------------------#

# pod documentation for help
if ( $help ) {
  exec ('perldoc',$0);
}

# Join items in array @species to make the binomial name
$organism = join ' ',@organism;

# Check that we have declared a species
unless ( $organism || $summary ) {
  print "# $0 " . gmtime( time()) ."\n\n";
  print "No organism declared. Aborting.\n";
  print " '-organism' option to set organism to process\n\n";
  exit(0);
}

#---------------------------------------------------------#

##
## Main options, multi-species reporting or single species processing
##

# '-summary' mode for summary of all species
if ( $summary ) {
	# print table header
	print "+---------------------------+----------+\n";
	print "| Organism                  | Loci     |\n";
	print "+---------------------------+----------+\n";

	# print organism data rows
        &parse_spreadsheet;

	# print table footer
        print "+---------------------------+----------+\n";
	# close filehandle for output file
}
elsif ( $write || $terminal ) {

	# print organism data rows
	if  ( $organism ) {
		&export_spreadsheet($organism);
	}
}

exit(0);

#---------------------------------------------------------#


sub parse_spreadsheet {

        my ($organism) = @_;
	my $key = "1y9HGZAvsQfniqJNAL75ZO5VJFKMcYo9gQpzCTErRaoY";
	my $curlcommand = "https://docs.google.com/spreadsheets/d/".$key."/export?exportFormat=tsv";
	my @f;

	# Get gene metadata via curl and process
	#---------------------------------------------------------#

	open (FILE, "curl --silent $curlcommand |") or print "// WARNING failed to open spreadsheet\n";
	while (<FILE>) {

		# ignore header lines
		next if (/Timestamp/);

#                print if $verbose;

		chomp;
		@f = split/\t/;

		# Summary file {get data}
		#---------------------------------------------------------#
		if ( $summary ) {
			$summary_data{Entries}++;
			$summary_data{Active}++         if ( $f[2] eq "ACTIVE" );

                        # Organism $f[1]
                        $species_count{$f[1]}++;        # increment species count
		}
                elsif ( $f[1] eq $organism ) {    # print loci name for organism
                        print OUTPUT "$f[2]\n" if $write;
                        print "$f[2]\n" if $terminal;
                }


	}
	close FILE;

# Timestamp	Species	Locus name	Contact email	Rationale for deletion

	# Summary file {write file}
	#---------------------------------------------------------#
	if ( $summary ) {
                foreach my $i (sort keys %species_count) {
			printf ("| %-25s | %8s |\n", $i,$species_count{$i});
		}
	}
	#_ end of summary reporting
}
#_ end of parse_spreadsheet subroutine

sub export_spreadsheet {

	my ($organism) = @_;
	my $organism4filename = $organism;
	$organism4filename =~ s/ /_/;

#	print "// Exporting \n";

        my $filename = "${organism4filename}_loci_2_delete.txt";

	if ($write) {
                print "// Exporting loci for $organism ($organism4filename)\n";
	        print "// Writing output to $filename\n";
	        open (OUTPUT, "> $filename");
                }
       &parse_spreadsheet($organism);

	close OUTPUT if ($write);

}
#_ end of export_spreadsheet subroutine

#-----------------------------------------------------------------------------------------------#

# Data hashes
#---------------------------------------------------------#
sub getdata {
	%organism2code = (
		'Aedes aegypti'             => 'AAEL',
		'Anopheles albimanus'       => 'AALB',
		'Anopheles arabiensis'      => 'AARA',
		'Anopheles atroparvus'      => 'AATE',
		'Anopheles christyi'        => 'ACHR',
		'Anopheles coluzzii'        => 'ACOM',
		'Anopheles culicifacies'    => 'ACUA',
		'Anopheles darlingi'        => 'ADAC',
		'Anopheles dirus'           => 'ADIR',
		'Anopheles epiroticus'      => 'AEPI',
		'Anopheles farauti'         => 'AFAR',
		'Anopheles funestus'        => 'AFUN',
		'Anopheles gambiae'         => 'AGAP',
		'Anopheles maculatus'       => 'AMAC',
		'Anopheles melas'           => 'AMEC',
		'Anopheles merus'           => 'AMEM',
		'Anopheles minimus'         => 'AMIN',
		'Anopheles quadriannulatus' => 'AQUA',
		'Anopheles sinensis'        => 'ASIN',
		'Anopheles stephensi'       => 'ASTE',
		'Biomphalaria glabrata'     => 'BGLB',
		'Culex quinquefasciatus'    => 'CPIJ',
		'Glossina austeni'          => 'GAUT',
		'Glossina brevipalpis'      => 'GBRI',
		'Glossina fuscipes'         => 'GFUI',
		'Glossina morsitans'        => 'GMOY',
		'Glossina pallidipes'       => 'GPAI',
		'Glossina palpalis'         => 'GPAP',
		'Ixodes scapularis'         => 'ISCW',
		'Lutzomyia longipalpis'     => 'LLOJ',
		'Musca domestica'           => 'MDOA',
		'Pediculus humanus'         => 'PHUM',
		'Phlebotomus papatasi'      => 'PPAI',
		'Rhodnius prolixus'         => 'RPRC',
		);

	%code2organism = reverse %organism2code;

}
#_ end of subroutine data

#---------------------------------------------------------#


__DATA__
Aedes aegypti
Anopheles albimanus
Anopheles arabiensis
Anopheles atroparvus
Anopheles christyi
Anopheles coluzzii
Anopheles culicifacies
Anopheles darlingi
Anopheles dirus
Anopheles epiroticus
Anopheles farauti
Anopheles funestus
Anopheles gambiae
Anopheles maculatus
Anopheles melas
Anopheles merus
Anopheles minimus
Anopheles quadriannulatus
Anopheles sinensis
Anopheles stephensi
Biomphalaria glabrata
Culex quinquefasciatus
Glossina austeni
Glossina brevipalpis
Glossina fuscipes
Glossina morsitans
Glossina pallidipes
Glossina palpalis
Ixodes scapularis
Lutzomyia longipalpis
Musca domestica
Pediculus humanus
Phlebotomus papatasi
Rhodnius prolixus
__END__


=pod

=head1 NAME - get_delete_candidates.pl

=head2 USAGE

This script writes list of candidate loci for deletion as part of a patch build.

=head3 Examples..

=over 2

=item -find updates for all species since a defined timepoint,

Use the I<-updated> option to find annotations that have been added or modifed since a particular date. If you add the I<-verbose> flag to get a list of the changes.

=over 2

=item All species,

C<-updated 20150101 -summary>

=item Single species,

C<-updated 20150101 -organism Anopheles gambiae>

=back

=item -write data to flatfile (for external use)

=over 2

=item Summary of all species

C<-summary>, delete candidate counts for all species

=item Single species,

C<-organism Anopheles gambiae -write>, export loci for single species

=back

=back

=head2 ARGUMENTS

B<get_delete_candidates.pl> arguments:


=head3 Selecting_organism

=over 4

=item -organism, Binomial name of organism to use as filter, e.g. Aedes aegypti

=back

=head3 Misc

=over 4

=item -verbose, verbose mode toggle on extra command line output.

=item -help, these help pages.

=back

=head3 Summary

=over 4

=item -summary, Summary of data broken down by organism, project and data type.

=item -write, Write all spreadsheet rows to output file for debugging purpose

=back

=head1 AUTHOR

Dan Lawson (lawson@ebi.ac.uk)

=cut
