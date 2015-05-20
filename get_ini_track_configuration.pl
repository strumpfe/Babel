#!/usr/bin/perl

## get_track_configuration.pl
## lawson@ebi.ac.uk
## 2013-12-01 v0.2

use strict;
use Getopt::Long;
use IO::Handle;
use Data::Dumper;
use Net::FTP;

# Hash name
our %SRA;

# Options
my $verbose;             # Verbose messages
my $moreverbose;         # More verbosity
my $debug;               # Write out data parsing debug information (very verbose)
my $help;                # Help documentation via POD
my $local;               # Read from local hash .dat file
my $fusion;              # Read from Google Fusion table
my $bam;                 # Write bam configuration
my $bigwig;              # Write bigWig configuration
my $species;             # Filter by species
my @species;             # Temp array to capture genus/species from command line. Concatenated into the $species var
my $summary;             # Track summary statistics
my $checkURLs;           # Validate source_URLs

# Google Fusion table ID & Key
my $fusiontable = "1tHQgMvCjvbZ36jg3Kgl32Y9eiVYtGfE8S_sYXls";
my $key         = $ENV{APIKEY};                                   # Get API key from environment variable or declare it on command line
my $dir         = $ENV{"HOME"} . "/DATA_FILES";                   # Define local directory for data files
my $data_file   = "VectorBase_SRA_info.dat";                      # Define output dat file for local copy of hash

#---------------------------------------------------------#
GetOptions (
    # Misc
    "verbose"       => \$verbose,
    "verboser"      => \$moreverbose,
    "help"          => \$help,
    # Data
    "local"         => \$local,
    "fusion"        => \$fusion,
    # File type
    "bam"           => \$bam,
    "bigwig"        => \$bigwig,
    # Search
    "species=s{2}"  => \@species,
    # General
    "summary"       => \$summary,
    "check"         => \$checkURLs,
    "key"           => \$key,
    );
#---------------------------------------------------------#

# Join items in array @species to make the binomial name
$species = join ' ',@species;

# Verbose option summary of script run
if ($moreverbose) {
#    print "# $0 " . gmtime( time()) ."\n";
    print "- Summary of data structure contents\n" if ($summary);
    print "- Track type: BAM\n" if ($bam);
    print "- Track type: bigWig\n" if ($bigwig);
    print "- Species   : $species\n";
}

# pod documentation for help
if ( $help ) {
  exec ('perldoc',$0);
}

# Check that we have a Google API key
unless ($key) {
  print "# $0 " . gmtime( time()) ."\n\n";
  print "No API key declared. Aborting.\n";
  print "Either set this as an environment variable (APIKEY)\n";
  print "Or set one using the following option:\n";
  print " '-key' option to set the Google Fusion API key\n\n";
  exit(0);
}

## Read SRA track information
## either from Fusion table or local file

# retrieve data from Google Fusion table
if ($fusion) {
  print "// Using data from Fusion table :: $fusiontable\n\n" if ( $moreverbose );
  &get_data_from_fusion;
}
# retrieve data from local hash structure
elsif ($local) {
  # Do we want to have abiliity to work with alternative hash structure, define via file option
  # retrieve the hash from the file.
  print "// Using data from hash :: $dir/$data_file\n\n" if ( $moreverbose );
  &get_data_from_local_file;
}
else {
  print "# $0 " . gmtime( time()) ."\n\n";
  print "No data source declared. Aborting.\n";
  print "Use one of the following options:\n";
  print " '-local' option to read from local hash structure\n";
  print " '-fusion' to read from Google Fusion table\n\n";
  print " '-help' for full documentation\n\n";
  exit(0);
}

unless ( ($species) or ($summary) ) {
  print "# $0 " . gmtime( time()) ."\n\n";
  print "No species declared. Aborting.\n";
  print "Use one of the following options:\n";
  print " '-species' option to declare which species to parse\n";
  exit(0);
}

# Check that species assignment is valid
if ( $species ) {
  my $found;

  while (<DATA>) {
    chomp;
    last if /^__END__$/;         # stop when __END__ is encountered
    next if /^\s*(#.*|\s*)$/;    # skip blank lines or comment lines
    print "// Attempting to match '$_' to $species\n" if ($moreverbose);
    if ($_ eq $species) { $found = 1; print "// Matched species $species\n" if ($moreverbose); }
  }
  unless ( $found > 0 ) {        # species name not matched
    print "# $0 " . gmtime( time()) ."\n\n";
    print "Couldn't match species name '$species'\n";
    print "Add the species binomial name to the __DATA__ secion at the foot of this script and re-run\n\n";
    exit(0);}
}


##
## Command option selected subroutines
##

&track_summary_stats if ($summary);

&write_bam_configuration($species)    if ($bam);
&write_bigwig_configuration($species) if ($bigwig);

exit(0);

#----------------------------------------------------------------------------------------------------------------------------------------------------------#

##
## Subroutines
##

##
## Write out BAM configuration
##

sub write_bam_configuration {

  my ($species) = @_;
  my $items; 

  # Repprt and check before proceeding
  print "// Write out bam configuration for the ini file for $species\n" if ($verbose);
  for my $i (sort keys %SRA) {                                         # $i = track name
    if ( $SRA{$i}->{species} eq $species ) { $items++; }
  }
  print "// Found $items tracks for $species\n" if ($verbose);
  if ( $items == 0 ) {                                                 # Species name is correct but no data available
    print "# $0 " . gmtime( time()) ."\n\n";
    print "No track data found for $species'\n";
    print "Add track data to the Google Fusion table and re-run\n\n";
    exit(0);
   }

  my @f;
  my $value;
  my %source_name;
  my %description;
  my %source_url;
  my $inter;

  # write header lines
  print "#############\n";
  print "# BAM CONFIG\n";
  print "#############\n";
  print "\n";
  print "[ENSEMBL_INTERNAL_BAM_SOURCES]\n";

  # write list of tracks
  foreach my $i (sort keys %SRA) {                     # $i = track name
    next unless ( $SRA{$i}->{species} eq $species );   # Only process those tracks from the selected species
    next unless ( $SRA{$i}->{source_type} eq "bam" );     # Only process those tracks with type "bam"
     print "// WARNING: Note sure how we got here unless $SRA{$i}->{species} == $species\n" if ($moreverbose);
    print "$i = rnaseq_align\n";
  }
  print "\n";

  # write track configuration
  foreach my $i (sort keys %SRA) {                     # $i = track name
    next unless ( $SRA{$i}->{species} eq $species );   # Only process those tracks from the selected species
    next unless ( $SRA{$i}->{source_type} eq "bam" );     # Only process those tracks with type "bam"
    print "[$i]\n";
    print "source_name    = " . $SRA{$i}->{source_name} . "\n";
    print "description    = $SRA{$i}->{description}\n";
    print "source_url     = $SRA{$i}->{source_url_bam}\n";
    print "source_type    = bam\n";
    print "display        = off\n\n";
  }
  print "\n";

} #_ end of write_bam_configuration subroutine

#-----------------------------------------------------------------------------------------------------------------------------------------------------------#

##
## Write out bigWig configuration
##

sub write_bigwig_configuration {

  my ($species) = @_;
  my $items; 

  # Repprt and check before proceeding
  print "// Write out bigwig configuration for the ini file for $species\n" if ($verbose);
  for my $i (sort keys %SRA) {       # $i = track name
    if ( $SRA{$i}->{species} eq $species ) { $items++; }
  }
  print "// Found $items tracks for $species\n" if ($verbose);
  if ( $items == 0 ) {       # Species name is correct but no data available
    print "# $0 " . gmtime( time()) ."\n\n";
    print "No track data found for $species'\n";
    print "Add track data to the Google Fusion table and re-run\n\n";
    exit(0);
   }

  my @f;
  my $value;
  my %source_name;
  my %description;
  my %source_url;
  my $inter;


  # write header lines
  print "#################\n";
  print "# bigWig CONFIG #\n";
  print "#################\n";
  print "\n";
  print "[ENSEMBL_INTERNAL_BIGWIG_SOURCES]\n";

  # write list of tracks
  foreach my $i (sort keys %SRA) {                           # $i = track name
    next unless ( $SRA{$i}->{species} eq $species );         # Only process those tracks from the selected species
    next unless ( $SRA{$i}->{source_type} eq "bigwig" );     # Only process those tracks with type "bigwig"
    print "$i = rnaseq_align\n";
  }
  print "\n";

  # write track configuration  (2nd iteration over whole list, should push to an array and read from that)
  foreach my $i (sort keys %SRA) {                           # $i = track name
    next unless ( $SRA{$i}->{species} eq $species );         # Only process those tracks from the selected species
    next unless ( $SRA{$i}->{source_type} eq "bigwig" );     # Only process those tracks with type "bigwig"
    print "[$i]\n";
    print "source_name    = $SRA{$i}->{source_name}\n";
    print "caption        = $SRA{$i}->{caption}\n";
    print "description    = $SRA{$i}->{description}\n";
    print "source_url     = $SRA{$i}->{source_url_bigwig}\n";
    print "source_type    = bigWig\n";
    print "display        = off\n";
    print "\n";
  }
  print "\n";

} #_ end of write_bigwig_configuration subroutine

#-----------------------------------------------------------------------------------------------------------------------------------------------------------#

##
## gat_data_from_fusion
##

sub get_data_from_fusion{

  my $track;
  my $description;
  my $value;
  my @f;
  my $start;

#  print "curl -s https://www.googleapis.com/fusiontables/v1/query? -d \"sql=SELECT%20Track_name,Species,Assembly,Source_name,Caption,Description,Source_URL_bam,Source_URL_bigwig,Source_type,Display,SRA_project,Release%20FROM%20$fusiontable&key=$key\" \n";

  open (FUSIONTABLE, "curl -s https://www.googleapis.com/fusiontables/v1/query? -d \"sql=SELECT%20Track_name,Species,Assembly,Source_name,Caption,Description,Source_URL_bam,Source_URL_bigwig,Source_type,Display,SRA_project,Release%20FROM%20$fusiontable&key=$key\" |");
  while (<FUSIONTABLE>) {
    if ( /rows/ ) {$start = 1;print "<<ON>>\n" if ($debug); push @f,"NULL"; next;}
    next unless ( $start );
    chomp;
    if ( $moreverbose ) { print "// Line: '$_'\n" if ($debug); }

    if    (/\[/)        { next;}
    elsif (/\]/)        { print "<<OFF>>\n\n" if ($debug);
     if ( $f[1] ne "" ) {
            $track = $f[1]; 
            print "// Writing data for track '$track'\n" if ( $debug );

            # Check for invalid chars in track name
            unless ( $track =~ m/^[a-zA-Z0-9\s\_]+$/) {
                print "// Track name contains non-standard characters: '$track'\n";
                exit(0);
            }

            # Check that track name is unique
            foreach my $i (sort keys %SRA){
                if ( $i eq $track ) { print "// Duplicate Track name found: $track. This will need resolving before continuing.\n\n"; exit(0);}
            }

            # Correct for URL encofing of fields with HTML tags  { Description }
            $description = $f[6]; $description =~ s/\\u003c/</g; $description =~ s/\\u003e/>/g;

            $SRA{$track}->{track}              = $track;
            $SRA{$track}->{species}            = $f[2];
            $SRA{$track}->{assembly}           = $f[3];
            $SRA{$track}->{source_name}        = $f[4];
            $SRA{$track}->{caption}            = $f[5];
            $SRA{$track}->{description}        = $description;
            $SRA{$track}->{source_url_bam}     = $f[7];
            $SRA{$track}->{source_url_bigwig}  = $f[8];
            $SRA{$track}->{source_type}        = $f[9];
            $SRA{$track}->{display}            = $f[10];
            $SRA{$track}->{sra_project}        = $f[11];
            $SRA{$track}->{released}           = $f[12];

            print "// $track,$f[2],$f[3],$f[4]\n\n<<ON>>\n" if ($debug);

            }
      @f=""; 
      next;
      }
    elsif (/(\".+\")/,) { $value = $1; $value =~ s/"//g; push @f,$value;  next; }
    elsif (/(\"\")/,)   { push @f,"NULL";}

  }
  close FUSIONTABLE;

  print "// Wrote " . scalar %SRA . " tracks to hash\n" if ($verbose);

  # write hash to dat file on local disk (for future reference)
  if ($local) {
    print "// Writing local copy of hash (as you declared the -local option)\n" if ($verbose);
    &write_data_to_local_file;
  }
}

#-----------------------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------------------#

##
## get_data_from_file
##

sub get_data_from_local_file {
    open (FH, "< $dir/$data_file") or die "$dir/$data_file : $!";
    undef $/;
    my $data = <FH>;
    eval $data;
    die if $@;
    $/ = "\n";
    close FH;
}


sub write_data_to_local_file {
    my $str =Data::Dumper->Dump( [\%SRA],[qw/*SRA/]);
    open (FH, "> $dir/$data_file") or die "Can't open $dir/$data_file as output .dat\n";
    print FH $str;
    close FH;
    print "// Wrote hash %SRA to '$dir/$data_file'\n\n" if ($verbose);
}

#-----------------------------------------------------------------------------------------------#

##
## Summary stats of SRA hash structure
##

sub track_summary_stats {

  # count variables
  my $no_tracks;                 # No of tracks
  my $no_released;               # Count of released (live) tracks
  my %species_count;             # Count of species names
  my %assembly_count;            # Count of assembly names
  my %project_count;             # Count of project by accession
  my %project_by_species;        # Count of project by species
  my %project_count_by_species;  # Count of project by species
  my $species_project;           # Concatenation for tracking species/project
  my $species;                   # temp var to capture project species designations
  my $project;                   # temp var to capture project species designations
  my %track_by_type;             # Count of tracks by type
  my $ftp;                       # FTP server object
  my $lookfile;

  # read through hash
  foreach my $i (sort keys %SRA) {
    print STDOUT "$SRA{$i}->{track}\t$SRA{$i}->{species}\n" if ( $verbose );

    # Tracks (sanity check)
    $no_tracks++;
    # Species
    $species_count{$SRA{$i}->{species}}++;
    # Assembly
    $assembly_count{$SRA{$i}->{assembly}}++;
    # Release
    $no_released++ if ( {$SRA{$i}->{release}} ne "");
    # Projects
    $project_count{$SRA{$i}->{sra_project}}++;
    $species_project = $SRA{$i}->{species} . '-' . $SRA{$i}->{sra_project};
    $project_by_species{$species_project}++;
    # Track types
    $track_by_type{$SRA{$i}->{source_type}}++;
  } #_ end of loop through hash

  # count some 
  my $no_species    = scalar (keys %species_count);
  my $no_assemblies = scalar (keys %assembly_count);
  my $no_projects   = scalar (keys %project_count);
  my $no_types      = scalar (keys %track_by_type);

  ##
  ## report
  ##
  print "\n";
  print "+---------------------------+\n";
  print "| SRA RNAseq tracks         |\n";
  print "+---------------------------+--------+\n";
  printf ("| No. of tracks             | %6s |\n", $no_tracks);
  printf ("| No. of tracks released    | %6s |\n", $no_released);
  print  "+---------------------------+--------+\n";
  printf ("| Species                   | %6s |\n", $no_species);
  print  "+---------------------------+--------+\n";
  foreach my $i (sort {$species_count{$b} <=> $species_count{$a} or "$a" cmp "$b"} keys %species_count) {
    next if ( $i eq "Species");
    printf ("| %-25s | %6s |\n", $i, $species_count{$i});
  }
  print  "+---------------------------+--------+\n";
  printf ("| Assemblies                | %6s |\n", $no_assemblies);
  print  "+---------------------------+--------+\n";
  foreach my $i (sort {$assembly_count{$b} <=> $assembly_count{$a} or "$a" cmp "$b"} keys %assembly_count) {
    next if ( $i eq "Assembly");
    printf ("| %-25s | %6s |\n", $i, $assembly_count{$i});
  }
  print   "+---------------------------+--------+\n";
  printf ("| Projects                  | %6s | (Projects can contain data from multiple species)\n", $no_projects);
  print   "+---------------------------+--------+\n";

  foreach my $i (sort {$project_by_species{$b} <=> $project_by_species{$a} or "$a" cmp "$b"} keys %project_by_species) {
    print   "| full list of species-project |\n" if ( $moreverbose );
    printf ("| %-25s | %6s |\n", $i, $project_by_species{$i}) if ( $moreverbose );
    ($species,$project) = split(/-/,$i);
    $project_count_by_species{$species}++; 
  }
  print   "| summary list of species-project |\n" if ( $moreverbose );
  foreach my $i (sort {$project_count_by_species{$b} <=> $project_count_by_species{$a} or "$a" cmp "$b"} keys %project_count_by_species) {
    printf ("| %-25s | %6s |\n", $i, $project_count_by_species{$i});
  }
  print "+---------------------------+--------+\n";
  printf ("| Data file types           | %6s |\n", $no_types);
  print   "+---------------------------+--------+\n";
  foreach my $i (sort {$track_by_type{$b} <=> $track_by_type{$a} or "$a" cmp "$b"} keys %track_by_type) {
    printf ("| %-25s | %6s |\n", $i, $track_by_type{$i});
  }
  print "+---------------------------+--------+\n\n";
 
  # Check existance of source_url links
  if ( $checkURLs ) {
    print "// Checking validity of source_url links\n// Files which don't exist or can't be accessed at the moment:\n";
    $ftp = Net::FTP->new("ftp.vectorbase.org", Debug => 0)   or die "Cannot connect to some.host.name: $@";
    $ftp->login("anonymous",'-anonymous@')                   or die "Cannot login ", $ftp->message;
   
    $, = "\n";
    foreach my $i (sort keys %SRA) {
      print "// Checking $SRA{$i}->{track}\t$SRA{$i}->{species}\t$SRA{$i}->{source_url}\n" if ($verbose);
      $lookfile = $SRA{$i}->{source_url}; $lookfile =~ s/ftp:\/\/ftp.vectorbase.org//g;
      my @files = $ftp->ls($lookfile);
      unless (scalar @files > 0) {
          print "// WARNING : File $SRA{$i}->{source_url} [$i]\n";
      }
    }
    print "\n";
  } #_ end of check link routine

  exit(0);   # Exit program here, don't want to return and parse files
}

#-----------------------------------------------------------------------------------------------#

__DATA__
Aedes aegypti
Anopheles albimanus
Anopheles arabiensis
Anopheles atroparvus
Anopheles christyi
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
Ixodes scapularis
Lutzomyia longipalpis
Musca domestica
Pediculus humanus
Phlebotomus papatasi
__END__
 
=pod
 
=head1 NAME - get_ini_track_configuraton.pl

=head2 USAGE 

This script writes the ini configuration section for BAM or bigWig files for the VectorBase project.

=head2 ARGUMENTS

B<get_ini_track_configuraton.pl> arguments:

=head3 Datafiles

=over 4

=item -local, Read data from local hash structure (as defined by variables $dir/$data_file) 

=item -fusion, Read data from Google fusion table

=item -fusion -local, Read data from Google fusion table and write hash to local disk for future reference (offline mode)

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
 
=item -summary, Summary of data broken down by organism, project and data type.

=item -check, Check that source_url paths are valid (i.e. that the files are present on the FTP site)

=item -help, these help pages

=back

=head1 AUTHOR

Dan Lawson (lawson@ebi.ac.uk)

=cut

