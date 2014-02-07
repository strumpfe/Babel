#!/usr/bin/perl

use strict;
use Getopt::Long;
use carp;

# variable for command-line options
my $gfffile;
my $coords;
my $feature;
my $slice;
my $slicefile;
my $window;
my $windowsize = 100000;
my $verbose;
my $help;

#---------------------------------------------------------#
&GetOptions(
    'gff:s'        => \$gfffile,
    'slicefile:s'  => \$slicefile,
    'coords:s'     => \$coords,
    'feature:s'    => \$feature,
    'slice'        => \$slice,
    'window'       => \$window,
    'windowsize:s' => \$windowsize,
    'verbose'      => \$verbose,
    'help'         => \$help,
    );
#---------------------------------------------------------#

# pod documentation for help
if ( $help ) {
  exec ('perldoc',$0);
}

# check that either -slice or -window mode has been selected
unless ( $slice || $window ) { 
	print "// You need to run in either -slice or -window mode\n";
	exit;
}

# check that a gff file have been provided
unless ( $gfffile ) { 
	print "// You haven't provided a gff file (-gff)\n";
	exit;
}

# check that either -coords or -slicefile have been provided
unless ( $coords || $slicefile ) { 
	print "// You haven't provided a slice (-coords) or list of slices (-slicefile)\n";
	exit;
}

# check format of coordinate slice if defined
if ( $coords && $slice ) { 
	unless ( $coords =~ /^\S+\:\d+\-\d+/ ) {
		print "// The slice (-coords) provided doesn't contain enough information region:start-end\n";
		exit;
	}
}

if ( $coords && $window ) { 
	if ( $coords =~ /^\S+\:\d+\-\d+/ ) {
		print "// The slice (-coords) provided contain too much information, only need a sequence name\n";
		exit;
	}
}

# check that -feature has been provided
unless ( $feature ) { 
	print "// You haven't provided a feature type (-feature)\n";
	exit;
}


# variables for main script
my $slice_count;     # single count for slice option
my %window_count;    # hash count for window option
my $seq;
my $start;
my $end;

my @f;
my $chunk;
my @slice2process;

# Queue slices for processing
if ($slicefile) {
	open (FILEofSLICES, "< $slicefile") or die "Failed to open $slicefile : $!";
	while (<FILEofSLICES>){
		chomp;	
		push(@slice2process,$_);
		print STDOUT "// Submitting $_ for processing\n" if ($verbose);
	}
	close FILEofSLICES;
}
else {
	push(@slice2process,$coords);
	print STDOUT "// Submitting $coords for processing\n" if ($verbose);
}

# run jobs
for (my $i = 0; $i < scalar(@slice2process); $i++) {
	# Convert 
	if ($slice) {
		($seq,$start,$end) = $slice2process[$i] =~ (/^(\S+):(\d+)-(\d+)/);	
		print "// Looking for genes within sequence $seq from $start to $end\n" if ($verbose);
		&findfeature_from_gff;
	}
	elsif ($window) {
		$seq = $slice2process[$i];
		print "// Looking for genes within sequence $seq with a window size of $windowsize\n" if ($verbose);			
		&findfeature_from_gff($seq);
	}
}

##
## Subroutines
##

sub findfeature_from_gff {

	print "// DEBUG info  Looking at $seq for $feature with coordinates [$start - $end]\n" if ($slice);
	print "// DEBUG info  Looking at $seq for $feature with window size $windowsize\n"    if ($window);



	# Loop through GFF3 file
	open (FILE, "< $gfffile") or die "Failed to open $gfffile : $!";
	while (<FILE>) {
	#	chomp;
		@f = split/\t/;

		# discard unless correct feature
		next unless ( $f[2] eq $feature);

		# discard unless correct seq_region
		next unless ( $f[0] eq $seq );

		if ($slice) {
			# discard if before start coord
			next if ( $f[3] < $start );
			# discard if before start coord
			last if ( $f[3] > $end );
			print if ($verbose);
			$slice_count++;
		}
		elsif ($window){
			$chunk = int($f[3]/$windowsize);
			print "//Feature starts at $f[3] == chunk $chunk\n" if ($verbose);
			$window_count{$chunk}++;
		}	
	}
	close FILE;

	if ($slice) {
		printf "$seq [%8s - %8s]\t%5s\n", $start,$end,$slice_count;
	}
	elsif ($window){
		foreach my $i (sort {$a<=>$b} keys %window_count) {
			printf "$seq [%8s - %8s]\t%5s\n", ($i*$windowsize),(($i+1)*$windowsize),$window_count{$i};
		}
	}
}

=pod
 
=head1 NAME - find_GFF_annotation.pl

=head2 USAGE 

This script counts features from a GFF3 file either for a slice or sequencial windows of user defined length (default is 100,000 bp).

=head2 ARGUMENTS

B<find_GFFF_annotation.pl> arguments:

=head3 Mandatory

=over 4

=item -gff, Path of the gff3 file from which data will be retreived

=item -feature, Feature type to count e.g. gene, mRNA, exon

=item -slice, Calculate feature counts in slice mode (must define seqregion and start/stop coordinates).

=item -window, Calculate feature counts in sliding window mode (must define seqregion only).

=back

=head3 Data

=over 4

=item -coords, Sequence region name (and start/end coordinates if -slice) for a single region/scaffold

=item -slicefile, <optional> list of slices or scaffolds to be processed, formatted appropriately for either -slice or -window

=back

=head3 Misc

=over 4

=item -windowsize, <optional> Define the size of the window (in bp), default is 100,000 bp

=back

=head3 Misc

=over 4

=item -verbose, Verbose mode toggle on extra command line output
 
=item -help, these help pages

=back

=head1 AUTHOR

Dan Lawson (lawson@ebi.ac.uk)

=cut



