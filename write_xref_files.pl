#!/usr/bin/perl
use strict;
use Getopt::Long;
use IO::Handle;


my $verbose;
my $verboser;
my $help;
my $summary;
my $full;
my $all,
my $organism;        # Code (4 letters) to designate which species to process
my @organism;
my $filetype;        # Xref file type to generate
# filetype options
my $symbol;          # write symbol file
my $synonym;         # write synonym file
my $citation;        # write citation file
my $description;     # write description file
# other
my $release;         # Release version
my $updated;         # updated?
my $protcod;         # protein-coding only
my $rna;             # ncRNA only
my $write;           # write output to file
my $terminal;        # write output to terminal
# Synonyms
my %synonym_data;
# Summary
my %summary_data;
my %sum_citation;
my %sum_families;
my %update_author;
my %update_species;
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
    "release=s"     => \$release,
    "updated=s"     => \$updated,
    # run mode
    "summary"       => \$summary,
    "full"          => \$full,
    "all"           => \$all,
    # Xref file
    #"filetype"      => \$filetype,
    "symbol"        => \$symbol,
    "synonym"       => \$synonym,
    "citation"      => \$citation,
    "description"   => \$description,
    "protcod"       => \$protcod,
    "rna"           => \$rna,
    # Output
    "write"         => \$write,
    "terminal"      => \$terminal,
    );
#---------------------------------------------------------#

# pod documentation for help
if ( $help ) {
  exec ('perldoc',$0);
}

# Check that we have a release version
unless ( $release || $updated) {
  print "# $0 " . gmtime( time()) ."\n\n";
  print "No release declared. Aborting.\n";
  print " '-release' option to set the VectorBase release version\n\n";
  exit(0);
}
unless ( ($release) =~ (/^VB\-\d+\-\d+$/) || $updated) {
  print "# $0 " . gmtime( time()) ."\n\n";
  print "Release string looks weird. It should be of the format VB-2001-01.\n";
  print " '-release' is set to '$release'\n";
  exit(0);
}

# Join items in array @species to make the binomial name
$organism = join ' ',@organism;

# Check that we have declared a species
unless ( $organism || $summary || $all) {
  print "# $0 " . gmtime( time()) ."\n\n";
  print "No organism declared. Aborting.\n";
  print " '-organism' option to set organism to process\n\n";
  exit(0);
}

# Check that we have declared a filetype
unless ( $symbol || $synonym || $citation || $description || $full || $summary || $updated ) {
  print "# $0 " . gmtime( time()) ."\n\n";
  print "No filetype declared. Aborting.\n";
  print "Choose one of the following:\n '-symbol', '-synonym', '-citation', '-description'\n\n";
  exit(0);
}
# Assign filetype
$filetype = "symbol"      if ( $symbol );
$filetype = "synonym"     if ( $synonym  );
$filetype = "citation"    if ( $citation );
$filetype = "description" if ( $description );

#---------------------------------------------------------#

##
## Main options, multi-species reporting or single species processing
##

# '-summary' mode for summary of all species
if ( $summary && !$updated) {
	# open filehandle for output file
	open (OUTPUT, "> ${release}_xref_data.txt") if ( $write );
	# print table header
	print "+---------------------+                           +--------------------------------+--------------------------------+\n";
	print "| Release: ".$release." |                           | Protein-coding CDS             | non-coding RNA                 |\n";
	print "+---------------------+-----+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+\n";
	print "| Organism                  | Rows     | Active   | CDS      | Symbol   | Describe | ncRNA    | Symbol   | Describe | Synonym  | Families | Uniq Fam | Citation | Papers   |\n";
	print "+---------------------------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+\n";

	# print organism data rows
	if  ( $organism ) {
		&parse_spreadsheet($organism);
	}
	else {
		while (<DATA>) {
    		chomp;
    		last if /^__END__$/;         # stop when __END__ is encountered
    		next if /^\s*(#.*|\s*)$/;    # skip blank lines or comment lines
			%summary_data = "";
			&parse_spreadsheet($_);
			%sum_citation = "";
			%sum_families = "";
		}
	}
	# print table footer
	print "+---------------------------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+----------+\n";
	# close filehandle for output file
	close OUTPUT if ( $write );
}
# '-full' mode exports all data types for a single organism
elsif ( $full || $all) {
#	my @runs = ('symbol','synoym','citation','description');


	# print organism data rows
	if  ( $organism ) {
		&export_spreadsheet($organism);
	}
	else {
		while (<DATA>) {
    		chomp;
    		last if /^__END__$/;         # stop when __END__ is encountered
    		next if /^\s*(#.*|\s*)$/;    # skip blank lines or comment lines

			undef %synonym_data;         # clear any legacy synonym data

			&export_spreadsheet($_);
		}
	}
}
# -updated mode, display rows modified since timestamp $updated
elsif ( $updated ) {

	# print organism data rows
	if  ( $organism ) {
		&parse_spreadsheet($organism);
	}
	else {
		while (<DATA>) {
    		chomp;
    		last if /^__END__$/;         # stop when __END__ is encountered
    		next if /^\s*(#.*|\s*)$/;    # skip blank lines or comment lines
#    		%summary_data = "";
			&parse_spreadsheet($_);
		}
	}
	print "\n";
	print "// Some overview stats for gene metadata additions since $updated\n//\n";
	print "// No. of updated symbols:   $summary_data{Symbol}\n";
	print "// No. of updated synonyms:  $summary_data{Synonym}\n";
	print "// No. of updated citations: $summary_data{Citation}\n";
	print "// Species:\n";
	foreach my $i (sort {$update_species{$b} <=> $update_species{$a} or "$a" cmp "$b"} keys %update_species) {printf ("// => %3s\t$i\n", $update_species{$i})}
	print "// Submitters:\n";
	foreach my $i (sort {$update_author{$b} <=> $update_author{$a} or "$a" cmp "$b"}   keys %update_author)  {printf ("// => %3s\t$i\n", $update_author{$i})}
	print "\n";
}
# default exports a single data type for a single organism
else {
	# open filehandle for output file
	my $filename = "${release}_".join('-',@organism)."_xref_${filetype}.txt";
	print "// Exporting $filetype for $organism\n";
	print "// Writing output to $filename\n";
	open (OUTPUT, "> $filename");
	# process organism data
	&parse_spreadsheet($organism);
	# close filehandle for output file
	close OUTPUT;
}

exit(0);

#---------------------------------------------------------#


sub parse_spreadsheet {

	my ($organism) = @_;
	my $key = $orgkeys{$organism2code{$organism}};
	my $curlcommand = "https://docs.google.com/spreadsheets/d/".$key."/export?exportFormat=tsv";
	my @f;

	# Get gene metadata via curl and process
	#---------------------------------------------------------#

	open (FILE, "curl --silent $curlcommand |") or print "// WARNING failed to open spreadsheet\n";
	while (<FILE>) {

		# ignore header lines
		next if (/VB_ID/);

		chomp;
		# write row to output file if -summary -write option for debugging purposes
		if ( ( $summary ) && ( $write ) ) {
			print OUTPUT "$_\n";
		}


		@f = split/\t/;

		# Symbol file {get data and write file}
		#---------------------------------------------------------#
		if ( $symbol) {
			# discard loci which aren't live (ACTIVE)
			next unless ( $f[2] eq "ACTIVE" );
			# process only those loci with a symbol
 			next unless ( $f[7] ne "");
 			# Protein-coding only?
 			if ( $protcod ) { next unless $f[3] eq "protein_coding";}
 			# ncRNA only?
 			if ( $rna ) { next unless $f[3] =~ /RNA/ ;}

			# write data to OUTPUT filehandle or STDOUT
 			if ( $terminal ) {
				print "$f[0]\t$f[1]\tgene\t$f[7]\t$f[10]\n";
			}
			else {
				print OUTPUT "$f[0]\t$f[1]\tgene\t$f[7]\t$f[10]\n";
			}
		}

		# Synonym file {get data}
		#---------------------------------------------------------#
		if ( $synonym) {
			# discard annotations which aren't synonyms
			next unless ( $f[2] eq "SYNONYM" );
			print "// Locus: $f[1] Symbol: $f[7] Synonym: $f[8]\n" if $verbose;
 			# store the synoym in hash
 			$synonym_data{$f[1]} = $synonym_data{$f[1]} . "$f[8], ";
 			print "// Assign $f[8] as synonym of $f[1] ($f[7])\n" if $verbose;
		}

		# Citation file {get data and write}
		#---------------------------------------------------------#
		if ( $citation) {
			# discard annotations which aren't synonyms
			next unless ( $f[2] eq "CITATION" );
			print "// Locus: $f[1] Symbol: $f[7] PubMed: $f[12] Supplied_by: $f[11]\n" if $verbose;
			# write data to OUTPUT filehandle or STDOUT
 			if ( $terminal ) {
		 		print "$f[12]\t$f[1]\tgene\t$f[7]\tSupplied by $f[11]\n";
			}
			else {
		 		print OUTPUT "$f[12]\t$f[1]\tgene\t$f[7]\tSupplied by $f[11]\n";
			}
	 	}

		# Description file {get data and write}
		#---------------------------------------------------------#
		if ( $description) {
			# discard annotations which aren't synonyms
			next unless ( $f[2] eq "ACTIVE" );
			# process only those loci with a symbol
	 		next if ( $f[7] ne "");
	 		next unless ( $f[10] ne "");

			print "// Locus: $f[1] Symbol: $f[7] PubMed: $f[12] Supplied_by: $f[11]\n" if $verbose;

			# write data to OUTPUT filehandle or STDOUT
 			if ( $terminal ) {
		 		print "$f[1]\t$f[10]\n";
			}
			else {
		 		print OUTPUT "$f[1]\t$f[10]\n";
			}
	 	}

		# Summary file {get data}
		#---------------------------------------------------------#
		if ( $summary && !$updated ) {
			$summary_data{Entries}++;
			$summary_data{Active}++         if ( $f[2] eq "ACTIVE" );

			$summary_data{CDS}++            if ( $f[3] eq "protein_coding") && ( $f[2] eq "ACTIVE" );   # Add requirement for row to have ACTIVE status
			$summary_data{ncRNA}++          if ( $f[3] =~ /RNA/)  && ( $f[2] eq "ACTIVE" );             # Add requirement for row to have ACTIVE status

			$summary_data{Symbol}++         if ( $f[7]  ne "" );
			$summary_data{Symbol_CDS}++     if ( ($f[7] ne "") && ($f[3] eq "protein_coding") );
			$summary_data{Symbol_ncRNA}++   if ( ($f[7] ne "") && ($f[3] =~ /RNA/) );

			$summary_data{Synonym}++        if ( $f[2] eq "SYNONYM" );

			$summary_data{Description}++    if ( $f[10]  ne "" );
			$summary_data{Desc_CDS}++       if ( ($f[10] ne "") && ($f[3] eq "protein_coding") );
			$summary_data{Desc_ncRNA}++     if ( ($f[10] ne "") && ($f[3] =~ /RNA/) );

			$summary_data{Family}++         if ( $f[9] ne "" );

			if ( $f[9] ne "" ) {
				$sum_families{$f[9]}++;
			}

			if ( $f[2] eq "CITATION" ) {
				$summary_data{Citation}++;
				$sum_citation{$f[12]}++;
 			}
		}
		# Updated summary
		#---------------------------------------------------------#
		if ( $updated ) {
			# $f[15]
			my $date = $f[15];
			$date =~ s/-//g;
			if ( $date > $updated) {
				printf ("[$date >> $updated] %16s %-25s %-25s [AC: $f[1] GN: %-10s DE: %-64s]\n", $f[0], $f[11], $f[13], $f[7], substr($f[10],0,64) ) if ( $verbose );

				# Warnings for incomplete curations..
				# Submitter
				if ( $f[11] eq "" ) { printf ("WARNING: $f[0] %-25s has no submitter name\n", $f[13]); }
				# Species
				if ( $f[13] eq "" ) { printf ("WARNING: $f[0] (Author: %-25s) has no species name\n", $f[11]); }

				# Count update type
				if ( $f[2] eq "ACTIVE")      {$summary_data{Symbol}++;}
				elsif ( $f[2] eq "SYNONYM")  {$summary_data{Synonym}++;}
				elsif ( $f[2] eq "CITATION") {$summary_data{Citation}++;}
				# Push update species & author
				$update_species{$f[13]}++;
				$update_author{$f[11]}++;
			}
		}
	}
	close FILE;

	##
	##  Tidy up writing synonym and summary files
	##

	# Synonym file {write file}
	#---------------------------------------------------------#
	if ( $synonym ) {
		my $names;
		foreach my $i (sort keys %synonym_data) {
			$names = $synonym_data{$i};
			chop $names; chop $names;
			# write data to OUTPUT filehandle or STDOUT
 			if ( $terminal ) {
				print "$i\t$names\n";
			}
			else {
				print OUTPUT "$i\t$names\n";
			}
		}
	}
	#_ end of synonym reporting


	# Summary file {write file}
	#---------------------------------------------------------#
	if ( $summary && !$updated) {
		my $no_papers   = scalar (keys %sum_citation) ; if ($no_papers) {$no_papers--;} # 2015-01-14 added -1 modifier
		my $no_families = scalar (keys %sum_families) ; if ($no_families) {$no_families--;} # 2015-01-14 added -1 modifier

		# if verbose write long list, else table row
		if ( $verbose ) {
			print "Organism:                  $organism\n";
			# spreadsheet details
			print "No. of rows:               $summary_data{Entries}\n";
			print "No. of active rows:        $summary_data{Active}\n";
			# protein-coding CDS
			print "No. of active CDS:         $summary_data{CDS}\n";
			print "No. of Symbols CDS:        $summary_data{Symbol_CDS}\n";
			print "No. of Description CDS:    $summary_data{Desc_CDS}\n";
			# ncRNAs
			print "No. of active ncRNA:       $summary_data{ncRNA}\n";
			print "No. of Symbols ncRNA:      $summary_data{Symbol_ncRNA}\n";
			print "No. of Description ncRNA:  $summary_data{Desc_ncRNA}\n";
			# synonyms
			print "No. of Synonyms:           $summary_data{Synonym}\n";
			# families
			print "No. of Family members:     $summary_data{Family}\n";
			print "No. of Families:           $no_families\n";

			foreach my $i (sort keys %sum_families){print "=> $sum_families{$i}\t$i\n";}

			# citations
			print "No. of Citations:          $summary_data{Citation}\n";
			print "No. of Papers:             $no_papers\n";
			print "\n";
		}
		else {
			printf ("| %-25s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s |\n", $organism,$summary_data{Entries},$summary_data{Active},$summary_data{CDS},$summary_data{Symbol_CDS},$summary_data{Desc_CDS},$summary_data{ncRNA},$summary_data{Symbol_ncRNA},$summary_data{Desc_ncRNA},$summary_data{Synonym},$summary_data{Family},$no_families,$summary_data{Citation},$no_papers);
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

	if ( $symbol || $full ) {
		my $filename = "${release}_${organism4filename}_xref_symbol.txt";
		print "// Exporting symbol for $organism ($organism4filename)\n";
		print "// Writing output to $filename\n";
		open (OUTPUT, "> $filename");
		&parse_spreadsheet($organism);
		close OUTPUT;
	}

	if ( $synonym || $full ) {
		my $filename = "${release}_${organism4filename}_xref_synonym.txt";
		print "// Exporting synonym for $organism ($organism4filename)\n";
		print "// Writing output to $filename\n";
		open (OUTPUT, "> $filename");
		&parse_spreadsheet($organism);
		close OUTPUT;
	}

	if ( $citation || $full ) {
		my $filename = "${release}_${organism4filename}_xref_citation.txt";
		print "// Exporting citation for $organism ($organism4filename)\n";
		print "// Writing output to $filename\n";
		open (OUTPUT, "> $filename");
		&parse_spreadsheet($organism);
		close OUTPUT;
	}

	if ( $description || $full ) {
		my $filename = "${release}_${organism4filename}_xref_description.txt";
		print "// Exporting description for $organism ($organism4filename)\n";
		print "// Writing output to $filename\n";
		open (OUTPUT, "> $filename");
		&parse_spreadsheet($organism);
		close OUTPUT;
	}
}
#_ end of export_spreadsheet subroutine

## Fields from spreadsheet
# Curation_ID
# Stable_ID
# Status
# Biotype
# Model
# From
# To
# Symbol
# Synonym
# GeneFamily
# Description
# Submitter
# Value
# Species
# Submitter comments
# Last update

#-----------------------------------------------------------------------------------------------#

# Data hashes
#---------------------------------------------------------#
sub getdata {
	%organism2code = (
		'Aedes aegypti'             => 'AAEL',
                'Aedes albopictus'          => 'AALF',
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
                'Sarcoptes scabiei'         => 'SSCA',
                'Stomoxys calcitrans'       => 'SCAU',
		);

	%code2organism = reverse %organism2code;

	%orgkeys = (
		'AAEL' => '1D6Ax2HJyWcnukIkPu2Zio7hZ_ajQiZbvzfb1DfINLwo',
                'AALF' => '1T8nVia4ATqIZKNi8FbMR9ikvEb9JDxGTppgkTUfLw3Q',
		'AALB' => '1IkszmWkEzJ5lTHw6ugOSGFV_0soguBIoTPOh7-CE51o',
		'AARA' => '1a5OtCuDUGxbrg_0rMcbDQHUHXtwbzO_g9vTm_K896UY',
		'AATE' => '1oxUeyJVjzrq4oe5LlVqHzEGB1Hd2_7_75sCGuWxEkT0',
		'ACHR' => '1QqhBJB4vlArR9jwaQvYSbxVGnljPRuBWQvUqLmibr7U',
		'ACOM' => '1Qc7GVQO4OXIPBuiIzITDPzjbKffsnobBg7s2NmeCr7Q',
		'ACUA' => '1Y9FjfL5bYHHSgDx26sRLLCMNVIlqqyT6T_ny9n0XKwQ',
		'ADAC' => '1dcOTVBAAtfV1Vksekr1BOZ0NKxUl8zqWiz12ULmQ2hQ',
		'ADIR' => '1haqU3Z0pctX3G1sMhVb8lk0GDtEvGObmODEcFwqMJwE',
		'AEPI' => '1US2FOLzgg1lXVkcmvFItp-wrQJ63FAISNfv63cMFdvc',
		'AFAR' => '1haKWlsjUJrKq9wsHnbMD72NwSlFsw1J_aeYsDBf-NpA',
		'AFUN' => '1u8cMMKSvwAn7Hm6NW6iELDHbxUc38K41buAySCtqdRo',
	 	'AGAP' => '1YbJ0JbFWnXhnDqGy63NQJIbkvjjg3Uuq1MoL5xF2_Oc',
	 	'AMAC' => '10m31MlBGhlSysA8vsrKxjz_2PcSuCUt4FdWVJ8myHP8',
	 	'AMEC' => '1_p3WGzyH6Q0XW6c-POXqvGrALPUxGILqTR3qmovQXEg',
	 	'AMEM' => '17qjOAt7UwqZkL3QOvYnNBUzYAeKS5pVTbc8wIycUjGE',
	 	'AMIN' => '1JhOKWHOAcAVkCzjkT_yDK8FNslBu5wiiadIle_jGu40',
		'AQUA' => '1jWDRiITV3_XP14CKfNogbDzADQCf3nctZzJKKMNmqUg',
		'ASIN' => '1NQL-9mmykZ0fWls6405mhCKSxIVfE2kvxaFRogb1tXY',
		'ASTE' => '1ZVMIQ_tr3zJBMuk-KdUcEXSIqVOjL8jANc_SKqQAlkU',
		'BGLB' => '1PPSMgx2-z6CPa4Cit5Kh1IUCa9i7Z6Rpnv4H5ks38UA',
		'CPIJ' => '1_6KHioJdjHYroDooP48FavcclKs1D2vqxPZvJzcTq_U',
		'GAUT' => '1VWQ8g_JXkz8L4XLvJJ5AVGwjhr4NP5am25-DUtJcjc0',
		'GBRI' => '1n8LzP0VOEUbgsFtkr2rX6FJNUEfr9CGCfKikrI0WdVo',
		'GFUI' => '131u1wnBlJaIlbLbmih2zu2q1lCZhPFK3oNfo9LD5gY0',
		'GMOY' => '1OytW4SzLunXCphzrX_mz6CwIsfn-zIzrt8eYmy6nbBk',
		'GPAI' => '1qh1tGOHr21wKHSq6r6JwNNLE-x1K4x7kDxSEXULpaus',
		'GPAP' => '1hRRFAQC6YJTxJAPDe5SuF0-I7se4jgicfPkcXPVU86E',
		'ISCW' => '10hoIT-Zhsk1byeB5nissrbPeKXGKO5ODt0rPwTzk8fY',
		'LLOJ' => '1COtYQkfRVbbAzKKNDuJgiQuebivOZbmf_z2RJJSoCJE',
		'MDOA' => '1ZNXGtXI_6L_qJjz8iMt02vNCFVdlD4Z0iisSqM2hw1A',
		'PHUM' => '1WeBgVThFg8mW8-xrqtfMSWmSEGk9iwKi-8RWZDTucEg',
		'PPAI' => '1qH5OdX08naKPWK21Eow7_l3K7CkoTy8MgtLX835rFPM',
		'RPRC' => '1ZDgS4pOzhKIwtS8h5Ho-Yxot1iZGu_9RWi-EWrJdWA4',
                'SSCA' => '1pUcRMAx2Ua21j87owwrowN9S5JMcZvAdrmfaY2piKEo',
                'SCAU' => '1xIVpy1SH4r7eyx4QB9gtQm091AEp2UBbKrWVG6N_qsc',
		);
}
#_ end of subroutine data

#---------------------------------------------------------#


__DATA__
Aedes aegypti
Aedes albopictus
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
Sarcoptes scabiei
Stomoxys calcitrans
__END__


=pod

=head1 NAME - write_xref_files.pl

=head2 USAGE

This script writes tsv files for use in VectorBase xref pipelines.

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

=item All species,

C<-release VB-2015-02 -symbol -write>, export single datatype for all species

C<-release VB-2015-02 -full -write -summary>, export all data types for all species

=item Single species,

C<-release VB-2015-02 -organism Anopheles gambiae -symbol -write>, export single datatype for single species

C<-release  VB-2015-02 -organism Anopheles gambiae -full -write>, export single datatype for single species

=back

=back

=head2 ARGUMENTS

B<write_xref_files.pl> arguments:

=head3 File_types

=over 4

=item -full, write all files (symbol|synonym|citation|description)

=item -symbol, write gene symbol file

=item -synonym, write gene synonym file

=item -citation, write PubMed citation file

=item -description, write gene description paper

=item -protcod, process only protein-coding loci (useful for symbol/description)

=back

=head3 Selecting_organism

=over 4

=item -organism, Binomial name of organism to use as filter, e.g. Aedes aegypti

=back

=head3 Misc

=over 4

=item -terminal, write output to STDOUT rather than file

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
