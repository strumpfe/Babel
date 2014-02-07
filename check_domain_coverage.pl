#!/usr/bin/perl


$verbose = 1;

# Get database names
while (<DATA>) {
    chomp;
    last if /^__END__$/;         # stop when __END__ is encountered
    next if /^\s*(#.*|\s*)$/;    # skip blank lines or comment lines
    $database{$_} = 1;
}


##
## Loop through each database and return the meta table
##

foreach $db (sort keys %database) {

	print "// Get data from database: $db\n" if ($verbose);

	open (GENES, "mysql -hmysql-eg-vectorbase -P4367 -uensro $db -e 'select count(*) from gene where biotype = \"protein_coding\";' |");
	while (<GENES>) {
		chomp;
		if ( /^\|\s+(\d+)\s+\|/ ) {
			$genes = $1;
			print "// No of protein_coding genes: $genes\n";
		}
	}
	close GENES;

	open (TRANS, "mysql -hmysql-eg-vectorbase -P4367 -uensro $db -e 'select count(*) from transcript where biotype = \"protein_coding\";' |");
	while (<TRANS>) {
		chomp;
		if ( /^\|\s+(\d+)\s+\|/ ) {
			$trans = $1;
			print "// No of protein_coding transcripts: $trans\n";
		}
	}
	close TRANS;

	open (DOMAINS, "mysql -hmysql-eg-vectorbase -P4367 -uensro $db -e 'SELECT COUNT(DISTINCT gene_id) FROM transcript t, translation tl, protein_feature p WHERE tl.transcript_id = t.transcript_id AND p.translation_id = tl.translation_id AND p.analysis_id IN (SELECT analysis_id FROM analysis_description WHERE display_label IN (\"PROSITE profiles\",\"Prints domain\",\"Pfam domain\",\"TIGRFAM domain\", \"Superfamily domain\", \"SMART domain\", \"Pfam\", \"scanprosite\", \"Superfamily\", \"Smart\", \"Panther\", \"Tigrfam\", \"Gene3d\",\"pfscan\", \"Prints\", \"PIRSF\", \"HAMAP\", \"SMART domain\", \"hmmpanther\", \"PROSITE patterns\", \"gene3d\", \"blastprodom\", \"PIRSF domain\", \"PANTHER domain\", \"Gene3D domain\"));' |");
	while (<DOMAINS>) {
		chomp;
		if ( /^\|\s+(\d+)\s+\|/ ) {
			$domains = $1;
			print "// No of transcripts with IPR domains: $domains\n";
			print "// % coverage of transcripts " . ($domains/$transcripts)*100 . "\n";
		}
	}
	close DOMAINS;







}
#_ end of database loop


__DATA__
aedes_aegypti_core_1312_73_2
anopheles_albimanus_core_1312_73_1
anopheles_arabiensis_core_1312_73_1
anopheles_atroparvus_core_1312_73_1
anopheles_christyi_core_1312_73_1
anopheles_culicifacies_core_1312_73_1
anopheles_darlingi_core_1312_73_2
anopheles_dirus_core_1312_73_1
anopheles_epiroticus_core_1312_73_1
anopheles_farauti_core_1312_73_1
anopheles_funestus_core_1312_73_1
anopheles_gambiaeM_core_1312_73_1
anopheles_gambiaeS_core_1312_73_1
anopheles_gambiae_core_1312_73_3
anopheles_maculatus_core_1312_73_1
anopheles_melas_core_1312_73_1
anopheles_merus_core_1312_73_1
anopheles_minimus_core_1312_73_1
anopheles_quadriannulatus_core_1312_73_1
anopheles_sinensis_core_1312_73_1
anopheles_stephensiI_core_1312_73_2
anopheles_stephensi_core_1312_73_1
biomphalaria_glabrata_core_1312_73_1
culex_quinquefasciatus_core_1312_73_2
drosophila_melanogaster_core_1312_73_546
glossina_morsitans_core_1312_73_1
ixodes_scapularis_core_1312_73_1
lutzomyia_longipalpis_core_1312_73_1
pediculus_humanus_core_1312_73_1
phlebotomus_papatasi_core_1312_73_1
rhodnius_prolixus_core_1312_73_1
__END__
