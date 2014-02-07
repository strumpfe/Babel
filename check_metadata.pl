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

	print "// Get meta table from database: $db\n" if ($verbose);

	open (LOOK, "mysql -hmysql-eg-vectorbase -P4367 -uensrw -pscr1b3vb $db -e 'select * from meta;' |");
	while (<LOOK>) {
		chomp;
		@f = split/\t/;
		if ( $f[2] eq "assembly.name" )     { print "Assembly\t(assembly.name)\t$f[3]\n"; }
		if ( $f[2] eq "genebuild.version" ) { print "Genebuild\t(genebuild.version)\t$f[3]\n"; }

#		print;
	}

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
