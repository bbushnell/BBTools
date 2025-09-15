#!/bin/bash
set -e

TAXPATH="auto"

#Fetch
wget -nv http://ftp.arb-silva.de/release_132/Exports/SILVA_132_SSURef_tax_silva_trunc.fasta.gz
wget -nv http://ftp.arb-silva.de/release_132/Exports/SILVA_132_LSURef_tax_silva_trunc.fasta.gz

#Process LSU
time reformat.sh in=SILVA_132_LSURef_tax_silva_trunc.fasta.gz out=lsu_utot.fa.gz fastawrap=4000 utot zl=9 qtrim=rl trimq=1 ow minlen=80 pigz=20
time clumpify.sh in=lsu_utot.fa.gz reorder groups=1 -Xmx31g fastawrap=4000 out=lsu_clumped.fa.gz zl=9 pigz=20
time rename.sh addprefix prefix="lcl\|LSU" in=lsu_clumped.fa.gz out=lsu_prefix.fa.gz pigz=20 zl=9 fastawrap=4000
time gi2taxid.sh -Xmx8g tree=auto in=lsu_prefix.fa.gz out=lsu_renamed.fa.gz silva zl=9 pigz=20 taxpath=$TAXPATH
time sortbyname.sh -Xmx31g in=lsu_renamed.fa.gz out=lsu_sorted.fa.gz ow taxa tree=auto fastawrap=4000 zl=9 pigz=20 allowtemp=f taxpath=$TAXPATH
time dedupe.sh in=lsu_sorted.fa.gz out=lsu_deduped100pct.fa.gz zl=9 pigz=20 ordered fastawrap=4000

#Process SSU
time reformat.sh in=SILVA_132_SSURef_tax_silva_trunc.fasta.gz out=ssu_utot.fa.gz fastawrap=4000 utot zl=9 qtrim=rl trimq=1 ow minlen=80 pigz=20
time clumpify.sh in=ssu_utot.fa.gz reorder groups=1 -Xmx31g fastawrap=4000 out=ssu_clumped.fa.gz zl=9 pigz=20
time rename.sh addprefix prefix="lcl\|SSU" in=ssu_clumped.fa.gz out=ssu_prefix.fa.gz pigz=20 zl=9 fastawrap=4000
time gi2taxid.sh -Xmx8g tree=auto in=ssu_prefix.fa.gz out=ssu_renamed.fa.gz silva zl=9 pigz=20 taxpath=$TAXPATH
time sortbyname.sh -Xmx31g in=ssu_renamed.fa.gz out=ssu_sorted.fa.gz ow taxa tree=auto fastawrap=4000 zl=9 pigz=20 allowtemp=f taxpath=$TAXPATH
time dedupe.sh in=ssu_sorted.fa.gz out=ssu_deduped_sorted.fa.gz zl=9 pigz=20 ordered fastawrap=4000
#time sortbyname.sh -Xmx31g in=ssu_deduped100pct.fa.gz out=ssu_deduped_sorted.fa.gz ow taxa tree=auto fastawrap=4000 zl=9 pigz=20 allowtemp=f

#Merge LSU and SSU
cat ssu_sorted.fa.gz lsu_sorted.fa.gz > both_sorted_temp.fa.gz
time sortbyname.sh -Xmx31g in=both_sorted_temp.fa.gz out=both_sorted.fa.gz ow taxa tree=auto fastawrap=4000 zl=9 pigz=20 allowtemp=f taxpath=$TAXPATH
rm both_sorted_temp.fa.gz
cat ssu_deduped_sorted.fa.gz lsu_deduped_sorted.fa.gz > both_deduped_sorted_temp.fa.gz
time sortbyname.sh -Xmx31g in=both_deduped_sorted_temp.fa.gz out=both_deduped_sorted.fa.gz ow taxa tree=auto fastawrap=4000 zl=9 pigz=20 allowtemp=f taxpath=$TAXPATH
rm both_deduped_sorted_temp.fa.gz
filterbytaxa.sh include=f in=both_deduped_sorted.fa.gz id=2323,256318,12908 tree=auto requirepresent=f out=both_deduped_sorted_no_unclassified_bacteria.fa.gz zl=9 fastawrap=9999 ow

#Sketch steps
time sketchblacklist.sh -Xmx16g in=both_deduped_sorted_no_unclassified_bacteria.fa.gz prefilter=f tree=auto taxa silva taxlevel=species ow out=blacklist_silva_species_500.sketch mincount=500
time sketchblacklist.sh -Xmx16g in=both_deduped_sorted_no_unclassified_bacteria.fa.gz prefilter=f tree=auto taxa silva taxlevel=genus ow out=blacklist_silva_genus_200.sketch mincount=200
time mergesketch.sh -Xmx1g in=blacklist_silva_species_500.sketch,blacklist_silva_genus_200.sketch out=blacklist_silva_merged.sketch
time sketch.sh files=31 out=both_taxa#.sketch in=both_deduped_sorted.fa.gz size=200 maxgenomefraction=0.1 -Xmx8g tree=auto mode=taxa ow silva blacklist=blacklist_silva_merged.sketch autosize ow

time sketch.sh files=31 out=both_seq#.sketch in=both_deduped_sorted.fa.gz size=200 maxgenomefraction=0.1 -Xmx8g tree=auto mode=sequence ow silva blacklist=blacklist_silva_merged.sketch autosize ow parsesubunit

cp blacklist_silva_merged.sketch /global/projectb/sandbox/gaag/bbtools/jgi-bbtools/resources/