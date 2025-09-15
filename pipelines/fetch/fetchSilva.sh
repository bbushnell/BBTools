#!/bin/bash
set -e

#Written by Brian Bushnell
#Last updated August 6, 2019

#Fetches and sketches Silva.
#Be sure taxonomy is updated first!
#To use this script outside of NERSC, modify $TAXPATH to point to your directory with the BBTools taxonomy data,
#e.g. TAXPATH="/path/to/taxonomy_directory/"
#Also, check Silva's archive (https://www.arb-silva.de/no_cache/download/archive/) for a newer version.

TAXPATH="auto"

#Fetch files.

#"_tax_silva_trunc" means it has taxonomic information, and it is truncated to retain only alignable 16S bases, not additional sequence.
wget -nv http://ftp.arb-silva.de/release_132/Exports/SILVA_132_SSURef_tax_silva_trunc.fasta.gz
wget -nv http://ftp.arb-silva.de/release_132/Exports/SILVA_132_LSURef_tax_silva_trunc.fasta.gz


#Process SSU

#This transforms U to T, trims leading or trailing Ns, discards sequences under 80bp, and changes the line wrap to 4000 to save space.
time reformat.sh in=SILVA_132_SSURef_tax_silva_trunc.fasta.gz out=ssu_utot.fa.gz fastawrap=4000 utot zl=9 qtrim=rl trimq=1 ow minlen=80 pigz=20

#Clumpify step is not strictly necessary
time clumpify.sh in=ssu_utot.fa.gz reorder groups=1 -Xmx31g fastawrap=4000 out=ssu_clumped.fa.gz zl=9 pigz=20 bgzip

#Prefix with the subunit identifier to allow unique names when SSU and LSU are merged
time rename.sh addprefix prefix="lcl\|SSU" in=ssu_clumped.fa.gz out=ssu_prefix.fa.gz pigz=20 bgzip zl=9 fastawrap=4000

#Rename with the prefix tid|number based on the TaxID.
time gi2taxid.sh -Xmx8g tree=auto gi=null accession=null in=ssu_prefix.fa.gz out=ssu_renamed.fa.gz silva zl=9 pigz=20 bgzip taxpath=$TAXPATH

#Sort by taxonomy to put records from the same organism together,a nd increase compression
time sortbyname.sh -Xmx31g in=ssu_renamed.fa.gz out=ssu_sorted.fa.gz ow taxa tree=auto fastawrap=4000 zl=9 pigz=20 bgzip allowtemp=f taxpath=$TAXPATH

#Remove duplicate or contained sequences
time dedupe.sh in=ssu_sorted.fa.gz out=ssu_deduped_sorted.fa.gz zl=9 pigz=20 bgzip ordered fastawrap=4000


#Process LSU

time reformat.sh in=SILVA_132_LSURef_tax_silva_trunc.fasta.gz out=lsu_utot.fa.gz fastawrap=4000 utot zl=9 qtrim=rl trimq=1 ow minlen=80 pigz=20 bgzip
time clumpify.sh in=lsu_utot.fa.gz reorder groups=1 -Xmx31g fastawrap=4000 out=lsu_clumped.fa.gz zl=9 pigz=20 bgzip
time rename.sh addprefix prefix="lcl\|LSU" in=lsu_clumped.fa.gz out=lsu_prefix.fa.gz pigz=20 bgzip zl=9 fastawrap=4000
time gi2taxid.sh -Xmx8g tree=auto gi=null accession=null in=lsu_prefix.fa.gz out=lsu_renamed.fa.gz silva zl=9 pigz=20 bgzip taxpath=$TAXPATH
time sortbyname.sh -Xmx31g in=lsu_renamed.fa.gz out=lsu_sorted.fa.gz ow taxa tree=auto fastawrap=4000 zl=9 pigz=20 bgzip allowtemp=f
time dedupe.sh in=lsu_sorted.fa.gz out=lsu_deduped100pct.fa.gz zl=9 pigz=20 bgzip ordered fastawrap=4000
#time sortbyname.sh -Xmx31g in=lsu_deduped_sorted.fa.gz out=lsu_deduped_sorted.fa.gz ow taxa tree=auto fastawrap=4000 zl=9 pigz=20 bgzip allowtemp=f taxpath=$TAXPATH


#Merge SSU and LSU

mergesorted.sh -Xmx31g ssu_sorted.fa.gz lsu_sorted.fa.gz bgzip unbgzip zl=9 out=both_sorted.fa.gz taxa tree=auto fastawrap=4000 minlen=60 ow taxpath=$TAXPATH
mergesorted.sh -Xmx31g ssu_deduped_sorted.fa.gz lsu_deduped_sorted.fa.gz bgzip unbgzip zl=9 out=both_deduped_sorted.fa.gz taxa tree=auto fastawrap=4000 minlen=60 ow taxpath=$TAXPATH

#Old version of merge
#cat ssu_sorted.fa.gz lsu_sorted.fa.gz > both_sorted_temp.fa.gz
#time sortbyname.sh -Xmx31g in=both_sorted_temp.fa.gz out=both_sorted.fa.gz ow taxa tree=auto fastawrap=4000 zl=9 pigz=20 bgzip allowtemp=f taxpath=$TAXPATH
#rm both_sorted_temp.fa.gz
#cat ssu_deduped_sorted.fa.gz lsu_deduped_sorted.fa.gz > both_deduped_sorted_temp.fa.gz
#time sortbyname.sh -Xmx31g in=both_deduped_sorted_temp.fa.gz out=both_deduped_sorted.fa.gz ow taxa tree=auto fastawrap=4000 zl=9 pigz=20 bgzip allowtemp=f taxpath=$TAXPATH
#rm both_deduped_sorted_temp.fa.gz


#Sketch steps

#Make a blacklist of kmers occuring in at least 500 different species.
#A blacklist is HIGHLY recommended for Silva or any ribo database.
#This needs to be copied to bbmap/resources/
time sketchblacklist.sh -Xmx31g in=both_deduped_sorted.fa.gz prefilter=f tree=auto taxa silva taxlevel=species ow out=blacklist_silva_species_500.sketch mincount=500 taxpath=$TAXPATH

#Make one sketch per sequence
time sketch.sh files=31 out=both_seq#.sketch in=both_deduped_sorted.fa.gz size=200 maxgenomefraction=0.1 -Xmx8g tree=auto mode=sequence ow silva blacklist=blacklist_silva_species_500.sketch autosize ow parsesubunit taxpath=$TAXPATH

#Make one sketch per TaxID
time sketch.sh files=31 out=both_taxa#.sketch in=both_deduped_sorted.fa.gz size=200 maxgenomefraction=0.1 -Xmx8g tree=auto mode=taxa ow silva blacklist=blacklist_silva_species_500.sketch autosize ow taxpath=$TAXPATH
