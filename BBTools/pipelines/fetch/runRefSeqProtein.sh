#!/bin/bash
#SBATCH -J sketch_refseq
#SBATCH -q genepool
#SBATCH -A gtrqc
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -t 71:00:00
#SBATCH --error=log_%j.err
#SBATCH --output=log_%j.out
#SBATCH --exclusive

set -e

#Written by Brian Bushnell
#Last updated August 19, 2019

#Filters prokaryotic clades, translates them to protein, and sketches them.
#To use this script outside of NERSC, modify $TAXPATH to point to your directory with the BBTools taxonomy data,
#e.g. TAXPATH="/path/to/taxonomy_directory/"

TAXPATH="auto"

#Ensure necessary executables are in your path
#module load pigz

mkdir prot
time filterbytaxa.sh in=sorted.fa.gz out=prot/prok.fa.gz fastawrap=4095 ids=Viruses,Bacteria,Archaea,plasmids tree=auto -Xmx16g include taxpath=$TAXPATH zl=9 requirepresent=f
cd prot

wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/protozoa/*genomic.fna.gz | gi2taxid.sh -Xmx1g in=stdin.fa.gz out=protozoa.fa.gz zl=9 server ow
wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/*genomic.fna.gz | gi2taxid.sh -Xmx1g in=stdin.fa.gz out=mito.fa.gz zl=9 server ow
wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/*genomic.fna.gz | gi2taxid.sh -Xmx1g in=stdin.fa.gz out=chloro.fa.gz zl=9 server ow

time callgenes.sh in=prok.fa.gz outa=prok.faa.gz -Xmx16g ow ordered=f zl=9
time callgenes.sh in=protozoa.fa.gz outa=protozoa.faa.gz -Xmx16g ow ordered=f zl=9
time callgenes.sh in=mito.fa.gz outa=mito.faa.gz -Xmx16g ow ordered=f zl=9
time callgenes.sh in=chloro.fa.gz outa=chloro.faa.gz -Xmx16g ow ordered=f zl=9

cat prok.faa.gz mito.faa.gz chloro.faa.gz protozoa.faa.gz > all.faa.gz

time sketchblacklist.sh -Xmx63g in=prok.faa.gz prepasses=1 tree=auto taxa taxlevel=family ow out=blacklist_prokprot_family_40.sketch mincount=40 k=9,12 sizemult=3 amino taxpath=$TAXPATH
time sketchblacklist.sh -Xmx63g in=prok.faa.gz prepasses=1 tree=auto taxa taxlevel=genus ow out=blacklist_prokprot_genus_80.sketch mincount=80 k=9,12 sizemult=3 amino taxpath=$TAXPATH
mergesketch.sh -Xmx1g in=blacklist_prokprot_genus_80.sketch,blacklist_prokprot_family_40.sketch out=blacklist_prokprot_merged.sketch amino name0=blacklist_prokprot_merged k=9,12 ow
time bbsketch.sh -Xmx63g in=all.faa.gz out=taxa#.sketch mode=taxa tree=auto files=31 ow unpigz minsize=200 prefilter autosize blacklist=blacklist_prokprot_merged.sketch k=9,12 depth sizemult=3 amino taxpath=$TAXPATH
