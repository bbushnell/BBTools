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
#Last updated August 7, 2019

#Fetches and sketches nt.
#Be sure taxonomy is updated first!
#To use this script outside of NERSC, modify $TAXPATH to point to your directory with the BBTools taxonomy data,
#e.g. TAXPATH="/path/to/taxonomy_directory/"

TAXPATH="auto"

#Ensure necessary executables are in your path
#module load pigz

#Fetch nt.
wget -q -O - ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz | gi2taxid.sh -Xmx1g in=stdin.fa.gz out=renamed.fa.gz pigz=32 unpigz bgzip zl=8 server ow shrinknames maxbadheaders=5000 badheaders=badHeaders.txt taxpath=$TAXPATH

#Sort by taxonomy.
#This makes sketching by taxa use much less memory because sketches can be written to disk as soon as they are finished.
time sortbyname.sh -Xmx96g in=renamed.fa.gz out=sorted.fa.gz ow taxa tree=auto fastawrap=1023 zl=9 pigz=32 minlen=60 bgzip unbgzip

#Make a blacklist of kmers occuring in at least 100 different genuses.
time sketchblacklist.sh -Xmx31g in=sorted.fa.gz prepasses=1 tree=auto taxa taxlevel=genus ow out=blacklist_nt_genus_100.sketch mincount=120 k=32,24 taxpath=$TAXPATH

#Generate 31 sketch files, with one sketch per species.
#Multiple files allow sketches to load faster on multicore systems. 
time bbsketch.sh -Xmx31g in=sorted.fa.gz out=taxa#.sketch mode=taxa tree=auto files=31 ow unpigz minsize=300 prefilter autosize blacklist=blacklist_nt_genus_100.sketch k=32,24 depth taxpath=$TAXPATH

#A query such as contigs.fa can now be compared to the new reference sketches like this:
#comparesketch.sh in=contigs.fa k=32,24 tree=auto taxa#.sketch blacklist=blacklist_nt_genus_100.sketch

#On NERSC systems, you can then set the default path to nt by pointing /global/projectb/sandbox/gaag/bbtools/nt/current at the path to the new sketches.
#Then you can use the default set of nt sketches like this:
#comparesketch.sh in=contigs.fa nt tree=auto
#That command automatically adds the default path to the sketches, the blacklist, and the correct values for K.
