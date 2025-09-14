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

#Sketches RefSeq.
#Be sure taxonomy is updated first!
#To use this script outside of NERSC, modify $TAXPATH to point to your directory with the BBTools taxonomy data,
#e.g. TAXPATH="/path/to/taxonomy_directory/"

TAXPATH="auto"

#Ensure necessary executables are in your path
#module load pigz

#Sort by taxonomy.
#This makes sketching by taxa use much less memory because sketches can be written to disk as soon as they are finished.
#If this stage still runs out of memory please increase the -Xmx flag or add the flag 'maxfiles=4'.
#Or you can just skip it; sorting should no longer be necessary.
time sortbyname.sh -Xmx100g in=renamed.fa.gz memmult=0.33 out=sorted.fa.gz zl=9 pigz=64 bgzip taxa tree=auto gi=ignore fastawrap=1023 minlen=60 readbufferlen=2 readbuffers=1 taxpath=$TAXPATH

#If sortbyname runs out of memory in the last stage (when sorted.fa.gz is being written), you can usually run this:
#mergesorted.sh -Xmx100g sort_temp* bgzip unbgzip zl=8 out=sorted.fa.gz taxa tree=auto fastawrap=1023 minlen=60 readbufferlen=2 readbuffers=1 ow taxpath=$TAXPATH


#Make a blacklist of kmers occuring in at least 140 different genuses.
time sketchblacklist.sh -Xmx63g in=sorted.fa.gz prepasses=1 tree=auto taxa taxlevel=species ow out=blacklist_refseq_genus_140.sketch mincount=140 k=32,24 sizemult=2 taxpath=$TAXPATH

#Generate 31 sketch files, with one sketch per species.
time bbsketch.sh -Xmx63g in=sorted.fa.gz out=taxa#.sketch mode=taxa tree=auto accession=null gi=null files=31 ow unpigz minsize=400 prefilter autosize blacklist=blacklist_refseq_genus_140.sketch k=32,24 depth sizemult=2 taxpath=$TAXPATH

#A query such as contigs.fa can now be compared to the new reference sketches like this:
#comparesketch.sh in=contigs.fa k=32,24 tree=auto taxa*.sketch blacklist=blacklist_refseq_genus_140 taxpath=$TAXPATH

#On NERSC systems, you can then set the default path to nt by pointing /global/projectb/sandbox/gaag/bbtools/refseq/current at the path to the new sketches.
#Then you can use the default set of refseq sketches like this:
#comparesketch.sh in=contigs.fa refseq tree=auto
#That command automatically adds the default path to the sketches, the blacklist, and the correct values for K.
