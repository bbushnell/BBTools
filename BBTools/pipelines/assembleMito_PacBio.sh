#!/bin/bash
set -e

#Written by Brian Bushnell and William Andreopoulos
#For error-corrected PacBio reads


#Setup

if [[ $NERSC_HOST == genepool ]]; then
	module load bbtools
	module load quast
	module load jgi-rqc-synbio
elif [[ $NERSC_HOST == denovo ]]; then
	#TODO
elif [[ $NERSC_HOST == cori ]]; then
	#TODO
fi


#First link reads as reads.fa in the working directory, then run the shellscript.

kmercountexact.sh in=reads.fa khist=khist_raw.txt peaks=peaks_raw.txt

primary=`grep "haploid_fold_coverage" peaks_raw.txt | sed "s/^.*\t//g"`
cutoff=$(( $primary * 3 ))

bbnorm.sh in=reads.fa out=highpass.fa pigz passes=1 bits=16 min=$cutoff target=9999999
reformat.sh in=highpass.fa out=highpass_gc.fa maxgc=0.4

#Optional machine-learning filtration here

kmercountexact.sh in=highpass_gc.fa khist=khist_124.txt k=124 peaks=peaks_124.txt smooth ow smoothradius=1 maxradius=1000 progressivemult=1.06 maxpeaks=16 prefilter=2

mitopeak=`grep "main_peak" peaks_124.txt | sed "s/^.*\t//g"`

upper=$((mitopeak * 5 / 3))
lower=$((mitopeak * 1 / 2))
mcs=$((mitopeak * 4 / 5))
mincov=$((mitopeak * 3 / 4))

tadwrapper.sh in=highpass_gc.fq.gz out=contigs_intermediate_%.fa k=93,124,155 outfinal=contigs_intermediate.fa prefilter=2 mincr=$lower maxcr=$upper mcs=$mcs mincov=$mincov

tadpole.sh in=filtered.fa out=contigs93.fa prefilter=2 mincr=$lower maxcr=$upper mcs=$mcs mincov=$mincov k=93
tadpole.sh in=filtered.fa out=contigs124.fa prefilter=2 mincr=$lower maxcr=$upper mcs=$mcs mincov=$mincov k=124
tadpole.sh in=filtered.fa out=contigs155.fa prefilter=2 mincr=$lower maxcr=$upper mcs=$mcs mincov=$mincov k=155

bbduk.sh in=highpass.fa ref=contigs_intermediate.fa outm=bbd005.fa k=31 mm=f mkf=0.05

tadpole.sh in=bbd005.fa out=contigs_bbd.fa prefilter=2 mincr=$((mitopeak * 3 / 7)) maxcr=$((upper * 2)) mcs=$mcs mincov=$mincov k=155

ln -s contigs_bbd.fa contigs.fa
