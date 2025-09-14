#!/bin/bash
set -e

#Written by Brian Bushnell and William Andreopoulos
#For Illumina reads


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


#First link reference as ref.fa and reads as reads.fq.gz

kmercountexact.sh in=reads.fq.gz khist=khist_raw.txt peaks=peaks_raw.txt

primary=`grep "haploid_fold_coverage" peaks_raw.txt | sed "s/^.*\t//g"`
cutoff=$(( $primary * 3 ))

bbnorm.sh in=reads.fq.gz out=highpass.fq.gz pigz passes=1 bits=16 min=$cutoff target=9999999
reformat.sh in=highpass.fq.gz out=highpass_gc.fq.gz maxgc=0.45

kmercountexact.sh in=highpass_gc.fq.gz khist=khist_100.txt k=100 peaks=peaks_100.txt smooth ow smoothradius=1 maxradius=1000 progressivemult=1.06 maxpeaks=16 prefilter=2

mitopeak=`grep "main_peak" peaks_100.txt | sed "s/^.*\t//g"`

upper=$((mitopeak * 6 / 3))
lower=$((mitopeak * 3 / 7))
mcs=$((mitopeak * 3 / 4))
mincov=$((mitopeak * 2 / 3))

tadwrapper.sh in=highpass_gc.fq.gz out=contigs_intermediate_%.fa k=78,100,120 outfinal=contigs_intermediate.fa prefilter=2 mincr=$lower maxcr=$upper mcs=$mcs mincov=$mincov

bbduk.sh in=highpass.fq.gz ref=contigs_intermediate.fa outm=bbd005.fq.gz k=31 mm=f mkf=0.05

tadpole.sh in=bbd005.fq.gz out=contigs_bbd.fa prefilter=2 mincr=$((mitopeak * 3 / 8)) maxcr=$((upper * 2)) mcs=$mcs mincov=$mincov k=100 bm1=6

ln -s contigs_bbd.fa contigs.fa
