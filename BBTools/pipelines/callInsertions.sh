#!/bin/bash
set -e

#Written by Brian Bushnell
#Last updated February 21, 2018

#This script is designed to call insertions longer than raw read length by lengthening reads prior to alignment.
#It was designed for 2x100bp reads with >=40x coverage.
#Some numbers and steps may need adjustment for different data.
#For large genomes, tadpole and bbmerge may need the flag "prefilter=2" to avoid running out of memory


#Reorder reads for speed of subsequent phases
clumpify.sh in=reads.fq.gz out=clumped.fq.gz dedupe optical

#Remove low-quality reads by position
filterbytile.sh in=clumped.fq.gz out=filtered_by_tile.fq.gz

#Trim adapters and discard reads with Ns
bbduk.sh in=filtered_by_tile.fq.gz out=trimmed.fq.gz maxns=0 ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=85 ref=adapters ftm=5 ordered

#Remove synthetic artifacts and spike-ins by kmer-matching.
bbduk.sh in=trimmed.fq.gz out=filtered.fq.gz k=31 ref=artifacts,phix ordered

#Error-correct phase 1
bbmerge.sh in=filtered.fq.gz out=ecco.fq.gz ecco mix vstrict ordered ihist=ihist_merge1.txt

#Error-correct phase 2
clumpify.sh in=ecco.fq.gz out=eccc.fq.gz passes=4 ecc unpair repair

#Error-correct phase 3
tadpole.sh in=eccc.fq.gz out=ecct.fq.gz ecc ordered

#Read extension
tadpole.sh in=ecct.fq.gz out=extended.fq.gz ordered mode=extend el=20 er=20 k=62

#Read merging
bbmerge-auto.sh in=extended.fq.gz out=merged.fq.gz outu=unmerged.fq.gz rem k=81 extend2=120 zl=8 ordered

#Alignment; only use merged reads
bbmap.sh in=merged.fq.gz out=merged.sam.gz slow bs=bs.sh pigz unpigz ref=reference.fa

#Variant-calling; ploidy may need adjustment.  For a large dataset "prefilter" may be needed.
#To call only insertions, "calldel=f callsub=f" can be added.  Not calling substitutions saves memory.
callvariants.sh in=merged.sam.gz out=vars.txt vcf=vars.vcf.gz ref=ref.fa ploidy=1

#Generate a bam file, if vieing in IGV is desired.
sh bs.sh

