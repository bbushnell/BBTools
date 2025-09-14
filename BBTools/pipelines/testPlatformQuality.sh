#!/bin/bash
set -e

#Written by Brian Bushnell
#Last updated February 21, 2018

#This script is designed to evaluate the quality of a new Illumina platform such as NovaSeq.
#Input is assumed to be interleaved and paired.
#This is just an example of a P.heparinus library.


#If necessary, first demultiplex by barcodes.
#demuxbyname.sh in=all.fq.gz delimiter=space suffixmode out=%.fq.gz

#Link input files
ln -s reads.fq.gz raw.fq.gz
ln -s P.heparinus.fa ref.fa

#Remove PhiX, unless that is the intended organism.
#PhiX interferes with adapter sequence detection since it uses different adapters.
bbduk.sh -Xmx1g in=raw.fq.gz out=filtered.fq.gz ref=phix k=31

#Optionally determine actual adapter sequence
#This may fail if there are not enough overlapping reads.
bbmerge.sh in=filtered.fq.gz outa=adapters.fa reads=1m strict

#Trim adapters
#"ref=adapters" will use BBMap's default set of adapters.
#You can alternately use the custom adapters discovered above, if the file contents look reasonable (you may want to trim a trailing poly-A or poly-C).
bbduk.sh -Xmx1g in=filtered.fq.gz out=trimmed.fq.gz ref=adapters k=23 mink=11 tbo tpe ktrim=r ow

#Map to reference
#A reference can be constructed with Tadpole if none is available
bbmap.sh -Xmx8g in=trimmed.fq.gz ref=ref.fa maxindel=2000 slow out=mapped.sam.gz covstats=covstats.txt covhist=covhist.txt 32bit bhist=bhist.txt qhist=qhist.txt aqhist=aqhist.txt bqhist=bqhist.txt lhist=lhist.txt ihist=ihist.txt ehist=ehist.txt qahist=qahist.txt indelhist=indelhist.txt mhist=mhist.txt gchist=gchist.txt idhist=idhist.txt delcov=f ow unpigz=t pigz=t zl=6

#Calculate duplicate rates, both optical and full
#You may need remove groups=1 if you don't have enough memory.
#For fair comparisons, it may be useful to subsample to a fixed number of reads or use the first X reads in a file when calculating dupe rates.
#These approaches give different results.  Using the "First X reads" is important for detecting optical and tile-edge duplicates.
clumpify.sh groups=1 in=trimmed.fq.gz dedupe 1>dedupe_normal.o 2>&1
clumpify.sh groups=1 in=trimmed.fq.gz dedupe optical dist=12k 1>dedupe_optical.o 2>&1

#Call variants
callvariants.sh in=mapped.sam.gz ref=ref.fa out=vars.txt vcf=vars.vcf ploidy=1 -Xmx8g ow

#Generate calibration matrices
calctruequality.sh in=mapped.sam.gz vcf=vars.vcf -Xmx8g

#Recalibrate quality scores
bbduk.sh ow in=mapped.sam.gz out=recal.sam.gz recalibrate -Xmx8g ow

#Generate quality statistics ignoring real variants
bbduk.sh ow in=mapped.sam.gz bhist=bhist.txt qhist=qhist.txt aqhist=aqhist.txt bqhist=bqhist.txt ehist=ehist.txt qahist=qahist.txt indelhist=indelhist.txt mhist=mhist.txt idhist=idhist.txt -Xmx8g vcf=vars.vcf ow

#Generate recalibrated statistics to see if recalibration can generate correct quality scores
bbduk.sh ow in=recal.sam.gz qhist=qhist2.txt aqhist=aqhist2.txt bqhist=bqhist2.txt qahist=qahist2.txt mhist=mhist2.txt idhist=idhist2.txt -Xmx8g vcf=vars.vcf ow

#See what is in the library
#You can iterate a few times by pulling out the primary organisms with BBDuk then rerunning SendSketch
sendsketch.sh in=trimmed.fq.gz reads=500k

#Map the reads to various references you found to see how much contamination there is.
seal.sh -Xmx31g in=trimmed.fq.gz ref=all_fused.fa stats=sealstats.txt ambig=toss clearzone=10 ow

