#!/bin/bash

##  Written by Brian Bushnell
##  Last modified April 11, 2020
##  Description:  Summarizes all the runs of "processCorona.sh" in this directory.
##  Then tars them if you want to send them somewhere.
##  Should be run only after all samples are processed individually.

##  Set minimum coverage for variant calls.
MINCOV=5

##  Specify the viral reference file.
##  NC_045512.fasta contains the SARS-CoV-2 genome, equivalent to bbmap/resources/Covid19_ref.fa
REF="NC_045512.fasta"

##  Call variants in multisample mode to solve things.
##  This reports the genotype of *all* samples at any position at which a variant is called in *any* sample.
callvariants.sh *_deduped_trimclip.sam.gz ref="$REF" multisample out=allVars.vcf ow -Xmx4g usebias=f strandedcov minstrandratio=0 maf=0.6 minreads="$MINCOV" mincov="$MINCOV" minedistmax=30 minedist=16 flagnearby

##  Make a summary of coverage at varous depth cutoffs for all libraries.
summarizecoverage.sh *basecov_border5.txt out=coverageSummary.txt

mkdir output
cp *.sh output
cp *.bam* output
cp *.txt output
cp *.vcf output
cp *genome*.fa output

rm results.tar
tar -cf results.tar output
