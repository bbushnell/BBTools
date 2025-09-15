#!/bin/bash

##  Written by Brian Bushnell
##  Last modified April 11, 2020
##  Description:  Outer wrapper for processing a bunch of Covid samples in a directory.
##  This is just a template that needs to be modified before use, according to the input file names.

echo "This template must be modifed before use according to your file names."
exit

##  Optionally, get rid of old files first, if the pipeline is being rerun.
rm *.sam.gz *.bam *.bai *.txt *_genome.fa *_adapters.fa *.vcf *.vcf.gz

##  Generate quality-score calibration matrices.
##  This only needs to be run on one sample.
sh ./recal.sh Sample1

##  Add a line like this for each interleaved PE file named, for example, Sample1.fq.gz
##  Alternately you could put some kind of loop here, depending on your naming convention.
sh processCorona.sh Sample1
sh processCorona.sh Sample2
##  etc.

##  Summarize the output if there are multiple libraries (optional).
sh makeSummary.sh 1>makeSummary.o 2>&1
