#!/bin/bash

##  Written by Brian Bushnell
##  Last modified May 21, 2020
##  Description:  Creates a recalibration matrix for Illumina data.

##  Usage:  recal.sh <prefix>
##  For example, "recal.sh Sample1" if the data is in Sample1.fq.gz

##  This script creates quality-score recalibration matrices for processing Illumina PE reads.
##  It needs to be run ONCE on a single library (preferably a large one) from a sequencing run.
##  Then the primary script will use the recalibration matrices for all of the libraries,
##  assuming all the processing is done in the same directory.
##  This script assumes input data is single-ended or paired and interleaved.
##  If data is paired in twin files, you can run "reformat.sh in1=r1.fq.gz in2=r2.fq.gz out=both.fq.gz" to interleave it.


##  Grab the sample name from the command line
NAME="$1"

##  Set this to read length
READLEN=150

##  Specify the viral reference file.
##  NC_045512.fasta contains the SARS-CoV-2 genome, equivalent to bbmap/resources/Covid19_ref.fa
REF="NC_045512.fasta"

##  Discover adapter sequence for this library based on read overlap.
##  This step should be skipped for single-ended reads.
bbmerge.sh in="$NAME".fq.gz outa="$NAME"_adapters.fa ow reads=1m

##  Adapter-trim and discard everything with adapter sequence so the reads are uniform length.
bbduk.sh -Xmx1g in="$NAME".fq.gz out=recal.fq.gz minlen="$READLEN" ktrim=r k=21 mink=9 hdist=2 hdist2=1 ref="$NAME"_adapters.fa altref=adapters ow tbo tpe

##  Note - if the reads are single-ended, use this command instead:
#bbduk.sh -Xmx1g in="$NAME".fq.gz out=trimmed.fq.gz minlen="$READLEN" ktrim=r k=21 mink=9 hdist=2 hdist2=1 ref=adapters ow

##  Map the reads with very high sensitivity.
bbmap.sh ref="$REF" in=recal.fq.gz outm=mapped.sam.gz vslow -Xmx6g ow

##  Discover true variants.
callvariants.sh in=mapped.sam.gz ref="$REF" out=recal.vcf -Xmx6g ow

##  Generate recalibration matrices.
callvariants.sh in=mapped.sam.gz ref="$REF" out=recal.vcf -Xmx6g ow usebias=f strandedcov minstrandratio=0

##  Now the recalibration matrices are stored in ./ref
##  BBDuk can be run with the 'recal' flag to recalibrate data (mapped or unmapped).
##  It should be run from the directory containing /ref
