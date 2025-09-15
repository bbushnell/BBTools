#!/bin/bash

##  Written by Brian Bushnell
##  Last modified April 11, 2020
##  Description:  Calls SARS-CoV-2 variants from Illumina amplicon data.
##                This script assumes input data is paired-end.

##  Usage:  processCorona.sh <prefix>
##  For example, "processCorona.sh Sample1" if the data is in Sample1.fq.gz

##  Grab the sample name from the command line.
NAME="$1"

##  Set minimum coverage for genotype calls.
##  Areas below this depth will be set to N in the consensus genome.
MINCOV=5

##  Specify the viral reference file.
##  NC_045512.fasta contains the SARS-CoV-2 genome, equivalent to bbmap/resources/Covid19_ref.fa
REF="NC_045512.fasta"

##  PCR primers
##  artic3 primers are in /bbmap/resources/artic3.fasta
PRIMERS="artic3.fasta"

##  This line is in case the script is being re-run, to clear the old output.
rm "$NAME"*.sam.gz "$NAME"*.bam "$NAME"*.bai "$NAME"*.txt "$NAME"*.fa "$NAME"*.vcf

##  If data is paired in twin files, interleave it into a single file.
##  Otherwise, skip this step.
##  In this case, the files are assumed to be named "Sample1_R1.fq.gz" and "Sample1_R2.fq.gz"
#reformat.sh in="$NAME"_R#.fq.gz out="$NAME".fq.gz

##  Split into Covid and non-Covid reads if this has not already been done.
##  This step can be skipped if non-Covid was already removed.
bbduk.sh ow -Xmx1g in="$NAME".fq.gz ref="$REF" outm="$NAME"_viral.fq.gz outu="$NAME"_nonviral.fq.gz k=25

##  Recalibrate quality scores prior to any trimming.
##  *** Requires a recalibration matrix in the working directory (see recal.sh for details). ***
##  This step is optional but useful for Illumina binned quality scores.
bbduk.sh in="$NAME"_viral.fq.gz out="$NAME"_recal.fq.gz recalibrate -Xmx1g ow

##  Discover adapter sequence for this library based on read overlap.
##  You can examine the adapters output file afterward if desired;
##  If there were too few short-insert pairs this step will fail (and you can just use the default Illumina adapters).
bbmerge.sh in="$NAME"_recal.fq.gz outa="$NAME"_adapters.fa ow reads=1m

##  Remove duplicates by sequence similarity.
##  This is more memory-efficient than dedupebymapping.
clumpify.sh in="$NAME"_recal.fq.gz out="$NAME"_clumped.fq.gz zl=9 dedupe s=2 passes=4 -Xmx31g

##  Perform adapter-trimming on the reads.
##  Also do quality trimming and filtering.
##  If desired, also do primer-trimming here by adding, e.g., 'ftl=20' to to trim the leftmost 20 bases.
##  If the prior adapter-detection step failed, use "ref=adapters"
bbduk.sh in="$NAME"_clumped.fq.gz out="$NAME"_trimmed.fq.gz minlen=60 ktrim=r k=21 mink=9 hdist=2 hdist2=1 ref="$NAME"_adapters.fa altref=adapters maq=14 qtrim=r trimq=10 maxns=0 tbo tpe ow -Xmx1g ftm=5

##  Trim artic3 primers, if present.
##  Disable this line if you are not amplifying with primers.
bbduk.sh in="$NAME"_trimmed.fq.gz out="$NAME"_trimmed2.fq.gz ref="$PRIMERS" ktrim=l restrictleft=30 k=22 hdist=3 qhdist=1 rcomp=f mm=f

##  Align reads to the reference.
##  Local flag is due to primer-amplification-induced anomalies at read ends;
##  for randomly-sheared, unamplified data, "local" should be omitted.
bbmap.sh ref="$REF" in="$NAME"_trimmed2.fq.gz outm="$NAME"_mapped.sam.gz nodisk local maxindel=500 -Xmx4g ow k=12

##  Deduplicate based on mapping coordinates.
##  Note that if you use single-ended amplicon data, you will lose most of your data here.
dedupebymapping.sh in="$NAME"_mapped.sam.gz out="$NAME"_deduped.sam.gz -Xmx31g ow

##  Remove junk reads with unsupported unique deletions; these are often chimeric.
filtersam.sh ref="$REF" ow in="$NAME"_deduped.sam.gz out="$NAME"_filtered.sam.gz mbad=1 del sub=f mbv=0 -Xmx4g

##  Remove junk reads with multiple unsupported unique substitutions; these are often junk, particularly on Novaseq.
##  This step is not essential but reduces noise.
filtersam.sh ref="$REF" ow in="$NAME"_filtered.sam.gz out="$NAME"_filtered2.sam.gz mbad=1 sub mbv=2 -Xmx4g

##  Trim soft-clipped bases.
bbduk.sh in="$NAME"_filtered2.sam.gz trimclip out="$NAME"_trimclip.sam.gz -Xmx1g ow

##  Call variants from the sam files.
##  The usebias=f/minstrandratio=0 flags are necessary due to amplicon-induced strand bias,
##  and should be removed if the data is exclusively shotgun/metagenomic or otherwise randomly fragmented.
callvariants.sh in="$NAME"_trimclip.sam.gz ref="$REF" out="$NAME"_vars.vcf -Xmx4g ow strandedcov usebias=f minstrandratio=0 maf=0.6 minreads="$MINCOV" mincov="$MINCOV" minedistmax=30 minedist=16 flagnearby

##  Calculate reduced coverage as per CallVariants defaults (ignoring outermost 5bp of reads).
pileup.sh in="$NAME"_trimclip.sam.gz basecov="$NAME"_basecov_border5.txt -Xmx4g ow border=5

##  Generate a mutant reference by applying the detected variants to the reference.
##  This is essentially the reference-guided assembly of the strain.
##  Also changes anything below depth MINCOV to N (via the mindepth flag).
##  Does not apply indels below MINCOV
applyvariants.sh in="$REF" out="$NAME"_genome.fa vcf="$NAME"_vars.vcf basecov="$NAME"_basecov_border5.txt ow mindepth="$MINCOV"

##  Make bam/bai files; requires samtools to be installed.
##  This step is only necessary for visualization, not variant-calling.
#samtools view -bShu "$NAME"_trimclip.sam.gz | samtools sort -m 2G -@ 3 - -o "$NAME"_sorted.bam
#samtools index "$NAME"_sorted.bam

##  At this point, "$NAME"_sorted.bam, "$NAME"_sorted.bam.bai, "$REF", and "$NAME"_vars.vcf can be used for visualization in IGV.

