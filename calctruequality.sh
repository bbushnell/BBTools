#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 8, 2024

Description:  Calculates observed quality scores from mapped sam/bam files.
Generates matrices for use in recalibrating quality scores.  By default, 
the matrices are written to /ref/qual/ in the current directory.

If you have multiple sam/bam files demultiplexed from a single sequencing run,
it is recommended to use all of them as input for increased statistical power.
Once the matrices are generated, recalibration can be done on mapped or
unmapped reads; you may get better results by recalibrating the fastq and 
remapping the calibrated reads.

Note!  Diploid organisms with a high heterozygousity rate will induce
inaccurate recalibration at the high end of the quality scale unless SNP
locations are masked or variations are called.  For example, recalibrating 
human reads mapped to an unmasked human reference would generate an 
expected maximal Q-score of roughly 30 due to the human 1/1000 SNP rate.
Variations can be ignored by using the callvars flag or providing
a file of variations.

Usage: calctruequality.sh in=<sam,sam,...sam> path=<directory>

Step 1.  Generate matrices as above.
Step 2.  Recalibrate reads (any kind of files):
bbduk.sh in=<file> out=<file> recalibrate


Parameters (and their defaults)

Input parameters:
in=<file,file>      Sam/bam file or comma-delimited list of files.  Alignments
                    must use = and X cigar symbols, or have MD tags, or
                    ref must be specified.
reads=-1            Stop after processing this many reads (if positive).
samstreamer=t       (ss) Load reads multithreaded to increase speed.
unpigz=t            Use pigz to decompress.

Output parameters:
overwrite=t         (ow) Set to true to allow overwriting of existing files.
path=.              Directory to write quality matrices (within /ref subdir).
write=t             Write matrices.
showstats=t         Print a summary.
pigz=f              Use pigz to compress.

Other parameters:
t=auto              Number of worker threads.
passes=2            Recalibration passes, 1 or 2.  2 is slower but gives more
                    accurate quality scores.
recalqmax=42        Adjust max quality scores tracked.  The actual highest
                    quality score allowed is recalqmax-1.
trackall=f          Track all available quality metrics and produce all
                    matrices, including the ones that are not selected for 
                    quality adjustment.  Reduces speed, but allows testing the
                    effects of different recalibration matrices.
indels=t            Include indels in quality calculations.
usetiles=f          Use per-tile quality statistics to generate matrices.
                    If this is true, the flag must also be used during
                    recalibration (e.g. in BBDuk).

Variation calling parameters:
varfile=<file>      Use the variants in this var file, instead of calling
                    variants.  The format can be produced by CallVariants.
vcf=<file>          Use the variants in this VCF file, instead of
                    calling variants.
callvars=f          Call SNPs, and do not count them as errors.
ploidy=1            Set the organism's ploidy.
ref=                Required for variation-calling.

*** 'Variant-Calling Cutoffs' flags in callvariants.sh are also supported ***

Matrix-selection parameters:
loadq102=           For each recalibration matrix, enable or disable that matrix with t/f.
                    You can specify pass1 or pass2 like this: loadq102_p1=f loadq102_p2=t.
                    The default is loadqbp_p1=t loadqbp_p2=t loadqb123_p=t.
clearmatrices=f     If true, clear all the existing matrix selections.  For example:
                    'clearmatrices loadqbp_p1'
                    This would ignore defaults and select only qbp for the first pass.

Avaliable matrix type parameters:
q102                Quality, leading quality, trailing quality.
qap                 Quality, average quality, position.
qbp                 Quality, current base, position.
q10                 Quality, leading quality.
q12                 Quality, trailing quality.
qb12                Quality, leading base, current base.
qb012               Quality, two leading bases, current base.
qb123               Quality, leading base, current base, trailing base.
qb234               Quality, current base, two trailing bases.
q12b12              Quality, trailing quality, leading base, current base.
qp                  Quality, position.
q                   Current quality score only.


Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
"
}

#This block allows symlinked shellscripts to correctly set classpath.
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx2g"
z2="-Xms2g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

calctruequality() {
	local CMD="java $EA $SIMD $EOOM $z $z2 -cp $CP jgi.CalcTrueQuality $@"
	echo $CMD >&2
	eval $CMD
}

calctruequality "$@"
