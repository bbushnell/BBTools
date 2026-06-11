#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified June 11, 2026

Description:  Inverts a VCF file produced by mutate.sh.
Swaps ref/alt alleles, flips INS/DEL types, and adjusts coordinates
from original-genome space to mutant-genome space.

This is useful for evaluating variant-calling accuracy using real reads.
The workflow is:
1. Start with a real genome and real reads from it.
2. Mutate the genome with mutate.sh to produce a mutant and a VCF.
3. Map the real reads to the mutant genome and call variants.
4. Invert the mutate.sh VCF with this tool.
5. Compare the inverted VCF (ground truth) to the called variants.

Example workflow:
mutate.sh in=genome.fa out=mutant.fa vcf=mutations.vcf id=0.95
bbmap.sh ref=mutant.fa in=real_reads.fq out=mapped.sam
callvariants.sh in=mapped.sam ref=mutant.fa out=called.vcf ploidy=1
invertvcf.sh in=mutations.vcf out=truth.vcf
gradevcf.sh in=called.vcf truth=truth.vcf ref=mutant.fa

Usage:  invertvcf.sh in=<input vcf> out=<output vcf>

I/O parameters:
in=<file>       Input VCF from mutate.sh.
out=<file>      Output inverted VCF.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
"
}

if [ -z "$1" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	usage
	exit
fi

resolveSymlinks(){
	SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=4g" "--xms=4g" "--percent=42" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP var2.InvertVCF $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
