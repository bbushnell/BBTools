#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 26, 2021

Description:  Mutates a reference by applying a set of variants.
When 2 variants overlap, the one with the higher allele count is used.

Usage:  applyvariants.sh in=<input file> vcf=<vcf file> out=<output file>

Standard parameters:
in=<file>       Reference fasta.
vcf=<file>      Variants.
basecov=<file>  Optional per-base coverage from BBMap or Pileup.
out=<file>      Output fasta.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression
                level; lower compression is faster.

Processing parameters:		
mincov=0        If positive and depth is below this, change ref to N.
                Requires a coverage file.
maxindel=-1     If positive, ignore indels longer than this.
noframeshifts=f Ignore indels that are not a multiple of 3 in length.

Renaming parameters:
name=           Optionally rename sequences to this.
addnumbers=f    Add _1 and so forth to ensure sequence names are unique.
prefix=t        Use the name as a prefix to the old name, instead of replacing
                the old name.
delimiter=_     Symbol to place between parts of the new name.
                For space or tab, use the literal word.

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
	SCRIPT="$0"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
	CP="$DIR/current/"
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=4g" "--xms=4g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP var2.ApplyVariants $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
