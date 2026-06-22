#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and UMP45
Last modified June 19, 2026

Description:  Performs set operations on BED files:
Union, intersection, and subtraction.  Reports base-pair coverage
statistics: bp covered by each input, shared bp, and bp unique to each.

BED coordinates are 0-based half-open; intervals are sorted and merged
per scaffold on load, so inputs may be unsorted or self-overlapping.

Usage:  bedset.sh in=<file,file,...> out=<file>

I/O parameters:
in=<file>       Input; at least 1 file (two or more for set operations).  BED
                files, or VCF files (auto-detected by .vcf/.vcf.gz and converted
                to padded variant-span intervals).
out=<file>      Output BED file (optional; stats always print to stderr).
overwrite=t     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Mode Parameters (choose one only):
subtract=t      Subtract all other files from the first file (default).
union=f         Make a union of all files.
intersection=f  Make an intersection of all files (region covered by ALL).

VCF Input Parameters (apply when an input is a .vcf/.vcf.gz file):
pad=0           Pad each variant's reference span by this many bp on each side,
                so deletions are fully covered with margin.
multiallelic=f  Keep only sites whose first-sample genotype is multiallelic
                (an allele index >=2) - e.g. to build a multiallelic-exclusion BED.

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

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
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

	parseJavaArgs "--xmx=4g" "--xms=4g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP var2.BedSet $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
