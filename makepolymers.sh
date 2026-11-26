#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 3, 2017

Description:  Creates repeated-motif FASTA sequences.
For each k in the requested range, this writes 4^k records: one for every
nucleotide k-mer, repeated until the sequence length is the smallest multiple of
k that is at least minlen+k-1 bases.
Can be used in conjunction with mutate.sh to generate low-complexity sequence.

Usage:  makepolymers.sh out=<output file> k=<repeat length> minlen=<sequence length>

I/O parameters:
out=<file>      Output FASTA file; required.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Processing parameters:
k=1             Length of repeated nucleotide units.
                To generate a sweep of multiple values of k,
                specify both mink and maxk.
                Valid range: 1<=mink<=maxk<=31.
minlen=31       Minimum sequence length; must be at least 1.
                Actual sequence length is the smallest multiple of k at least
                minlen+k-1, exposing every starting phase of a minlen-base
                window in each repeated motif.

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

	parseJavaArgs "--xmx=4g" "--xms=4g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.MakePolymers $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" jgi.MakePolymers "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
