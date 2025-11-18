#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified July 30, 2024

Description:  Selects a subset of files from an all-to-all identity
comparison.  The subset will contain exactly X files with maximal
pairwise ANI, or all files with at most Y pairwise identity.
This program is similar to representative.sh but does not use taxonomy.

Input should be in 3+ column TSV format (first 3 are required):
(query, ref, ANI)
...as produced by CompareSketch when run like this:
comparesketch.sh ata format=3 includeself perfile records=99999 *.fasta

Usage:  picksubset.sh in=<file> out=<file> invalid=<file> files=<number>

Parameters:
in=             Input file comparing all-to-all comparisons.
out=            Output file for the list of files to retain.
invalid=        Output file for the list of files to discard.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
files=0         Number of files to retain.
ani=0           Maximum pairwise ANI allowed, expressed as a percent.
NOTE: files or ani, or both, must be set.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will
                specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                The max is typically around 85% of physical memory.
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

	parseJavaArgs "--xmx=2g" "--xms=2g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.PickSubset $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
