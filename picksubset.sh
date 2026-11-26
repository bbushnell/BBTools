#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified July 30, 2024

Description:  Greedily selects a diverse subset from an all-to-all
identity comparison.  It removes one endpoint of the highest-ANI
remaining edge until at most files= nodes remain, the highest remaining
ANI is at most ani=, or all supplied edges have been processed.
This program does not use taxonomy or promise a globally optimal subset.

Input should be in 3+ column TSV format (first 3 are required):
(query, ref, ANI)
...as produced by CompareSketch when run like this:
comparesketch.sh ata format=3 includeself perfile records=99999 *.fasta

Query-column names define selectable nodes.  A reference consulted
during selection or final reporting must appear among query rows
actually read, or the run reports the missing reference and query.
The lines= limit can truncate coverage.  #Query/Ref/ANI or
#Query/Ref/AAI headers identify optional KID and SSU columns by
name; headerless input uses legacy KID and SSU positions, columns
10 and 12 when present.

Usage:  picksubset.sh in=<file> out=<file> invalid=<file> files=<number>

Parameters:
in=             Input file comparing all-to-all comparisons.
out=            Output file for the list of files to retain.
invalid=        Output file for the list of files to discard.
overwrite=t     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
files=0         Stop after retaining at most this many files.
ani=0           Stop when the next pairwise ANI would be at most this percent.
lines=0         Stop after this many nonempty physical input records,
                including comments and headers; 0 or negative is unlimited.
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

	parseJavaArgs "--xmx=2g" "--xms=2g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.PickSubset $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" jgi.PickSubset "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
