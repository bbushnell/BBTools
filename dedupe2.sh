#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Jonathan Rood
Last modified September 15, 2015

Dedupe2 is identical to Dedupe except it supports hashing unlimited kmer
prefixes and suffixes per sequence.  Dedupe supports at most 2 of each,
but uses slightly more memory.  You can manually set the number of kmers to
hash per read with the numaffixmaps (nam) flag.  Dedupe will automatically
call Dedupe2 if necessary (if nam=3 or higher) so this script is no longer
necessary.

For documentation, please consult dedupe.sh; syntax is identical.
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

	parseJavaArgs "--xmx=1g" "--xms=3200m" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.Dedupe2 $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
