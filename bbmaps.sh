#!/bin/bash

usage(){
echo "
BBMapS — Streamer/Writer-based BBMap (bbmapnova).
Functionally equivalent to bbmap.sh: same alignment engine (BBMapThread/
AbstractMapThread, unchanged), same flag surface, same output. The only
difference is internal I/O plumbing — worker threads call stream.Streamer/
stream.Writer directly instead of the classic ConcurrentReadInputStream/
ConcurrentReadOutputStream pair. See align2/BBMapS.java.

Usage:  bbmaps.sh ref=<fasta> in=<reads> out=<sam> t=32

Same flag surface as bbmap.sh (build=, in=, in2=, ref=, t=, out=, etc.).
Run bbmap.sh -h for the full flag reference.
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

	parseJavaArgs "--xmx=3200m" "--xms=3200m" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP align2.BBMapS build=1 overwrite=true fastareadlen=500 $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
