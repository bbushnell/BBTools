#!/bin/bash

usage(){
echo "
Written by Brian Bushnell (tool authored by UMP45)
Last modified June 12, 2026

Description:  Reads a neural network and writes it back out, purely to convert
the on-disk format.  The output coding follows the current default (A48), so an
old decimal .bbnet becomes A48 -- lossless for the stored weights and roughly
30% smaller.  The weights are unchanged; this is a format round-trip, not
re-training.

Usage:  netconvert.sh in=<old.bbnet> out=<new.bbnet>

in=<file>       Input network (.bbnet); decimal or A48 (auto-detected).
out=<file>      Output network (.bbnet).
a48=t           Write A48 coding (default).  Set a48=f to write decimal.
overwrite=t     (ow) Permit overwriting the output file.

Java Parameters:
-Xmx            Set Java memory usage; e.g. -Xmx1g.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
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

	parseJavaArgs "--xmx=2000m" "--xms=2000m" "--percent=10" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ml.NetConvert $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
