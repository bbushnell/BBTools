#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 10, 2018

Description:  Merges .pgm files.

Usage:  mergepgm.sh in=x.pgm,y.pgm out=z.pgm

File parameters:
in=<file,file>  A pgm file or comma-delimited list of pgm files.
out=<file>      Output filename.
normalize=f     Merge proportionally to base counts, so small models
                have equal weight to large models.  Normalization happens
                before applying the @ multiplier.
@ symbol        Input filenames in the form of 'x.pgm@0.1' will have
                a multiplier applied to that model prior to merging.

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

	parseJavaArgs "--xmx=200m" "--xms=200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prok.PGMTools $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
