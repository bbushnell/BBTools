#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified December 19, 2018
Description:  Generates a length histogram of input reads.

Usage:	readlength.sh in=<input file>

Parameters:
in=<file>    	The 'in=' flag is needed only if the input file is not the first parameter.  'in=stdin.fq' will pipe from standard in.
in2=<file>   	Use this if 2nd read of pairs are in a different file.
out=<file>   	Write the histogram to this file.  Default is stdout.
bin=10       	Set the histogram bin size.
max=80000    	Set the max read length to track.
round=f      	Places reads in the closest bin, rather than the highest bin of at least readlength.
nzo=t        	(nonzeroonly) Do not print empty bins.
reads=-1     	If nonnegative, stop after this many reads.

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
	CP="$DIR/current/"
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=400m" "--xms=400m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.MakeLengthHistogram $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
