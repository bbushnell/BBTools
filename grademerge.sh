#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 24, 2026

Description:  Grades merged-read lengths against insert sizes encoded in headers.
Length agreement is checked; sequence bases are not compared. Recognizes BBTools
SYN headers, insert= headers, and RandomReadsMG headers.

Usage:  grademerge.sh in=<file>
        grademerge.sh <readable-input-file>

Parameters:
in=<file>       Merged-read input file, or 'stdin'. Use in= when a filename
                contains '='; a positional input cannot contain '='. Put a
                positional filename first, or use in= after other options.
raw=<file>      Optional unmerged input for additional total and mergeable counts.
                raw1 is an alias. A '#' in this filename expands to mate 1 and 2.
raw2=<file>     Optional second raw-read file; used only with raw or raw1.
verbose=f      Print additional diagnostics.

Raw-input counts refer to primary records returned by the reader. A positive
encoded insert smaller than the total length of the read and any mate contributes to the
mergeable count. Merged-read length comparisons use the in= records.

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

	parseJavaArgs "--xmx=200m" "--xms=200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.GradeMergedReads $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" jgi.GradeMergedReads "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
