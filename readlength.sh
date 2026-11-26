#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 24, 2026
Description:  Generates a length histogram of input reads.

Usage:	readlength.sh in=<input file>

Parameters:
in=<file>       Input sequence file; in=stdin.fq reads standard input.
                The first bare argument can also specify the input.
in2=<file>      Optional second mate file.
out=<file>      Histogram output; default stdout. Aliases: hist, lhist.
bin=10          Positive integer bin width in bases; labels are multiples of bin.
max=80000       Cap bins at floor(max/bin)*bin; longer reads remain in that last bin.
round=f         f: round lengths down; t: nearest bin, halfway ties upward.
nzo=t           Omit zero-count bins. With f, print empty bins in the scanned range.
reads=-1        Positive input fragment limit (reads or pairs); negative unlimited.
                Zero is unsupported. The histogram counts each mate separately.
overwrite=t     Allow replacement of an existing histogram file.
append=f        Append a complete report (raw/gzip only); requires overwrite=f.
                Appending does not merge the reports' bins.

Min, max and average use actual lengths. Median, mode and population standard
deviation use binned/capped lengths. Row base totals retain actual bases.
Cumulative columns include the current and all larger bins.

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

	parseJavaArgs "--xmx=400m" "--xms=400m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.MakeLengthHistogram $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" jgi.MakeLengthHistogram "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
