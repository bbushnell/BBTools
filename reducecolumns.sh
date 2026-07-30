#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 24, 2025

Description:  Extracts specific columns from a tab-delimited vector file,
producing a new file with only the selected columns.  Useful for removing
input dimensions from ML training data (e.g., dropping codon one-hots to
test whether they contribute signal).  Writes a new #dims header
automatically.  Processes files line-by-line with low memory overhead.

Usage:  reducecolumns.sh <in> <out> <col0> <col1> <col2> ...

Parameters:
in=<file>       Input tab-delimited vector file (with #dims header).
out=<file>      Output file with only the selected columns.
col0,col1,...   Zero-based column indices to keep, in order.
                The LAST column listed should be the target/output column.
                All other columns are treated as inputs.
                Supports three formats:
                  5       Single column (0-indexed).
                  0-8     Range: columns 0 through 8 inclusive.
                  17+     Open range: column 17 through the last column.
                Example: 0-8 11 15 17+  keeps columns 0-8, 11, 15, and
                17 onward.  The last column in the expanded list becomes
                the target.

Input format:
Tab-delimited with a '#dims <inputs> <outputs>' header line.
Lines starting with '#' are skipped (the tool writes its own header).

Output format:
Same tab-delimited format with updated '#dims N 1' header, where N is
the number of selected columns minus 1 (inputs only; last column = output).

Example:
  To keep columns 0-7, 32-33, and 34 from a 35-column file:
  reducecolumns.sh vectors.tsv reduced.tsv 0-7 32-34

  To keep everything from column 10 onward:
  reducecolumns.sh vectors.tsv reduced.tsv 10+

Java Parameters:
-Xmx            Set memory usage (default 2g, usually sufficient).

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

	parseJavaArgs "--xmx=2g" "--xms=2g" "--percent=24" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ml.ReduceColumns $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
