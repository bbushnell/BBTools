#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 16, 2026

Description:  Validates a directory of paired genome (.fna) + annotation (.gff)
              files.  For each basename with both an .fna(.gz) and a .gff(.gz):
                1) Checks every GFF feature line has exactly 9 tab-separated
                   fields; a short count flags a truncated/partial record.
                2) Cross-checks sequence IDs BOTH directions: every seqid the
                   GFF references (feature seqids + ##sequence-region pragmas)
                   must appear as an FNA header, and vice-versa.
              Also flags unpaired files and empty files.  Multithreaded across
              pairs.  Complements a gzip-integrity check (validategz.sh): gzip
              -t catches truncated compression, this catches content problems
              (partial GFF records, FNA/GFF version mismatch) in files that
              decompress cleanly.

Usage:  validatepairs.sh in=<directory> [out=report.txt]

Parameters:
in=<dir,dir>    Input director(ies) to scan for .fna/.gff pairs.  A bare
                directory path may also be given as a positional argument.
out=<file>      Write the per-pair FAIL report here (default: stderr).
printpass=f     Also print a line for each passing pair.
bidirectional=t Report FNA seqids absent from the GFF (the 'vice-versa'
                direction).  Set f to only report GFF seqids absent from the
                FNA (the always-an-error direction).  Aliases: fnaingff.
maxbad=10       Max example seqids/line-numbers to list per failure reason.
t=              Worker threads (default: all available).
ow=t            Overwrite the output file if it exists.

Java Parameters:
-Xmx            Set Java memory usage, overriding autodetection.
-eoom           Exit if an out-of-memory exception occurs.  Requires Java 8u92+.
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

	parseJavaArgs "--xmx=4g" "--xms=4g" "--percent=42" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP gff.ValidatePairs $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
