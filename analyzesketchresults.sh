#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified December 19, 2019

Description:  Analyzes sketch results from query, ref, ani format.

Usage:  analyzesketchresults.sh in=<file> out=<outfile>


Parameters and their defaults:

in=<file>           Required input file of Sketch results in 3column format.
in2=<file>          Optional second input file of Sketch results in amino mode.
out=stdout.txt      Output file for summary of per-tax-level averages.
outaccuracy=<file>  Output file for accuracy results; requires query taxIDs and printcal.
outmap=<file>       Output file for ANI vs AAI.  Requires in2.
bbsketch            Parse BBSketch output format (default).
mash                Parse Mash output format.  Files should be named like this:
                    tid_511145_Escherichia_coli_str._K-12_substr._MG1655.fa.gz
blast               Parse Blast output format (TODO).

ow=f                (overwrite) Overwrites files that already exist.
app=f               (append) Append to files that already exist.

Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.

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

	parseJavaArgs "--xmx=300m" "--xms=300m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP sketch.AnalyzeSketchResults $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
