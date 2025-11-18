#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 27, 2017

Description:  Prints sequence gc content once per interval.

Usage:  plotgc.sh in=<input file> out=<output file>

Parameters:
in=<file>       Input file. in=stdin.fa will pipe from stdin.
out=<file>      Output file.  out=stdout will pipe to stdout.
interval=1000   Interval length.
offset=0        Position offset.  For 1-based indexing use offset=1.
psb=t           (printshortbins) Print gc content for the last bin of a contig
                even when shorter than interval.

Java Parameters:

-Xmx            This will set Java's memory usage, overriding automatic
                memory detection. -Xmx20g will 
                specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  
                The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
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

	parseJavaArgs "--xmx=1400m" "--xms=1400m" "--percent=42" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP driver.PlotGC $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
