#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified Jan 7, 2020

Description:  Shrinks sketches to a smaller fixed length.

Please read bbmap/docs/guides/BBSketchGuide.txt for more information.

Usage:       subsketch.sh in=file.sketch out=sub.sketch size=1000 autosize=f
Bulk usage:  subsketch.sh in=big#.sketch out=small#.sketch sizemult=0.5

Standard parameters:
in=<file>       Input sketch file containing one or more sketches.
out=<file>      Output sketch file.
size=10000      Size of sketches to generate, if autosize=f.
autosize=t      Autosize sketches based on genome size.
sizemult=1      Adjust default sketch autosize by this factor.
blacklist=      Apply a blacklist to the sketch before resizing.
files=31        If the output filename contains a # symbol,
                spread the output across this many files, replacing
                the # with a number.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

For more detailed information, please read /bbmap/docs/guides/BBSketchGuide.txt.
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

	parseJavaArgs "--xmx=4g" "--xms=4g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP sketch.SubSketch $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
