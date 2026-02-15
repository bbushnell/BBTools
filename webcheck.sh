#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified December 19, 2018

Description:  Parses a webcheck log.
Input is expected to look like this:
Tue Apr 26 16:40:09 2016|https://rqc.jgi-psf.org/control/|200 OK|0.61

Usage:  webcheck.sh <input files>


Standard parameters:
in=<file>       Primary input.  Can use a wildcard (*) if 'in=' is omitted.
out=<file>      Summary output; optional.
fail=<file>     Output of failing lines; optional.
invalid=<file>  Output of misformatted lines; optional.
extendedstats=f (es) Print more stats.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
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

	parseJavaArgs "--xmx=1g" "--xms=1g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP driver.ProcessWebcheck $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
