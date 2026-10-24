#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 24, 2026

Description:  Reorders reads randomly, keeping pairs together.
Unlike Shuffle, Shuffle2 can write temp files to handle large datasets.
Temporary files are created in the working directory.

Usage:  shuffle2.sh in=<file> out=<file>

Standard parameters:
in=<file>       The 'in=' flag is needed if the input file is not the first parameter.  'in=stdin' will pipe from standard in.
in2=<file>      Use this if 2nd read of pairs are in a different file.
out=<file>      The 'out=' flag is needed if the output file is not the second parameter.  'out=stdout' will pipe to standard out.
out2=<file>     Use this to write 2nd read of pairs to a different file.
overwrite=t     (ow) Set to false to force the program to abort rather than overwrite an existing file.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster.
int=auto        (interleaved) Set to t or f to override interleaving autodetection.

Processing parameters:
seed=-1         A nonnegative seed gives repeatable shuffling for identical
                input, heap size, and settings. Negative seeds use fresh randomness.
deterministic=f Set to t as shorthand for seed=1, or f for seed=-1.
allowtemp=t     Allow temporary files when the in-memory threshold is reached.
memmult=0.35    Fraction of available memory used per input block; must be in (0,1].
maxfiles=16     Maximum input files per merge; must be at least 2.
delete=t        Delete successfully merged temporary files.
minlen=0        Retain a pair when either mate reaches this minimum length.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
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

	parseJavaArgs "--xmx=2g" "--xms=2g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP sort.Shuffle2 $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" sort.Shuffle2 "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
