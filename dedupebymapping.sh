#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified April 4, 2020

Description:  Deduplicates mapped reads based on pair mapping coordinates.

Usage:   dedupebymapping.sh in=<file> out=<file>

Parameters:
in=<file>           The 'in=' flag is needed if the input file is not the
                    first parameter.  'in=stdin' will pipe from standard in.
out=<file>          The 'out=' flag is needed if the output file is not the
                    second parameter.  'out=stdout' will pipe to standard out.
overwrite=t         (ow) Set to false to force the program to abort rather
                    than overwrite an existing file.
ziplevel=2          (zl) Set to 1 (lowest) through 9 (max) to change
                    compression level; lower compression is faster.
keepunmapped=t      (ku) Keep unmapped reads.  This refers to unmapped
                    single-ended reads or pairs with both unmapped.
keepsingletons=t    (ks) Keep all pairs in which only one read mapped.  If
                    false, duplicate singletons will be discarded.
ignorepairorder=f   (ipo) If true, consider reverse-complementary pairs
                    as duplicates.

Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an
                    out-of-memory exception occurs.  Requires Java 8u92+.
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

	parseJavaArgs "--xmx=3g" "--xms=3g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.DedupeByMapping $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
