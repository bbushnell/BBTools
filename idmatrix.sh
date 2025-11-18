#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 25, 2014

Description:  Generates an identity matrix via all-to-all alignment.

*** WARNING: This program may produce incorrect results in some cirumstances.
*** It is not advisable to use until fixed.

Usage:	idmatrix.sh in=<file> out=<file>

Parameters:
in=<file>       File containing reads. in=stdin.fa will pipe from stdin.
out=<file>      Matrix output. out=stdout will pipe to stdout.
threads=auto    (t) Set number of threads to use; default is number of
                logical processors.
percent=f       Output identity as percent rather than a fraction.
edits=          Allow at most this much edit distance.  Default is the
                length of the longest input sequence. Lower is faster.
width=          Alignment bandwidth, lower is faster.  Default: 2*edits+1.
usejni=f        (jni) Do alignments faster, in C code.  Requires
                compiling the C code; details are in /jni/README.txt.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding automatic
                memory detection. -Xmx20g will specify
                20 gigs of RAM, and -Xmx200m will specify 200 megs.
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

	parseJavaArgs "--xmx=3200m" "--xms=3200m" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.IdentityMatrix $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
