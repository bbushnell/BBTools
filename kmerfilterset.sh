#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 24, 2019

Description:  Generates a set of kmers such that every input sequence will
contain at least one kmer in the set.  This is a greedy algorithm which
retains the top X most common kmers each pass, and removes the sequences
matching those kmers, so subsequent passes are faster.

This will not generate an optimally small set but the output will be
quite small.  The file size may be further decreased with kcompress.sh.

Usage:  kmerfilterset.sh in=<input file> out=<output file> k=<integer>

File parameters:
in=<file>       Primary input.
out=<file>      Primary output.
temp=<file>     Temporary file pattern (optional).  Must contain # symbol.
initial=<file>  Initial kmer set (optional).  This can be used to accelerate
                the process.

Processing parameters:
k=31            Kmer length.
rcomp=t         Consider forward and reverse-complement kmers identical.
minkpp=1        (minkmersperpass) Retain at least this many kmers per pass.
                Higher is faster but results in a larger set.
maxkpp=2        (maxkmersperpass) Retain at most this many kmers per pass;
                0 means unlimited.
mincount=1      Ignore kmers seen fewer than this many times in this pass.
maxpasses=3000  Don't run more than this many passes.
maxns=BIG       Ignore sequences with more than this many Ns.
minlen=0        Ignore sequences shorter than this.

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
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=1000m" "--xms=1000m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.KmerFilterSetMaker $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
