#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified July 16, 2018

Description:  Compresses sequence data into a fasta file containing each kmer
exactly once.  Allows arbitrary kmer set operations via multiple passes.

Usage:  kcompress.sh in=<reads> out=<contigs> min=<1> max=<2147483647>

Input parameters:
in=<file>           Primary input file for reads to use as kmer data.
in2=<file>          Second input file for paired data.
reads=-1            Only process this number of reads, then quit (-1 means all).

Output parameters:
out=<file>          Write contigs (in contig mode).
showstats=t         Print assembly statistics after writing contigs.
fuse=0              Fuse output sequences into chunks at least this long,
                    padded with 1 N between sequences.

Prefiltering parameters:
prefilter=0         If set to a positive integer, use a countmin sketch
                    to ignore kmers with depth of that value or lower.
prehashes=2         Number of hashes for prefilter.
prefiltersize=0.2   (pff) Fraction of memory to use for prefilter.
minprobprefilter=t  (mpp) Use minprob for the prefilter.
prepasses=1         Use this many prefiltering passes; higher be more thorough
                    if the filter is very full.  Set to 'auto' to iteratively 
                    prefilter until the remaining kmers will fit in memory.

Hashing parameters:
k=31                Kmer length (1 to 31).
prealloc=t          Pre-allocate memory rather than dynamically growing; 
                    faster and more memory-efficient.  A float fraction (0-1)
                    may be specified; default is 1.
minprob=0.5         Ignore kmers with overall probability of correctness below this.
minprobmain=t       (mpm) Use minprob for the primary kmer counts.
threads=X           Spawn X threads (default is number of logical processors).

Assembly parameters:
mincount=1          (min) Only retain kmers that occur at least this many times.
maxcount=BIG        (max) Only retain kmers that occur at most this many times.
requiresamecount    (rsc) Only build contigs from kmers with exactly the same count.
rcomp=t             Store forward and reverse kmers together.  Setting this to
                    false will only use forward kmers.


Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.
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

	parseJavaArgs "--xmx=14g" "--xms=14g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP assemble.KmerCompressor $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
