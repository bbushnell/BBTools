#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 28, 2026

Description:  Estimates cardinality of unique kmers in sequence data.
See also kmercountmulti.sh.

Usage:  loglog.sh in=<file> k=<31>

Parameters:
in=<file>       (in1) Input file, or comma-delimited list of files.
in2=<file>      Optional second file for paired reads.
k=31            Use this kmer length for counting.
buckets=2048    Use this many buckets for counting; higher decreases
                variance, for large datasets.  Must be a power of 2.
seed=-1         Use this seed for hash functions.  A negative number forces
                a random seed.
minprob=0       Set to a value between 0 and 1 to exclude kmers with a lower
                probability of being correct.


Shortcuts:
The # symbol will be substituted for 1 and 2.
For example:
loglog.sh in=read#.fq
...is equivalent to:
loglog.sh in1=read1.fq in2=read2.fq

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

Supported input formats are fastq, fasta, scarf, sam, and bam.
Supported compression formats are gzip and bz2.
To read from stdin, set 'in=stdin'.  The format should be specified with an extension, like 'in=stdin.fq.gz'

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

	parseJavaArgs "--xmx=200m" "--xms=200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP cardinality.LogLogWrapper $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
