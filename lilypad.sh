#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified September 13, 2019

Description:  Uses mapped paired reads to generate scaffolds from contigs.
Designed for use with ordinary paired-end Illumina libraries.

Usage:  lilypad.sh in=mapped.sam ref=contigs.fa out=scaffolds.fa

Standard parameters:
in=<file>       Reads mapped to the reference; should be sam or bam.
ref=<file>      Reference; may be fasta or fastq.
out=<file>      Modified reference; should be fasta.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Processing parameters:
gap=10          Pad gaps with a minimum of this many Ns.
mindepth=4      Minimum spanning read pairs to join contigs.
maxinsert=3000  Maximum allowed insert size for proper pairs.
mincontig=200   Ignore contigs under this length if there is a
                longer alternative.
minwr=0.8       (minWeightRatio) Minimum fraction of outgoing edges
                pointing to the same contig.  Lower values will increase
                continuity at a risk of misassemblies.
minsr=0.8       (minStrandRatio) Minimum fraction of outgoing edges
                indicating the same orientation.  Lower values will increase
                continuity at a possible risk of inversions.
passes=8        More passes may increase continuity.
samestrandpairs=f  Read pairs map to the same strand.  Currently untested.

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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP consensus.Lilypad $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
