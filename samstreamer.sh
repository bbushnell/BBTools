#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 15, 2025

Description:  Interconverts sam, bam, fasta, or fastq rapidly.
Sam and bam input also allow filtering options; bam allows bai generation.

Usage:  samstreamer.sh in=<file> out=<file>
        samstreamer.sh <in> <out>
Examples:
samstreamer.sh reads.sam.gz mapped.bam unmapped=f
samstreamer.sh sorted.bam sorted.bai
samstreamer.sh sorted.bam reads.fq.gz 

Filtering parameters:
minpos=         Ignore alignments not overlapping this range.
maxpos=         Ignore alignments not overlapping this range.
minmapq=        Ignore alignments with mapq below this.
maxmapq=        Ignore alignments with mapq above this.
minid=0.0       Ignore alignments with identity below this.
maxid=1.0       Ignore alignments with identity above this.
contigs=        Comma-delimited list of contig names to include. These 
                should have no spaces, or underscores instead of spaces.
                If present, this will be a whitelist.
mapped=t        Include mapped reads.
unmapped=t      Include unmapped reads.
mappedonly=     If true, include only mapped reads.
unmappedonly=   If true, only include unmapped reads.
secondary=t     Include secondary alignments.
supplementary=t Include supplementary alignments.
lengthzero=t    Include alignments without bases.
invert=f        Invert sam filters.
ordered=t       Keep reads in input order.
duplicate=t     Include reads marked as duplicate.
qfail=t         Include reads marked as qfail.
ref=<file>      Optional reference file.


Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
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

	parseJavaArgs "--xmx=2g" "--xms=2g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP stream.SamStreamerWrapper $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
