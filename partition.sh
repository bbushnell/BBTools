#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 14, 2026

Description:  Splits a sequence file evenly into multiple files.

Usage:  partition.sh in=<file> out=<outfile> ways=<number>

in2 and out2 are for paired reads and are optional.
If input is paired and out2 is not specified, data will be written interleaved.
Output filenames MUST contain a '%' symbol.  This will be replaced by a number.

Parameters and their defaults:

in=<file>       Input file.
out=<file>      Output file pattern (containing a % symbol, like 'part%.fa').
in2, out2       Optional flags for use with twin fastq files.
ways=-1         The number of output files to create; must be positive.
pacbio=f        Set to true to keep PacBio subreads together.  Only works in
                count mode.
int=f           (interleaved) Determines whether INPUT file is considered interleaved.
zl=4            (ziplevel) Set compression level, 1 (low) to 9 (max).

Mode parameters:
partitionby=<mode>    Partition by metric: count, bp, gc, hh, caga, length, depth
                      count: Round-robin (default)
                      bp: Balance by number of base pairs
                      gc/hh/caga: Split by composition metrics
                      length: Split by sequence length
                      depth: Split by coverage depth
cutoff=<x,y,z>  Custom partition cutoffs (auto-sets ways to cutoffs+1)
cov=<file>      A coverage file from covmaker or pileup, or a sam or bam file,
                used in depth mode; if unset, depth will be parsed from contig
                headers in Tadpole, SPAdes, or MetaHipMer format.


Depth mode options:
cov=<file>      A coverage file from covmaker or pileup, or a sam or bam file,
                used in depth mode; if unset, depth will be parsed from contig
                headers in Tadpole, SPAdes, or MetaHipMer format.

Depth mode examples:
  partition.sh in=contigs.fa partitionby=depth ways=4
    # Auto-balance by depth parsed from contig headers
  partition.sh in=contigs.fa partitionby=depth cov=coverage.txt ways=4
    # Use external coverage file (pileup or covmaker format)
  partition.sh in=contigs.fa partitionby=depth cov=reads.bam ways=3
    # Calculate depth from BAM alignments
  partition.sh in=contigs.fa partitionby=depth cutoff=10,50,200
    # Custom cutoffs creating 4 partitions: <10, 10-50, 50-200, >200


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

	parseJavaArgs "--xmx=1g" "--xms=256m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP scalar.PartitionReads $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
