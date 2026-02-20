#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Contributor: Noire
Last modified February 9, 2026

Description:  Processes multiple genome files individually through Tadpole assembler
while preserving taxonomic ID labels. Reads genome files with taxID in filename
(pattern: tid_<number>_...) or FASTA headers. By default, uses append mode where
each genome's contigs are appended directly to the output file with unique contig IDs.

This tool eliminates the need for coassembly, preventing chimeric contigs and
simplifying the workflow for metagenomic binning evaluation datasets.

Usage:  reassemble.sh in=<files> out=<contigs> k=<kmer>

Standard parameters:
in=<file>           Input files. Comma-delimited, directories, and wildcards supported.
out=<file>          Output file for assembled contigs.
k=<int>             Kmer length for Tadpole assembly (required).

Reassemble-specific parameters:
failfast=f          Abort on first failure (default: false, continue processing).
tempdir=<path>      Use temporary files instead of append mode. If specified, each genome
                    assembles to a temp file, then all are concatenated. If null (default),
                    Tadpole appends directly to output file (more efficient).
delete=t            Delete temporary files after concatenation (default: true).
                    Only relevant if tempdir is specified.
verbose=f           Verbose logging (default: false).

All other parameters are passed through to Tadpole. Common Tadpole parameters:
mcs=1               minCountSeed (default: 1 in code for sparse genomes).
mce=1               minCountExtend (default: 1 in code for sparse genomes).
mincontig=1         Minimum contig length (default: 1 in code).
prefilter=0         Prefilter level.
mode=contig         Assembly mode (contig/extend/correct).

Usage examples:

# Basic usage with directory input
reassemble.sh in=genomes/ out=assembled.fa k=155

# Comma-delimited file list
reassemble.sh in=tid_123.fa,tid_456.fa out=output.fa k=31

# With custom parameters and fail-fast mode
reassemble.sh in=genomes/*.fa out=output.fa k=155 mcs=2 mce=2 mincontig=200 failfast=t

Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM. The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an out-of-memory exception occurs.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

if [ "$1" = "help" ] || [ "$1" = "-help" ] || [ "$1" = "--help" ] || [ "$1" = "-h" ] || [ "$1" = "--h" ]; then
	usage
	exit 0
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

	parseJavaArgs "--xmx=8g" "--xms=8g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP assemble.Reassemble $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
