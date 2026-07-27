#!/bin/bash

usage(){
echo "
Written by Ady
Last modified July 26, 2026

Description:  Builds DynamicDemiLog sketches from FASTA/FASTQ files
and writes them as A48-encoded TSV.  Computes GC content during hashing.
Perfile mode supports multithreading for processing many files in parallel.

Usage:  ddlwriter.sh in=a.fa,b.fa out=ddls.tsv.gz [k=31] [buckets=2048]

Parameters:
in=<file>       Input file(s), comma-delimited.
out=<file>      Output DDL TSV file.
k=31            K-mer length for hashing.
buckets=2048    Number of DDL buckets.
seed=12345      Hash seed.
mode=perfile    One DDL per input file (default).  Multithreaded.
                persequence: one DDL per sequence/contig.
                pertid: merge sequences sharing the same taxonomy ID.
                Each value also works as a bare flag, e.g. 'pertid'.
parsetaxid=t    Extract taxonomy IDs from filenames (tid_NNNN) or
                headers (tid|NNNN).  Set false for anonymous mode.
lineage=t       Attach a taxonomic lineage string to each record, using
                the tax tree.  Costs a one-time tree load here, but saves
                every downstream reader from loading the tree at all.
kmers=f         Also store the source kmers, one per bucket, in the
                output.  Adds a data line per record; readers restore them.
blacklist=<file>  Ignore kmers present in this file while sketching.  The
                filename is recorded in the output header.
exponent=5      Exponent bit width within each 16-bit bucket value; the
                mantissa gets the remaining 16-exponent bits.  Written to
                the output header and restored automatically on load.
                5 suits genomes (an extra mantissa bit, and the cardinalities
                it gives up are unreachable); use 6 for read files.
threads=auto    Number of threads for perfile mode (auto = all cores).
                Peaks at ~32 cores for bgzipped RefSeq input.
overwrite=f     Overwrite existing output.
verbose=f       Print additional progress information.
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

	parseJavaArgs "--xmx=4g" "--xms=400m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ddl.DDLWriter $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
