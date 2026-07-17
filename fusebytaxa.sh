#!/bin/bash

usage(){
echo "
Written by Noire
Last modified July 16, 2026

Description:  Fuses contigs sharing a taxonomy ID (from 'tid|TAXID|...' headers)
into N-padded genome-scale sequences, emitting only fused sequences at least
'floor' bp long.  Turns a combined RefSeq clade file (mostly short rRNA markers
plus some assemblies) into assembly-attempt sequences for contig-classification
benchmarking, dropping marker-only taxa.

Requires taxonomically-grouped input (contigs of the same taxid adjacent), which
RefSeq clade files usually already are; otherwise sort first with sortbyname.sh taxa=t.

Usage:  fusebytaxa.sh in=<file> out=<file> floor=1m pad=10 maxlen=50m

Parameters:
in=<file>       Input sequences (fasta/fastq, may be gzipped).
out=<file>      Output fused sequences.
floor=1m        Only emit fused sequences at least this long.  Set 0 to emit all
                (e.g. for viruses, which are small).
pad=10          Number of Ns inserted between fused contigs.
maxlen=50m      Split a taxon's fused sequence into chunks no longer than this.
overwrite=f     Overwrite existing output.
-Xmx            Set Java memory; streaming, so little is needed.
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

	parseJavaArgs "--xmx=2g" "--xms=400m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP synth.FuseByTaxa $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
