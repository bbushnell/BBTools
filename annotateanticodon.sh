#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Noire
Last modified August 13, 2026

Description:  Annotates tRNA features in a GFF with the anticodon triplet.
NCBI GFFs give the anticodon by genomic position (anticodon=(pos:X..Y)) and the
amino acid by name (product=tRNA-Xxx), but not the triplet as text.  This tool
reads the 3 bases at that position from the paired genome fasta, reverse-
complements them on the minus strand, converts DNA->RNA (T->U), and appends
Note=tRNA-Xxx(YYY) to the feature.  All lines pass through unchanged except tRNA
lines, which gain the Note= attribute (idempotent; already-annotated tRNAs are
left as-is).  This is the format TrnaConsensusBuilder.parseAnticodon reads first,
so annotated GFF flows through cutgff -> TrnaConsensusBuilder unchanged.

Usage:  annotateanticodon.sh in=<fna file> gff=<gff file> out=<gff file>

gff= is optional; the gff filename is assumed from the fasta name if omitted.
For many genomes, pass multiple fna and an output directory:

annotateanticodon.sh in=archaea4/tid_1.fna.gz,archaea4/tid_2.fna.gz outdir=staging/

File Parameters:
in=<file>           Input FNA (fasta) file(s); comma-delimited for multiple.
gff=<file>          Input GFF file(s) (optional; inferred from fna name).
out=<file>          Output (annotated) GFF file.  Single input only.
outdir=<dir>        Output directory; writes one annotated GFF per input,
                    named after its input GFF.  Required for multiple inputs.
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

	parseJavaArgs "--xmx=1g" "--xms=1g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP gff.AnnotateAnticodon $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
