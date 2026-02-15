#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 4, 2025

Description:  Finds orfs and calls genes in unspliced prokaryotes.
This includes bacteria, archaea, viruses, and mitochondria.
Can also predict 16S, 18S, 23S, 5S, and tRNAs.

Usage:  callgenes.sh in=contigs.fa out=calls.gff outa=aminos.faa out16S=16S.fa

File parameters:
in=<file>       A fasta file; the only required parameter.
out=<file>      Output gff file.
outa=<file>     Amino acid output.
out16s=<file>   16S output.
model=<file>    A pgm file or comma-delimited list.
                If unspecified a default model will be used.
stats=stderr    Stats output (may be stderr, stdin, a file, or null).
hist=null       Gene length histogram.
compareto=      Optional reference gff file to compare with the gene calls.
                'auto' will name it based on the input file name.

Formatting parameters:
json=false      Print stats in JSON.
binlen=21       Histogram bin length.
bins=1000       Maximum histogram bins.
pz=f            (printzero) Print histogram lines with zero count.



Other parameters:
minlen=60       Don't call genes shorter than this.
trd=f           (trimreaddescription) Set to true to trim read headers after
                the first whitespace.  Necessary for IGV.
merge=f         For paired reads, merge before calling.
detranslate=f   Output canonical nucleotide sequences instead of amino acids.
recode=f        Re-encode nucleotide sequences over called genes, leaving
                non-coding regions unchanged.

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
	CP="$DIR/current/"
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=6g" "--xms=6g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prok.CallGenes $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
