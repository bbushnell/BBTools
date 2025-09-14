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

#This block allows symlinked shellscripts to correctly set classpath.
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx6g"
z2="-Xms6g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
}
calcXmx "$@"

function callgenes() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP prok.CallGenes $@"
	#Too long to echo sometimes since wildcards can be expanded
	#echo $CMD >&2
	eval $CMD
}

callgenes "$@"
