#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified September 27, 2018

Description:  Generates a prokaryotic gene model (.pkm) for gene calling.
Input is fasta and gff files.
The .pkm file may be used by CallGenes.

Usage:  analyzegenes.sh in=x.fa gff=x.gff out=x.pgm

File parameters:
in=<file>       A fasta file or comma-delimited list of fasta files.
gff=<file>      A gff file or comma-delimited list.  This is optional;
                if present, it must match the number of fasta files.
                If absent, a fasta file 'foo.fasta' will imply the
                presence of 'foo.gff'.
out=<file>      Output pgm file.

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

z="-Xmx2g"
z2="-Xms2g"
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

function analyze() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP prok.AnalyzeGenes $@"
	#echo $CMD >&2
	eval $CMD
}

analyze "$@"
