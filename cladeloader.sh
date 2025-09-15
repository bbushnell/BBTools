#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 4, 2025

Description:  Loads fasta files and writes clade files.

Usage: cladeloader.sh in=contigs.fa out=clades.clade

Parameters:
in=<file,file>  Fasta files with tid in headers.
out=<file>      Output file.
maxk=5          Limit max kmer length (range 3-5).
a48             Output counts in ASCII-48 instead of decimal.
16s=<file,file> Optional tax-labeled file of 16S sequences.
18s=<file,file> Optional tax-labeled file of 16S sequences.
replaceribo     Set true if existing ssu should be replaced by new ones.
usetree=f       Load a taxonomic tree to generate lineage strings.
aligner=quantum Options include ssa2, glocal, drifting, banded, crosscut.

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

z="-Xmx4g"
z2="-Xms4g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 4000m 42
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

vector() {
	local CMD="java $EA $EOOM $z -cp $CP bin.CladeLoader $@"
	echo $CMD >&2
	eval $CMD
}

vector "$@"
