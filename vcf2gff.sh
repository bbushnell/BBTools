#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 13, 2019

Description:  Generates a GFF3 from a VCF.

Usage:  vcf2gff.sh in=<vcf file> out=<gff file>

Parameters:
in=<file>       Input VCF file.
out=<file>      Output GFF file.
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

z="-Xmx200m"
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

gff() {
	local CMD="java $EA $EOOM $z -cp $CP gff.GffLine $@"
#	echo $CMD >&2
	eval $CMD
}

gff "$@"
