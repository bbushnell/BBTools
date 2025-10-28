#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 21, 2015

Description:  Counts GC content of reads or scaffolds.

Usage:  countgc in=<input> out=<output> format=<format>

Input may be stdin or a fasta or fastq file, compressed or uncompressed.
Output (which is optional) is tab-delimited.
Parameters:
format=1:   name   length   A   C   G   T   N
format=2:   name   GC
format=4:   name   length   GC
Note that in format 1, A+C+G+T=1 even when N is nonzero.

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

z="-Xmx120m"
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

countgc() {
	local CMD="java $EA $SIMD $EOOM $z -cp $CP jgi.CountGC $@"
	echo $CMD >&2
	eval $CMD
}

countgc "$@"
