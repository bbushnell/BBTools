#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 2, 2017

Description:  Remove Smart Bell adapters from PacBio reads.

Usage:        removesmartbell in=<input> out=<output> split=t

Input may be fasta or fastq, compressed or uncompressed (not H5 files).

Parameters:
in=file         Specify the input file, or stdin.
out=file        Specify the output file, or stdout.
adapter=        Specify the adapter sequence (default is normal SmrtBell).
split=t            t: Splits reads at adapters.
                   f: Masks adapters with X symbols.

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

z="-Xmx400m"
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

removesmartbell() {
	local CMD="java $EA $EOOM $z -cp $CP pacbio.RemoveAdapters2 $@"
	echo $CMD >&2
	eval $CMD
}

removesmartbell "$@"
