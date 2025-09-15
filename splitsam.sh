#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 9, 2015

Description:  Splits a sam file into three files:
Plus-mapped reads, Minus-mapped reads, and Unmapped.
If 'header' is the 5th argument, header lines will be included.

Usage:  splitsam <input> <plus output> <minus output> <unmapped output>

Input may be stdin or a sam file, raw or gzipped.
Outputs must be sam files, and may be gzipped.
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

function splitsam() {
	local CMD="java $EA $EOOM -Xmx128m -cp $CP jgi.SplitSamFile $@"
	echo $CMD
	eval $CMD
}

splitsam "$@"
