#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 21, 2015

Description:  BLASTs queries against reference, and filters out hits with
              scores less than 'cutoff'.  The score is taken from column 12
              of the BLAST output.  The specific BLAST command is:
              blastall -p blastn -i QUERY -d REFERENCE -e 0.00001 -m 8

Usage:  estherfilter.sh <query> <reference> <cutoff>

For example:

estherfilter.sh reads.fasta genes.fasta 1000 > results.txt

'fasta' can be used as a fourth argument to get output in Fasta format.  Requires more memory.

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

z="-Xmx3200m"
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

estherfilter() {
	local CMD="java $EA $EOOM $z -cp $CP driver.EstherFilter $@"
	echo $CMD >&2
	eval $CMD
}

estherfilter "$@"
