#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 5, 2022

Description:  Parses verbose output from bloomfilter.sh for a specific paper.
Irrelevant for most people, but useful for reproducing published results.
You use it to parse output from bloomfilter.sh and tabulate it.

Usage:  bloomfilterparser.sh in=<input file> out=<output file>

...where the input file is whatever bloomfilter.sh prints to the screen.  E.G.
in=slurm-3249652.out out=summary.txt

You get details of calls to increment() if you add the verbose flag.

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

z="-Xmx300m"
z="-Xms300m"
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

parsebloom() {
	local CMD="java $EA $EOOM $z -cp $CP bloom.ParseBloomFilter $@"
	echo $CMD >&2
	eval $CMD
}

parsebloom "$@"
