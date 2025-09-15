#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 5, 2020

Description:  Processes data.  (WIP)

Usage:  runhmm.sh in=<file> out=<file>

Parameters and their defaults:

ow=f                    (overwrite) Overwrites files that already exist.

Processing Parameters:

None yet!


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

function runhmm() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP hmm.HMMSearchReport $@"
	echo $CMD >&2
	eval $CMD
}

runhmm "$@"
