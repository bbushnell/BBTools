#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified June 6, 2016

Description:  Summarizes the output of GradeMerge for comparing 
read-merging performance.

Usage:  summarizemerge.sh in=<file>

Parameters:
in=<file>       A file containing GradeMerge output.

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

summarizemerge() {
	local CMD="java $EA $EOOM $z -cp $CP driver.ProcessSpeed $@"
#	echo $CMD >&2
	eval $CMD
}

summarizemerge "$@"
