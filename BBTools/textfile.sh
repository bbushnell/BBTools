#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 17, 2015

Description:  Displays contents of a text file.
Start line and stop line are zero-based.  Start is inclusive,
stop is exclusive.

Usage:  textfile.sh <file> <start line> <stop line>

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

function tf() {
	local CMD="java $EA $EOOM -Xmx120m -cp $CP fileIO.TextFile $@"
	echo $CMD >&2
	eval $CMD
}

tf "$@"
