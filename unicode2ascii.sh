#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 17, 2017

Description:  Replaces unicode and control characters with printable ascii characters.
WARNING - this does not work in many cases, and is not recommended!
It is only retained because there is some situation in which it is needed.

Usage:  unicode2ascii.sh in=<file> out=<file>

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

function unicode2ascii() {
	local CMD="java $EA $EOOM $z -cp $CP jgi.UnicodeToAscii $@"
	echo $CMD >&2
	eval $CMD
}

unicode2ascii "$@"
