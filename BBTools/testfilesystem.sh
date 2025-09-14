#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified December 11, 2017

Description:  Logs filesystem performance by creating, deleting,
and copying files.

Usage:  testfilesystem.sh <in> <out> <log> <size> <ways> <interval in seconds>

'in' should contain the # symbol if ways>1.

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

z="-Xmx50m"
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

function testfs() {
	local CMD="java $EA $EOOM $z -cp $CP jgi.TestFilesystem $@"
	echo $CMD >&2
	eval $CMD
}

testfs "$@"
