#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 12, 2019

Description:  Compares CDS, rRNA, and tRNA lines in gff files.

Usage:  comparegff.sh in=<input gff> ref=<reference gff>

Standard parameters:
in=<file>       Query gff.
ref=<file>      Reference gff.

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

z="-Xmx1g"
z2="-Xms1g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 1000m 42
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

comparegff() {
	local CMD="java $EA $EOOM $z -cp $CP gff.CompareGff $@"
	echo $CMD >&2
	eval $CMD
}

comparegff "$@"
