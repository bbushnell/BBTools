#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 21, 2015

Description:    Merges coverage stats lines (from pileup) for the same OTU,
                according to some custom naming scheme.

Usage:          mergeOTUs.sh in=<file> out=<file>

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

function mergeOTUs() {
	local CMD="java $EA $EOOM $z -cp $CP driver.MergeCoverageOTU $@"
	echo $CMD >&2
	eval $CMD
}

mergeOTUs "$@"
