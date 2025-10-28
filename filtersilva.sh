#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified Jan 27, 2020

Description:  Removes .

Usage:           filtersilva.sh in=x.fa out=y.fa

Standard parameters:
in=<file>       Input file.
out=<file>      Output file.

Additional file parameters:
tree=auto       Path to TaxTree.


Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-da             Disable assertions.

For more detailed information, please read /bbmap/docs/guides/BBSketchGuide.txt.
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

z="-Xmx4g"
z2="-Xms4g"
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
	#freeRam 3200m 84
	#z="-Xmx${RAM}m"
	#z2="-Xms${RAM}m"
}
calcXmx "$@"

filtersilva() {
	local CMD="java $EA $SIMD $EOOM $z -cp $CP prok.FilterSilva $@"
	echo $CMD >&2
	eval $CMD
}

filtersilva "$@"
