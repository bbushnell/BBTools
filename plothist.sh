#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 5, 2024

Description:  Generates histograms from a tile dump.
Also works on other 2D numeric matrices with a header.
Output files are automatically named from the header columns.

Usage:  plothist.sh in=<input file> bins=<number>

Parameters:
in=<file>       Input dump file.
bins=1000       Bins per histogram.
overwrite=t     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

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
z2="-Xms300m"
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

plothist() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP hiseq.PlotHist $@"
	echo $CMD >&2
	eval $CMD
}

plothist "$@"
