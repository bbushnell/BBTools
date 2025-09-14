#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 14, 2024

Description:  Processes a tile dump from FilterByTile.

Usage:  tiledump.sh in=<input file> out=<output file>

Standard parameters:
in=<file>       Input dump file.
out=<file>      Output dump file.
overwrite=t     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Processing parameters:
x=-1            Widen tiles to at least this X width.
y=-1            Widen tiles to at least this Y width.
reads=-1        Widen tiles to at least this average number of reads.
alignedreads=250  Average aligned reads per tile for error rate calibration.

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
	freeRam 4000m 48
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

tiledump() {
	local CMD="java $EA $EOOM $z -cp $CP hiseq.TileDump $@"
	echo $CMD >&2
	eval $CMD
}

tiledump "$@"
