#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified December 6, 2018

Description:  Gathers statistics on Kapa spike-in rates.

Usage:  kapastats.sh in=<input file> out=<output file>

Parameters:
in=<file>       TSV file of plate IDs, one ID per line.
out=<file>      Primary output, or read 1 output.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
raw=f           Output raw observations rather than statistics.

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

z="-Xmx200m"
z2="-Xms200m"
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

kapastats() {
	local CMD="java $EA $EOOM $z -cp $CP jgi.GatherKapaStats $@"
	echo $CMD >&2
	eval $CMD
}

kapastats "$@"
