#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 27, 2017

Description:  Prints sequence gc content once per interval.

Usage:  plotgc.sh in=<input file> out=<output file>

Parameters:
in=<file>       Input file. in=stdin.fa will pipe from stdin.
out=<file>      Output file.  out=stdout will pipe to stdout.
interval=1000   Interval length.
offset=0        Position offset.  For 1-based indexing use offset=1.
psb=t           (printshortbins) Print gc content for the last bin of a contig
                even when shorter than interval.

Java Parameters:

-Xmx            This will set Java's memory usage, overriding automatic
                memory detection. -Xmx20g will 
                specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  
                The max is typically 85% of physical memory.
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
JNI="-Djava.library.path=""$DIR""jni/"
JNI=""

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
	freeRam 1400m 42
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

plotgc() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP driver.PlotGC $@"
	echo $CMD >&2
	eval $CMD
}

plotgc "$@"
