#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified March 5, 2018

Description:  Benchmarks a disk with multithreaded I/O.

Usage:  diskbench.sh path=<path> data=<8g> passes=<2> threads=<>

Parameters:
path=           Location to read and write.
data=8g         Number of bytes to process per pass.
threads=        Number of threads to use.  By default, all logical threads.
                In RW mode the number of active threads is doubled.
mode=rw         I/O mode:
                   r:  Test read speed only.
                   w:  Test write speed only.
                   rw: Test read and write speed simultaneously.

Processing parameters:
None yet!

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

z="-Xmx400m"
z2="-Xms400m"
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

diskbench() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP fun.DiskBench $@"
	echo $CMD >&2
	eval $CMD
}

diskbench "$@"
