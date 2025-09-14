#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 3, 2017

Description:  Creates polymer sequences.
Can be used in conjunction with mutate.sh to generate low-complexity sequence.

Usage:  makepolymers.sh out=<output file> k=<repeat length> minlen=<sequence length>

I/O parameters:
out=<file>      Output genome.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Processing parameters:
k=1             Length of repeating polymeric units.
                To generate a sweep of multiple values of k,
                specify both mink and maxk.
minlen=31       Ensure sequences are at least this long.
                Specifically, minlen=X will ensure sequences are long enough
                that all possible kmers of length X are present.

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
	freeRam 4000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

makepolymers() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP jgi.MakePolymers $@"
	echo $CMD >&2
	eval $CMD
}

makepolymers "$@"
