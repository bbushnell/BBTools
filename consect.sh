#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 25, 2016

Description:  Generates the conservative consensus of multiple
error-correction tools.  Corrections will be accepted only
if all tools agree.  This tool is designed for substitutions only,
not indel corrections.

Usage:  consect.sh in=<file,file,file,...> out=<file>

Standard parameters:
in=             A comma-delimited list of files; minimum of 3.
                All files must have reads in the same order.
                The first file must contain the uncorrected reads.
                All additional files must contain corrected reads.
out=<file>      Output of consensus reads.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Processing parameters:
None yet!

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
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

z="-Xmx2g"
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

consect() {
	local CMD="java $EA $EOOM $z -cp $CP jgi.Consect $@"
	echo $CMD >&2
	eval $CMD
}

consect "$@"
