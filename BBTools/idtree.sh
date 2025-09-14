#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified July 1, 2016

Description:  Makes a Newick tree from an identity matrix.
Intended for use with matrices created by idmatrix.sh.

Usage:  idtree.sh in=<input file> out=<output file>

Standard parameters:
in=<file>       Identity matrix in TSV format.
out=<file>      Newick tree output.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

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

z="-Xmx800m"
z2="-Xms800m"
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
	freeRam 800m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

idtree() {
	local CMD="java $EA $EOOM $z -cp $CP tax.IDTree $@"
	echo $CMD >&2
	eval $CMD
}

idtree "$@"
