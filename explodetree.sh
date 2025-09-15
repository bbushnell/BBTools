#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified December 13, 2017

Description:   Constructs a directory and file tree of sequences
corresponding to a taxonomic tree.

Usage:  explodetree.sh in=<file> out=<path> tree=<file>

Parameters:
in=             A fasta file annotated with taxonomic data in headers, 
                such as modified RefSeq.
out=            (path) Location to write the tree.
tree=           Location of taxtree file.

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

z="-Xmx2000m"
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

function explodetree() {
	local CMD="java $EA $EOOM $z -cp $CP tax.ExplodeTree $@"
	echo $CMD >&2
	eval $CMD
}

explodetree "$@"
