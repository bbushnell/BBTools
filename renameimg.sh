#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 23, 2017

Description:  Renames img records to be prefixed by their id.
This is for internal JGI use and has no external utility.

Usage:  renameimg.sh in=auto out=renamed.fa.gz

Parameters:
in=         3-column tsv with imgID, taxID, and file path.
            These files will have their sequences renamed and concatenated.
img=        Optional, if a different (presumably bigger) file will be used for taxonomic assignment.
            For example, in could be a subset of img, potentially with incorrect taxIDs.

Java Parameters:
-Xmx        This will set Java's memory usage, overriding autodetection.
            -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom       This flag will cause the process to exit if an out-of-memory exception occurs.  Requires Java 8u92+.
-da         Disable assertions.

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

z="-Xmx1g"
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

function rename() {
	local CMD="java $EA $EOOM $z -cp $CP tax.RenameIMG $@"
	echo $CMD >&2
	eval $CMD
}

rename "$@"
