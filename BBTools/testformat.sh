#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 4, 2016

Description:  Tests file extensions and contents to determine format,
quality, compression, interleaving, and read length.  More than one file
may be specified.  Note that ASCII-33 (sanger) and ASCII-64 
(old Illumina/Solexa) cannot always be differentiated.

Usage:  testformat.sh <file>

See also:  testformat2.sh, stats.sh

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

z="-Xmx120m"
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

testformat() {
	local CMD="java $EA $EOOM $z -cp $CP fileIO.FileFormat $@"
#	echo $CMD >&2
	eval $CMD
}

testformat "$@"
