#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified April 25, 2019

Description:  Compresses or decompresses files based on extensions.
This only exists because the syntax and default behavior of many
compression utilities is unintuitive; it is just a wrapper, and
relies on existing executables in the command line (pigz, lbzip, etc.)
Does not delete the input file.
Does not untar files.

Usage:  unzip.sh in=<file> out=<file>

Parameters:
in=<file>       Input file.
out=<file>      Output file for good reads.
zl=             Set the compression level; 0-9 or 11.

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

z="-Xmx80m"
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

unzip() {
	local CMD="java $EA $EOOM $z -cp $CP jgi.Unzip $@"
#	echo $CMD >&2
	eval $CMD
}

unzip "$@"
