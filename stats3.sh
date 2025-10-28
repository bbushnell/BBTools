#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 21, 2025

Description:  In progress.
Generates some assembly stats for multiple files.

Usage:        stats3.sh in=file
Or:           stats3.sh in=file,file
Or:           stats3.sh file file file

Parameters:
in=file         Specify the input fasta file(s), or stdin.
                Multiple files can be listed without a 'in=' flag.
out=stdout      Destination of primary output; may be directed to a file.

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

stats() {
	local CMD="java $EA $SIMD $EOOM $z -cp $CP jgi.AssemblyStats3 $@"
#	echo $CMD >&2
	eval $CMD
}

stats "$@"
