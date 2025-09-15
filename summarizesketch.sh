#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 25, 2018

Description:  Summarizes the output of BBSketch. 

Usage:  summarizesketch.sh in=<file,file...> out=<file>

You can alternately run 'summarizesketch.sh *.txt out=out.txt'

Parameters:
in=<file>       A list of stats files, or a text file containing one stats file name per line.
out=<file>      Destination for summary.
tree=           A TaxTree file.
level=genus     Ignore contaminants with the same taxonomy as the primary hit at this level.
unique=f        Use the contaminant with the most unique hits rather than highest score.

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

summarizesketch() {
	local CMD="java $EA $EOOM $z -cp $CP sketch.SummarizeSketchStats $@"
#	echo $CMD >&2
	eval $CMD
}

summarizesketch "$@"
