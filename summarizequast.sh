#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 17, 2016

Description:  Summarizes the output of multiple Quast reports for
making box plots.

Usage:  summarizequast.sh */quast/report.tsv

Parameters:
out=stdout      Destination for summary.
required=       A required substring in assembly names for filtering.
normalize=t     Normalize each metric to the average per report.
box=t           Print only 5 points per metric for box plots.

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

summarizequast() {
	local CMD="java $EA $EOOM $z -cp $CP driver.SummarizeQuast $@"
#	echo $CMD >&2
	eval $CMD
}

summarizequast "$@"
