#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified April 5, 2020

Description:  Summarizes coverage information from basecov files
              created by pileup.sh.  They should be named like
              'sample1_basecov.txt' but other naming styles are fine too.

Usage:        summarizecoverage.sh *basecov.txt out=<output file>

Parameters:
in=<file>           'in=' is not necessary.  Any filename used as a
                    parameter will be assumed to be an input basecov file.
out=<file>          Write the summary here.  Default is stdout.
reflen=-1           If positive, use this as the total reference length.
                    Otherwise, assume basecov files report every ref base.

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

z="-Xmx200m"
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

summarize() {
	local CMD="java $EA $EOOM $z -cp $CP covid.SummarizeCoverage $@"
#	echo $CMD >&2
	eval $CMD
}

summarize "$@"
