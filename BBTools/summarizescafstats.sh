#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 26, 2015

Description:  Summarizes the scafstats output of BBMap for evaluation
of cross-contamination.  The intended use is to map multiple libraries or 
assemblies, of different multiplexed organisms, to a concatenated reference 
containing one fused scaffold per organism.  This will convert all of the 
resulting stats files (one per library) to a single text file, with multiple 
columns, indicating how much of the input hit the primary versus nonprimary 
scaffolds.

Usage:  summarizescafstats.sh in=<file,file...> out=<file>

You can alternatively use a wildcard, like this:
summarizescafstats.sh scafstats_*.txt out=summary.txt

Parameters:
in=<file>       A list of stats files, or a text file containing one stats file name per line.
out=<file>      Destination for summary.

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

summarizescafstats() {
	local CMD="java $EA $EOOM $z -cp $CP driver.SummarizeCoverage $@"
#	echo $CMD >&2
	eval $CMD
}

summarizescafstats "$@"
