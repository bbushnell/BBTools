#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 22, 2024

Description:  Filters a sam file to select only reads with substitution errors
for bases with quality scores in a certain interval.  Used for manually 
examining specific reads that may have incorrectly calibrated quality scores.

Usage:  filtersubs.sh in=<file> out=<file> minq=<number> maxq=<number>

Parameters:
in=<file>       Input sam or bam file.
out=<file>      Output file.
minq=0          Keep only reads with substitutions of at least this quality.
maxq=99         Keep only reads with substitutions of at most this quality.
countindels=t   Also keep reads with indels in the quality range. 
minsubs=1       Require at least this many substitutions.
minclips=0      Discard reads with more clip operations than this.
maxclips=-1     If nonnegative, discard reads with more clip operations.
keepperfect=f   Also keep error-free reads.

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

filtersubs() {
	local CMD="java $EA $EOOM $z -cp $CP jgi.FilterReadsWithSubs $@"
#	echo $CMD >&2
	eval $CMD
}

filtersubs "$@"
