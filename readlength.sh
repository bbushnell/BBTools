#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified December 19, 2018
Description:  Generates a length histogram of input reads.

Usage:	readlength.sh in=<input file>

Parameters:
in=<file>    	The 'in=' flag is needed only if the input file is not the first parameter.  'in=stdin.fq' will pipe from standard in.
in2=<file>   	Use this if 2nd read of pairs are in a different file.
out=<file>   	Write the histogram to this file.  Default is stdout.
bin=10       	Set the histogram bin size.
max=80000    	Set the max read length to track.
round=f      	Places reads in the closest bin, rather than the highest bin of at least readlength.
nzo=t        	(nonzeroonly) Do not print empty bins.
reads=-1     	If nonnegative, stop after this many reads.

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
z2="-Xms400m"
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
	local CMD="java $EA $SIMD $EOOM $z -cp $CP jgi.MakeLengthHistogram $@"
#	echo $CMD >&2
	eval $CMD
}

stats "$@"
