#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified Oct 11, 2023

Description:  Scores sequences using a neural network.  Only the initial Xbp
are used, for sequences longer than the network size.

Usage:  scoresequence.sh in=<sequences> out=<renamed sequences> net=<net file>

Input may be fasta or fastq, compressed or uncompressed.

Standard parameters:
in=<file>       Sequence data.
out=<file>      Sequences renamed with their scores.
net=<file>      Network file to apply to the sequences.
hist=<file>     Histogram of scores (x100, so 0-1 maps to 0-100).
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Processing parameters:
rcomp=f         Use the max score of a sequence and its reverse complement.
parse=f         Parse sequence headers for 'result=' to determine whether
                they are positive or negative examples.
annotate=t      Rename output reads by appending 'score='.
filter=f        Retain only reads above or below a cutoff.  Setting the cutoff
                or highpass flag will automatically set this to true.
cutoff=0.5      Score cutoff for filtering; scores mostly range from 0 to 1.
highpass=t      Retain sequences ABOVE cutoff if true, else BELOW cutoff.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

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
z2="-Xms1g"
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

score() {
	local CMD="java $EA $EOOM $z $SIMD -cp $CP ml.ScoreSequence $@"
	echo $CMD >&2
	eval $CMD
}

score "$@"
