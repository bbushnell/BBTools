#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified Oct 30, 2023

Description:  Scores sequences using a neural network.  This is similar
to scoresequence.sh but multithreaded and with more filtering options.

Usage:  netfilter.sh in=<sequences> out=<pass> outu=<fail> net=<net file>

Input may be fasta or fastq, compressed or uncompressed.

Standard parameters:
in=<file>       Input sequences.
out=<file>      Sequences passing the filter.
outu=<file>     Sequences failing the filter.
net=<file>      Network file to apply to the sequences.
hist=<file>     Histogram of scores (x100, so 0-1 maps to 0-100).
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Processing parameters:
rcomp=f         Use the max score of a sequence and its reverse complement.
parse=f         Parse sequence headers for 'result=' to determine whether
                they are positive or negative examples.  Positive examples
                should be annotated with result=1, and negative with result=0.
annotate=f      Rename output sequences by appending 'score='.
filter=t        Retain only reads above or below a cutoff.  Setting the cutoff
                or highpass flag will automatically set this to true.
cutoff=auto     Score cutoff for filtering; scores mostly range from 0 to 1.
                'auto' will use the cutoff embedded in the network.
highpass=t      Retain sequences ABOVE cutoff if true, else BELOW cutoff.
scoremode=      single (default): Apply the network once to each sequence.
                  If the sequence is longer than the network's inputs, use the
                  first X bases.
                average: Apply the network to the sequence once every 
                  'stepsize' bases, and use the average.  May be slow.
                max: Like average, but uses the maximum score.
                min: Like average, but uses the minimum score.
pairmode=       average (default): For paired reads, average the two scores.
                max: Use the higher of the two scores.
                min: Use the lower of the two scores.
stepsize=1      If scoremode is other than 'single', score a window every
                this many bases.  The window width is defined by the network.
                Higher values of stepsize are faster.
overlap=        This can be set instead of stepsize; if either flag is used it
                will override the other.  Setting overlap will make windows
                overlap that much, so 'overlap=0' is equivalent to 'stepsize=W'
                where W is the width of the network in bases.

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

z="-Xmx2g"
z2="-Xms2g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 2000m 20
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

score() {
	local CMD="java $EA $EOOM $z $z2 $SIMD -cp $CP ml.NetFilter $@"
	echo $CMD >&2
	eval $CMD
}

score "$@"
