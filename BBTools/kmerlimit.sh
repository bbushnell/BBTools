#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified July 31, 2018

Description:  Stops producing reads when the unique kmer limit is reached.
This is approximate.  If the input has been Clumpified, the order should be
randomized first with shuffle2.sh or sortbyname.sh with the flowcell flag.

Differences between versions:
kmerlimit.sh uses 1 pass and outputs all reads until a limit is hit,
meaning the input reads should be in random order with respect to sequence.
kmerlimit2.sh uses 2 passes and randomly subsamples from the file, so
it works with reads in any order.

Usage:  kmerlimit.sh in=<input file> out=<output file> limit=<number>

Standard parameters:
in=<file>       Primary input, or read 1 input.
in2=<file>      Read 2 input if reads are in two files.
out=<file>      Primary output, or read 1 output.
out2=<file>     Read 2 output if reads are in two files.
overwrite=t     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression
                level; lower compression is faster.

Processing parameters:
k=31            Kmer length, 1-32.
limit=          The number of unique kmers to produce.
mincount=1      Ignore kmers seen fewer than this many times.
minqual=0       Ignore bases with quality below this.
minprob=0.2     Ignore kmers with correctness probability below this.

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

z="-Xmx1000m"
z2="-Xms1000m"
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

kmerlimit() {
	local CMD="java $EA $EOOM $z -cp $CP sketch.KmerLimit $@"
	echo $CMD >&2
	eval $CMD
}

kmerlimit "$@"
