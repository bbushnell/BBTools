#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 20, 2020

Description:  Estimates cardinality of unique kmers in sequence data.
Processes multiple kmer lengths simultaneously to produce a histogram.

Usage:  kmercountmulti.sh in=<file> sweep=<20,100,8> out=<histogram output>

Parameters:
in=<file>           (in1) Input file, or comma-delimited list of files.
in2=<file>          Optional second file for paired reads.
out=<file>          Histogram output.  Default is stdout.
k=                  Comma-delimited list of kmer lengths to use.
sweep=min,max,incr  Use incremented kmer values from min to max. For example,
                    sweep=20,26,2 is equivalent to k=20,22,24,26.
buckets=2048        Use this many buckets for counting; higher decreases
                    variance, for large datasets.  Must be a power of 2.
seed=-1             Use this seed for hash functions.  
                    A negative number forces a random seed.
minprob=0           Set to a value between 0 and 1 to exclude kmers with a 
                    lower probability of being correct.
hashes=1            Use this many hash functions.  More hashes yield greater
                    accuracy, but H hashes takes H times as long.
stdev=f             Print standard deviations.

Shortcuts:
The # symbol will be substituted for 1 and 2.
For example:
kmercountmulti.sh in=read#.fq
...is equivalent to:
kmercountmulti.sh in1=read1.fq in2=read2.fq

Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.

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

z="-Xmx500m"
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

function kmercountmulti() {
	local CMD="java $EA $EOOM $z -cp $CP jgi.KmerCountMulti $@"
	echo $CMD >&2
	eval $CMD
}

kmercountmulti "$@"
