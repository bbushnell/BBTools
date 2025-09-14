#!/bin/bash

usage(){
echo "
Written by Jasper Toscani Field
Last modified June 11, 2020

Description:  Counts positional occurrences of reference kmers in reads.

Usage:  kmerposition.sh in=<input file> out=<output file> ref=<reference file> k=<kmer length>

Input may be fasta or fastq, compressed or uncompressed.

Standard parameters:
in=<file>       Primary input, or read 1 input.
in2=<file>      Read 2 input if reads are in two files.
ref=<file>      Reference file.
out=<file>      Primary output, statistics on found kmers.

Processing parameters:
k=19            Kmer length.
rcomp=t         If true, also match for reverse-complements.

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

z="-Xmx200m"
z2="-Xms200m"
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

kmerposition3() {
	local CMD="java $EA $EOOM $z -cp $CP jasper.KmerPosition3 $@"
	echo $CMD >&2
	eval $CMD
}

kmerposition3 "$@"
