#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified March 28, 2018

Description:  Finds and trims junctions in mapped Hi-C reads.
For the purpose of reporting junction motifs, this requires paired-end reads,
because only improper pairs will be considered as possibly containing 
junctions.  However, all reads that map with soft-clipping will be trimmed
on the 3' (right) end, regardless of pairing status.

Usage:  processhi-c.sh in=<mapped reads> out=<trimmed reads>

Parameters:
in=<file>       A sam/bam file containing mapped reads.
out=<file>      Output file of trimmed reads.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
printkmers=t    Generate files with kmer counts at junction sites.
junctions=junctions_k%.txt  File pattern for junction output.

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

process() {
	local CMD="java $EA $EOOM $z -cp $CP jgi.FindHiCJunctions $@"
	echo $CMD >&2
	eval $CMD
}

process "$@"
