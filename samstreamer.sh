#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 29, 2025

Description:  Converts sam/bam to sam, bam, fasta, or fastq rapidly.
Multithreaded.

Usage:  streamsam.sh in=<file> out=<file>

Filtering parameters:
minpos=         Ignore alignments not overlapping this range.
maxpos=         Ignore alignments not overlapping this range.
minmapq=        Ignore alignments with mapq below this.
maxmapq=        Ignore alignments with mapq above this.
contigs=        Comma-delimited list of contig names to include. These 
                should have no spaces, or underscores instead of spaces.
mapped=t        Include mapped reads.
unmapped=t      Include unmapped reads.
secondary=f     Include secondary alignments.
supplimentary=t Include supplimentary alignments.
lengthzero=f    Include alignments without bases.
invert=f        Invert sam filters.
ordered=t       Keep reads in input order.
ref=<file>      Optional reference file.
simd            Add this flag for turbo speed.  Requires Java 17+ and AVX2.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
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
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

calcXmx () {
    # Source the new scripts
    source "$DIR""/memdetect.sh"
    source "$DIR""/javasetup.sh"
    
    parseJavaArgs "--mem=2g" "--mode=fixed" "$@"
    
    # Set environment paths
    setEnvironment
}
calcXmx "$@"

samstreamer() {
	local CMD="java $EA $EOOM $z $SIMD -cp $CP stream.SamStreamerWrapper $@"
	echo $CMD >&2
	eval $CMD
}

samstreamer "$@"
