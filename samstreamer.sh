#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 30, 2025

Description:  Converts sam/bam to sam, bam, fasta, or fastq rapidly.
Allows some filtering operations.  Multithreaded.

Usage:  samstreamer.sh in=<file> out=<file> <other arguments>
or
samstreamer.sh <input_file> <output_file> <other arguments>
e.g.
samstreamer.sh mapped.bam mapped.sam.gz minid=0.95 simd


File parameters
in=             Input file, must be sam, sam.gz, or bam.
                To stream, use in=stdin.sam with the proper extension.
out=            Output file, optional, any sequence format allowed.

Filtering parameters:
minpos=         Ignore alignments not overlapping this range.
maxpos=         Ignore alignments not overlapping this range.
minmapq=        Ignore alignments with mapq below this.
maxmapq=        Ignore alignments with mapq above this.
minid=0.0       Ignore alignments with identity below this.
maxid=1.0       Ignore alignments with identity above this.
contigs=        Comma-delimited list of contig names to include. These 
                should have no spaces, or underscores instead of spaces.
mapped=t        Include mapped reads.
unmapped=t      Include unmapped reads.
secondary=t     Include secondary alignments.
                Recommended false for fasta/fastq.
supplimentary=t Include supplimentary alignments.
                Recommended false for fasta/fastq.
lengthzero=t    Include alignments without bases.
                Recommended false for fasta/fastq.
invert=f        Invert sam filters.
ordered=t       Keep reads in input order.
ref=<file>      Optional reference file.
simd            Add this flag for turbo speed.  Requires Java 17+ and AVX2,
                or other 256-bit vector instruction sets.

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
    
    parseJavaArgs "--mem=1g" "--mode=fixed" "$@"
    
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
