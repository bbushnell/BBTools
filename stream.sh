#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 13, 2025

Description:  Converts between sam, bam, fasta, fastq.
              Supports subsampling, paired files, and multithreading.

Usage:  stream.sh in=<file> out=<file> <other arguments>
or
stream.sh <input_file> <output_file> <other arguments>
e.g.
stream.sh mapped.bam mapped.sam.gz
stream.sh in=reads.fq out=subset.fq samplerate=0.1

File parameters:
in=<file>       Primary input file, type detected from extension.
in2=<file>      Secondary input file for paired reads.
out=<file>      Primary output file, optional, type based on extension.
out2=<file>     Secondary output file for paired reads.
                Note: Use # symbol for auto-numbering, e.g. reads_#.fq

Processing parameters:
samplerate=1.0  Fraction of reads to keep (0.0 to 1.0).
sampleseed=17   Random seed for subsampling (-1 for random).
reads=-1        Quit after processing this many reads (-1 = all).
ordered=t       Maintain input order in output.

Threading parameters:
threadsin=-1    Reader threads (-1 = auto).
threadsout=-1   Writer threads (-1 = auto).

Other parameters:
simd            Add this flag for turbo speed. Requires Java 17+ and AVX2,
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

streamer() {
	local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP stream.StreamerWrapper $@"
	echo $CMD >&2
	eval $CMD
}

streamer "$@"