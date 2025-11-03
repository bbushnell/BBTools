#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified Novemeber 2, 2025

Description:  Converts bteween sam, bam, fasta, fastq.

Usage:  stream.sh in=<file> out=<file> <other arguments>
or
stream.sh <input_file> <output_file> <other arguments>
e.g.
stream.sh mapped.bam mapped.sam.gz

File parameters
in=             Input file, must have correct extension.
                To stream, use e.g. in=stdin.sam.
out=            Output file, optional, type based on extension.
                To stream, use e.g. out=stdout.fq
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
	local CMD="java $EA $EOOM $z $SIMD -cp $CP stream.FastqWriter $@"
	echo $CMD >&2
	eval $CMD
}

streamer "$@"
