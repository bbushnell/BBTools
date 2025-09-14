#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 18, 2025

Description:  Aligns a query sequence to a reference using multiple aligners.
Outputs the identity, rstart and rstop positions, time, and #loops.

Usage:
testaligners.sh <query> <ref>
testaligners.sh <query> <ref> <iterations> <threads> <simd>

Parameters:
query           A literal nucleotide sequence or fasta file.
ref             A literal nucleotide sequence or fasta file.
iterations      Optional integer for benchmarking multiple iterations.
threads         Number of parallel instances to use.
simd            Enable SIMD operations; requires AVX-256 and Java 17+.

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

calcXmx () {
    # Source the new scripts
    source "$DIR""/memdetect.sh"
    source "$DIR""/javasetup.sh"
    
    parseJavaArgs "--mem=2000m" "--mode=fixed" "$@"
    
    # Set environment paths
    setEnvironment
}
calcXmx "$@"

align() {
	local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP aligner.Test $@"
	#echo $CMD >&2
	eval $CMD
}

align "$@"
