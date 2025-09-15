#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 30, 2025

Description:  Statistical analysis tool that calculates the Average Nucleotide
Identity (ANI) between random DNA sequences. Generates pairs of random 
sequences of specified lengths, aligns them, and produces a histogram of 
identity distributions. This demonstrates that random sequences converge to
varying, length-dependent identity approaching roughly approximately 55%,
with standard deviation decreasing with length, providing a baseline for
evaluating the significance of real sequence alignments.

Usage:  alignrandom.sh <start> <mult> <steps> <iters> <buckets> <maxloops> <output>

Positional Parameters and Defaults (optional, ordered, without the name):
start=10        Starting sequence length for analysis
mult=10         Length multiplier between intervals (each step: length*=mult)
steps=4         Number of length intervals to test 
iters=200       Number of random sequence pairs to align per interval
buckets=100     Number of histogram bins for identity distribution
maxloops=max    Maximum total alignments to prevent excessive runtime
output=stdout   Output file for ANI histogram results

Example:
alignrandom.sh 20 5 6 500
Tests lengths 20, 100, 500, 2500, 12500, 62500 with 500 iterations each.

Output:
Produces a tab-delimited histogram showing the distribution of alignment
identities for random sequence pairs at each tested length.

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
    
    parseJavaArgs "--mem=8g" "--mode=fixed" "$@" simd
    
    # Set environment paths
    setEnvironment
}
calcXmx "$@" simd

align() {
	local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP aligner.AlignRandom $@"
	#echo $CMD >&2
	eval $CMD
}

align "$@"
