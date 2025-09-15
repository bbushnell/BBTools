#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 31, 2025

Description:  Tests multiple aligners using random sequences.
The sequences have variable pairwise ANI, and each
ANI level is tested multiple times for average accuracy
and loop count.
Outputs the identity, rstart and rstop positions, time, and #loops.
Note that the 'design' ANI is approximate and will not match
the measured ANI.

Usage:
testaligners2.sh iterations=30 maxani=100 minani=90 step=2

Parameters:
length=40k      Length of sequences.
iterations=32   Iterations to average; higher is more accurate.
maxani=80       Max ANI to model.
minani=30       Min ANI to model.
step=2          ANI step size.
sinewaves=0     Sinewave count to model variable conservation.
threads=        Parallel alignments; default is logical cores.
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
    
    parseJavaArgs "--mem=3200m" "--mode=fixed" "$@"
    
    # Set environment paths
    setEnvironment
}
calcXmx "$@"

align() {
	local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP aligner.TestAlignerSuite $@"
	#echo $CMD >&2
	eval $CMD
}

align "$@"
