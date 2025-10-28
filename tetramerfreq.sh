#!/bin/bash

usage(){
echo "
Written by Shijie Yao and Brian Bushnell
Last modified April 25, 2025

Description: DNA Tetramer analysis.
DNA tetramers are counted for each sub-sequence of window size in the sequence.  
The window slides along the sequence by the step length.
Sub-sequence shorter than the window size is ignored. Tetramers containing N are ignored. 

Usage: tetramerfreq.sh in=<input file> out=<output file> step=500 window=2000

Input may be fasta or fastq, compressed or uncompressed.

Standard parameters:
in=<file>       DNA sequence input file 
out=<file>      Output file name
step/s=INT      Step size (default 500) 
window/w=INT    Window size (default 2kb); <=0 turns windowing off (e.g. short reads)
short=T/F       Print lines for sequences shorter than window (default F)
k=INT           Kmer length (default 4)
gc              Print a GC column in the output.
float           Output kmer frequencies instead of counts.
comp            Output GC-compensated kmer frequencies.

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

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
    # Source the new scripts
    source "$DIR""/memdetect.sh"
    source "$DIR""/javasetup.sh"
    
    # Parse Java arguments with tool-specific defaults
    # Use auto mode with 42% of available RAM, minimum 1000MB
    parseJavaArgs "--mem=1000m" "--percent=42" "--mode=auto" "$@"
    
    # Set environment paths
    setEnvironment
    
    # Set the Java memory parameters
    z="-Xmx${RAM}m"
    z2="-Xms${RAM}m"
}
calcXmx "$@"

tetramerfreq() {
	local CMD="java $EA $SIMD $EOOM $XMX -cp $CP jgi.TetramerFrequencies $@"
	echo $CMD >&2
	eval $CMD
}

tetramerfreq "$@"
