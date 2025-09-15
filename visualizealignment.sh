#!/bin/bash

usage(){
echo "
Shell script written by Brian Bushnell
Java code written by Claude.
Last modified May 4, 2025

Description:  Converts a text exploration map from some aligners to an image.
Supports Quantum, Banded, Drifting, Glocal, WaveFront, and MSA9. 

Usage:
visualizealignment.sh <map>
or
visualizealignment.sh <map> <image>

Parameters:
map             Text file of score-space from an aligner.
image           Output name, context sensitive; supports png, bmp, jpg.
                Image name is optional; if absent, .txt will be replaced
                by .png in the input filename.

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
    
    # Parse Java arguments with tool-specific defaults
    # Use auto mode with 84% of available RAM, minimum 4000MB
    parseJavaArgs "--mem=2000m" "--percent=42" "--mode=auto" "$@"
    
    # Set environment paths
    setEnvironment
    
    # Set the Java memory parameters
    z="-Xmx${RAM}m"
    z2="-Xms${RAM}m"
}
calcXmx "$@"

convert() {
	local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP aligner.VisualizationConverter $@"
	#echo $CMD >&2
	eval $CMD
}

convert "$@"
