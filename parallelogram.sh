#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 4, 2025

Description:  Converts a parallelogram-shaped alignment visualization to a rectangle.
This tool transforms the output from CrossCutAligner so it can be properly
visualized by visualizealignment.sh. The transformation shifts coordinates
to create a rectangular matrix from the parallelogram pattern.

Usage:
parallelogram.sh <input_map> <output_map>

Parameters:
input_map       Input text file containing parallelogram-shaped matrix data.
output_map      Output text file with rectangular matrix data.

Example workflow:
crosscutaligner.sh ATCGATCG GCATGCTA map1.txt
parallelogram.sh map1.txt map2.txt
visualizealignment.sh map2.txt alignment.png

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

    parseJavaArgs "--mem=200m" "--mode=fixed" "$@"

    # Set environment paths
    setEnvironment
}
calcXmx "$@"

transform() {
	if [[ $# -eq 0 ]]; then
		usage
		return
	fi
	local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP aligner.Parallelogram $@"
	#echo $CMD >&2
	eval $CMD
}

transform "$@"