#!/bin/bash

usage(){
echo "
Written by Chloe
Last modified October 18, 2025

Description:  Converts BAM (Binary Alignment/Map) files to SAM 
(Sequence Alignment/Map) text format. Reads BGZF-compressed BAM files 
and outputs tab-delimited SAM format.

Usage:  bamlinestreamer.sh <input.bam> <output.sam>

Standard parameters:
in=<file>        Input BAM file (first positional argument).
out=<file>       Output SAM file (second positional argument).

Java Parameters:
-Xmx             This will set Java's memory usage, overriding autodetection.
                 -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                 specify 200 megs. The max is typically 85% of physical memory.
-eoom            This flag will cause the process to exit if an out-of-memory
                 exception occurs.  Requires Java 8u92+.
-da              Disable assertions.

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

calcXmx () {
    # Source the new scripts
    source "$DIR""/memdetect.sh"
    source "$DIR""/javasetup.sh"

    parseJavaArgs "--mem=8g" "--mode=fixed" "$@"

    # Set environment paths
    setEnvironment
}
calcXmx "$@"

bamlinestreamer() {
	if [[ $# -eq 0 ]]; then
		usage
		return
	fi
	local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP stream.bam.Bam2Sam $@"
	#echo $CMD >&2
	eval $CMD
}

bamlinestreamer "$@"
