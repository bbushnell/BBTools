#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 4, 2019

Description:  Discards all but the best copy of a ribosomal gene per TaxID.
Gene sequences should be named like this: >tid|123|whatever
Sequences are selected based on the number of fully defined bases.

Usage:  keepbestcopy.sh in=<input file> out=<output file> rate=<float>

Input may be fasta or fastq, compressed or uncompressed.

Standard parameters:
in=<file>       Input sequences.
out=<file>      Output sequences.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression
                level; lower compression is faster.

Processing parameters:
maxlen=1600     Prefer sequences shorter than this.

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

z="-Xmx8g"
z2="-Xms8g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 8g 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

keepbest() {
	local CMD="java $EA $SIMD $EOOM $z -cp $CP jgi.KeepBestCopy $@"
	echo $CMD >&2
	eval $CMD
}

keepbest "$@"
