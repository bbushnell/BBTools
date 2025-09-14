#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 3, 2014

Description:  Transforms interleaved phylip to fasta.

Usage:   phylip2fasta.sh in=<input> out=<output>

Input may be stdin or an interleaved phylip file, compressed or uncompressed.

Input Parameters:
in=<file>       The input phylip file; this is the only required parameter.
unpigz=true     Decompress with pigz for faster decompression.

Output Parameters:
out=<file>      Fasta output destination.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
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

z="-Xmx1g"
z2="-Xms1g"
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
	freeRam 800m 82
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

convert() {
	local CMD="java $EA $EOOM $z -cp $CP jgi.PhylipToFasta $@"
	echo $CMD >&2
	eval $CMD
}

convert "$@"
