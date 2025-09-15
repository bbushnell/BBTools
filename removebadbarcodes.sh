#!/bin/bash

usage(){
echo "
Written by Brian Bushnell.
Last modified March 5, 2018

Description:  Removes reads with barcodes containing non-ACGT bases.
Read headers must be in standard Illumina format.

Usage:  removebadbarcodes.sh in=<file> out=<file>

Parameters:
in=<file>       Input reads; required parameter.
out=<file>      Destination for good reads; optional.
ziplevel=2      (zl) Compression level for gzip output.
pigz=f          Spawn a pigz (parallel gzip) process for faster 
                compression than Java.  Requires pigz to be installed.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx800m will specify 800 megs.
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
JNI="-Djava.library.path=""$DIR""jni/"
JNI=""

z="-Xmx200m"
z2="-Xms200m"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
}
calcXmx "$@"


removebadbarcodes() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP jgi.RemoveBadBarcodes $@"
	echo $CMD >&2
	eval $CMD
}

removebadbarcodes "$@"
