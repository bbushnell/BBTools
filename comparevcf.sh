#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified April 30, 2019

Description:  Performs set operations on VCF files:
Union, intersection, and subtraction.

Usage:  comparevcf.sh in=<file,file,...> out=<file>

I/O parameters:
in=<file>       Input; must be at least 2 files.
out=<file>      Output file.
ref=<file>      Reference file; optional.  Usually not needed.
shist=<file>    (scorehist) Output for variant score histogram.
overwrite=f     (ow) Set to false to force the program to abort rather than
bgzip=f         Use bgzip for gzip compression.

Mode Parameters (choose one only):
subtract=t      Subtract all other files from the first file.
union=f         Make a union of all files.
intersection=f  Make an intersection of all files.

Processing Parameters:
addsamples=t    Include all samples in the output lines. (TODO)
splitalleles=f  Split multi-allelic lines into multiple lines.
splitsubs=f     Split multi-base substitutions into SNPs.
canonize=t      Trim variations down to a canonical representation.

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

z="-Xmx4g"
z2="-Xms4g"
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
	freeRam 4000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

comparevcf() {
	local CMD="java $EA $SIMD $EOOM $z $z2 -cp $CP var2.CompareVCF $@"
	echo $CMD >&2
	eval $CMD
}

comparevcf "$@"
