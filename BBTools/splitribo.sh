#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 22, 2020

Description:  Splits a file of various rRNAs into one file per type
(16S, 18S, 5S, 23s).

Usage:  splitribo.sh in=<file,file> out=<pattern>

Standard parameters:
in=<file>       Input file.
out=<pattern>   Output file pattern, such as out_#.fa.  The # symbol is
                required and will be substituted by the type name, such as
                16S, to make out_16S.fa, for example.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
ziplevel=9      (zl) Set to 1 (lowest) through 9 (max) to change compression
                level; lower compression is faster.

types=16S,18S,5S,23S,m16S,m18S,p16S
                Align to these sequences.  Fewer types is faster.  m16S
                and m18S are mitochondrial; p16S is plastid (chloroplast).

Processing parameters:
minid=0.59      Ignore alignments with identity lower than this to a 
                consensus sequences.
refineid=0.70   Refine score by aligning to clade-specific consensus if
                the best alignment to a universal consensus is below this.

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
	freeRam 4000m 42
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

mergeribo() {
	local CMD="java $EA $EOOM $z -cp $CP prok.SplitRibo $@"
	echo $CMD >&2
	eval $CMD
}

mergeribo "$@"
