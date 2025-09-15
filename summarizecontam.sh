#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified March 19, 2018

Description:  Summarizes monthly contam files into a single file.
This is for internal JGI use.

Usage:  summarizecontam.sh <input files> out=<output file>

Parameters:
in=<file,file>  Input contam summary files, comma-delimited.
                Alternately, file arguments (from a * expansion) will be 
                considered input files.
out=<file>      Output.
tree=auto       Taxtree file location (optional).
overwrite=t     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Filter Parameters (passing all required to pass):
minreads=0      Ignore records with fewer reads than this.
minsequnits=0   Ignore records with fewer seq units than this.

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
	freeRam 1000m 24
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

process() {
	local CMD="java $EA $EOOM $z -cp $CP driver.SummarizeContamReport $@"
	echo $CMD >&2
	eval $CMD
}

process "$@"
