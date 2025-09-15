#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 16, 2025

Description:  Calculates EST (expressed sequence tags) capture by an assembly from a sam file.
Designed to use BBMap output generated with these flags: k=13 maxindel=100000 customtag ordered

Usage:          bbest.sh in=<sam file> out=<stats file>

Parameters:
in=<file>       Specify a sam file (or stdin) containing mapped ests.
out=<file>      Specify the output stats file (default is stdout).
ref=<file>      Specify the reference file (optional).
est=<file>      Specify the est fasta file (optional).
fraction=0.98   Min fraction of bases mapped to ref to be 
                considered 'all mapped'.

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

z="-Xmx120m"
z2="-Xms120m"
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

function bbest() {
	local CMD="java $EA $EOOM $z -cp $CP jgi.SamToEst $@"
#	echo $CMD >&2
	eval $CMD
}

bbest "$@"
