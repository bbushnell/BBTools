#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 10, 2018

Description:  Merges .pgm files.

Usage:  mergepgm.sh in=x.pgm,y.pgm out=z.pgm

File parameters:
in=<file,file>  A pgm file or comma-delimited list of pgm files.
out=<file>      Output filename.
normalize=f     Merge proportionally to base counts, so small models
                have equal weight to large models.  Normalization happens
                before applying the @ multiplier.
@ symbol        Input filenames in the form of 'x.pgm@0.1' will have
                a multiplier applied to that model prior to merging.

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

function mergepgm() {
	local CMD="java $EA $EOOM $z -cp $CP prok.PGMTools $@"
	echo $CMD >&2
	eval $CMD
}

mergepgm "$@"
