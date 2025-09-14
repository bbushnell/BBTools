#!/bin/bash

usage(){
echo "
Description:  Generates a histogram of kmer counts for the input reads or assemblies.

Usage:	khist.sh in=<input> hist=<histogram output>

Please see bbnorm.sh for more information.  
All the flags are the same, only the parameters (near the bottom of this file) differ.
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

z="-Xmx31g"
z2="-Xms31g"
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
	freeRam 31000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

khist() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP jgi.KmerNormalize bits=32 ecc=f passes=1 keepall dr=f prefilter hist=stdout minprob=0 minqual=0 mindepth=0 minkmers=1 hashes=3 $@"
	echo $CMD >&2
	eval $CMD
}

khist "$@"
