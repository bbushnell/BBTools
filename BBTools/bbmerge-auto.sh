#!/bin/bash

usage(){
echo "
bbmerge-auto.sh is a wrapper for BBMerge that attempts to use all available
memory, instead of a fixed amount.  This is for use with the Tadpole options
of error-correction (ecct) and extension, which require more memory.
For merging by overlap only, please use bbmerge.sh.  If you set memory
manually with the -Xmx flag, bbmerge.sh and bbmerge-auto.sh are equivalent.

For information about usage and parameters, please run bbmerge.sh.
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

z="-Xmx14g"
z2="-Xms14g"
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
	freeRam 15000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

function merge() {
	local CMD="java $EA $EOOM $z $z2 $JNI -cp $CP jgi.BBMerge $@"
	echo $CMD >&2
	eval $CMD
}

merge "$@"
