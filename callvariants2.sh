#!/bin/bash

#Running callvariants2.sh is equivalent to running callvariants.sh with the "multi" flag.
#See callvariants.sh for usage information.

#callvariants2 is intended for multiple sam/bam files, one from each sample, which should have variants called independently; the point is that allele frequencies will be reported for ALL samples at locations where ANY sample has a variant called.
#callvariants2 is NOT a better version of callvariants, it's the same, just designed for multisample processing.
#If you have only 1 sample (regardless of how many sam/bam files there are) you should use callvariants.sh without the "multi" flag.

usage(){
bash "$DIR"callvariants.sh
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

callvariants2() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP var2.CallVariants2 $@"
	echo $CMD >&2
	eval $CMD
}

callvariants2 "$@"
