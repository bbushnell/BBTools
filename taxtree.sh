#!/bin/bash

usage(){
echo "
Written by Brian Bushnell.
Last modified Jan 7, 2020

Description:  Creates tree.taxtree from names.dmp and nodes.dmp.
These are in taxdmp.zip available at ftp://ftp.ncbi.nih.gov/pub/taxonomy/
The taxtree file is needed for programs that can deal with taxonomy, 
like Seal and SortByTaxa.

Usage:  taxtree.sh names.dmp nodes.dmp merged.dmp tree.taxtree.gz

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM.  The max is typically 85% of physical memory.
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

z="-Xmx2g"
z2="-Xms2g"
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
	freeRam 2000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"


taxtree() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP tax.TaxTree $@"
	echo $CMD >&2
	eval $CMD
}

taxtree "$@"
