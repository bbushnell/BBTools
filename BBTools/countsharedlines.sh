#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified September 15, 2015

Description:  Counts the number of lines shared between sets of files.
One output file will be printed for each input file.  For example,
an output file for a file in the 'in1' set will contain one line per
file in the 'in2' set, indicating how many lines are shared.

Usage:  countsharedlines.sh in1=<file,file...> in2=<file,file...>

Parameters:
include=f       Set to 'true' to include the filtered names rather than excluding them.
prefix=f        Allow matching of only the line's prefix (all characters up to first whitespace).
case=t          (casesensitive) Match case also.
ow=t            (overwrite) Overwrites files that already exist.
app=f           (append) Append to files that already exist.
zl=4            (ziplevel) Set compression level, 1 (low) to 9 (max).

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
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

z="-Xmx800m"
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
	freeRam 800m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

function countsharedlines() {
	local CMD="java $EA $EOOM $z -cp $CP driver.CountSharedLines $@"
	echo $CMD >&2
	eval $CMD
}

countsharedlines "$@"
