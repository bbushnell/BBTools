#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified December 19, 2018

Description:  Parses a webcheck log.
Input is expected to look like this:
Tue Apr 26 16:40:09 2016|https://rqc.jgi-psf.org/control/|200 OK|0.61

Usage:  webcheck.sh <input files>


Standard parameters:
in=<file>       Primary input.  Can use a wildcard (*) if 'in=' is omitted.
out=<file>      Summary output; optional.
fail=<file>     Output of failing lines; optional.
invalid=<file>  Output of misformatted lines; optional.
extendedstats=f (es) Print more stats.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

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

webcheck() {
	local CMD="java $EA $EOOM $z -cp $CP driver.ProcessWebcheck $@"
	echo $CMD >&2
	eval $CMD
}

webcheck "$@"
