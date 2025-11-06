#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Isla
Last modified November 6, 2025

Description:  Runs any BBTools Java class with Java Flight Recorder profiling.

Usage:  profile.sh <classname> <arguments> profile=<output.jfr>
e.g.
profile.sh stream.StreamerWrapper in=foo.sam profile=profile.jfr
profile.sh align2.BBMap in=reads.fq ref=genome.fa profile=mapping.jfr -Xmx8g

Parameters:
profile=<file>  Output JFR file (required).
<classname>     Fully qualified Java class to run (first non-flag argument).
-Xmx<size>      Java heap size (optional, default 2g).

All other parameters are passed to the target class.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
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

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

CP="$DIR""current/"

calcXmx () {
    # Source the new scripts
    source "$DIR""/memdetect.sh"
    source "$DIR""/javasetup.sh"
    
    parseJavaArgs "--mem=2g" "--mode=fixed" "$@"
    
    # Set environment paths
    setEnvironment
}
calcXmx "$@"

# Extract profile output file
PROFILE="profile.jfr"
CLASSNAME=""
MAXSIZE="2g"
ARGS=()

for arg in "$@"; do
	if [[ $arg == profile=* ]]; then
		PROFILE="${arg#profile=}"
	elif [[ -z "$CLASSNAME" ]] && [[ ! $arg == -* ]] && [[ ! $arg == *=* ]]; then
		# First non-flag, non-key=value argument is the class name
		CLASSNAME="$arg"
	elif [[ $arg == maxsize=* ]]; then
		MAXSIZE="${arg#maxsize=}"
	else
		# Pass everything else to the target class
		ARGS+=("$arg")
	fi
done

if [[ -z "$PROFILE" ]]; then
	echo "Error: profile=<filename.jfr> parameter required" >&2
	exit 1
fi

if [[ -z "$CLASSNAME" ]]; then
	echo "Error: No class name specified" >&2
	exit 1
fi

profiler() {
	local JFR_OPTS="-Xlog:jfr*=off -XX:StartFlightRecording:settings=profile,jdk.CPUTimeSample#enabled=true,filename=$PROFILE,maxsize=$MAXSIZE,dumponexit=true"
	local CMD="java $EA $EOOM $XMX $XMS $SIMD $JFR_OPTS -cp $CP $CLASSNAME ${ARGS[@]}"
	echo "$CMD" >&2
	eval "$CMD"
}
profiler
