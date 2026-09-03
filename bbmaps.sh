#!/bin/bash

usage(){
echo "
BBMapS — EXPERIMENTAL Streamer/Writer-based BBMap (bbmapnova Phase 1)
First slice: t=1 plumbing only. RouteWriter/BBSplitterInvoker are still
NoOp placeholders, so no real SAM/FASTQ output is written yet — this is a
smoke-test driver for the dispatcher/worker/coordinator lifecycle, not a
replacement for bbmap.sh. See align2/BBMapS.java and
/mnt/c/playground/Nowi/plans/BBMapUpgrade_Phase1_CoordinatorDesign_Nowi.md.

Usage:  bbmaps.sh ref=<fasta> in=<reads> out=<sam> t=1

Same flag surface as bbmap.sh where implemented (build=, in=, in2=, ref=,
t=, out= — out= is accepted but not yet honored: NoOp writers discard it).
Run bbmap.sh -h for the full flag reference; unimplemented flags here are
silently inherited from AbstractMapper's parser but most have no effect
until Yaoyao's real RouteWriter/BBSplitterInvoker land.
"
}

if [ -z "$1" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	usage
	exit
fi

resolveSymlinks(){
	SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=3200m" "--xms=3200m" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP align2.BBMapS build=1 overwrite=true fastareadlen=500 $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
