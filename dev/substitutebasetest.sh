#!/bin/bash

usage(){
echo "
Written by Nepgear
Last modified August 2026

Description:  Exhaustive verification for ukmer.Kmer#substituteBase (Item 3c, Fix #3).
For several deterministic sequences (mixed bases, poly-T, a palindromic/canonical-flip-
prone repeat, a realistic mixed sequence), every position x every alternate base, checks
substituteBase() against an independently-rebuilt Kmer (never derived from the in-place
object), then checks the restore is bit-identical to the untouched original. Runs under
both PACKED=true (k=32,33,63,64,95,127) and a representative PACKED=false (legacy
symmetric layout) subset (kbig=62,93), since the primitive derives its geometry from the
live object's k/perWordK/kbig fields rather than assuming a fixed word width.

This is a committed regression harness (current/ukmer/SubstituteBaseTest.java), pending
Brian's decision on whether it becomes permanent -- not a scratch-only file. Exits
non-zero on any failure (the test itself throws), so the shell exit status is
authoritative.

Usage:  substitutebasetest.sh

Expected output: 'Total checks: 9104' / 'Failures: 0' / 'ALL PASS', exit code 0.

Contact: bbushnell@lbl.gov
"
}

if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
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
	DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=200m" "--xms=200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ukmer.SubstituteBaseTest"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
