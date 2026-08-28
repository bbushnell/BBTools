#!/bin/bash

set -e

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

runTest(){
	local class="$1"
	java -ea --add-modules jdk.incubator.vector -cp "$CP" "$class"
}

resolveSymlinks
runTest assemble.BubblePopperUnitTest
runTest assemble.PathPreservingBubbleSimplifierSpec
runTest assemble.ReadThreadedXResolverUnitTest
runTest assemble.CrossKTipOverlapperUnitTest
