#!/bin/bash

usage(){
echo "
Written by Citan
Last modified September 2026

Description:  Runs CallGenesSixsConfigTest -- focused configuration checks for the
default-off paired 6S ncRNA bundle (sixs=/ssrs=/6s=, subordinate to ncrna=). Verifies
default-off state, family-name aliasing, per-family idPass/idBorderline/hbmPass defaults,
boundary offsets, gate subordination, and override-requires-gate behavior, plus a
resource-backed constructed-family check (testShippedBundle()) that requires the six
staged sixs_rf00013_*/sixs_rf01685_* resources and fails loudly if they are absent.
Self-contained assertion-based test, no args.

Usage:  testsixsconfig.sh
"
}

pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

resolveSymlinks(){
	if [ -d "$DIR""current/" ]; then
		CP="$DIR""current/"
	elif [ -d "$DIR""../current/" ]; then
		CP="$DIR""../current/"
	fi
}

setEnv(){
	#DIR must be reassigned to the BBTools ROOT before sourcing javasetup.sh -- see the fix note in
	#testtrnakmerindex.sh (javasetup.sh reads DIR itself to find memdetect.sh).
	DIR="$DIR""../"
	. "$DIR""javasetup.sh"

	parseJavaArgs "--xmx=1g" "--xms=1g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prok.CallGenesSixsConfigTest $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
