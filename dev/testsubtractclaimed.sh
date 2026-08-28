#!/bin/bash

usage(){
echo "
Written by G11
Last modified August 2026

Description:  Runs TrnaCallerSubtractClaimedTest -- focused regression coverage for the
2026-08-27 subtractClaimed fix (interior-claim splits must keep BOTH remainders, not just
the right one). Self-contained assertion-based test, no args.

Usage:  testsubtractclaimed.sh
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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prok.TrnaCallerSubtractClaimedTest $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
