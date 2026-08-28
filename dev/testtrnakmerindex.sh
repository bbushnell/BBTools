#!/bin/bash

usage(){
echo "
Written by G11
Last modified August 2026

Description:  Runs TrnaKmerIndexBoundTest -- focused regression coverage for the 2026-08-27
unique-query-k-mer dedup fix in TrnaKmerIndex.shortlist() and its qKmers denominator (unique
query-k-mer count, not the raw positional count). Self-contained assertion-based test, no args.

Usage:  testtrnakmerindex.sh
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
	#DIR must be reassigned to the BBTools ROOT (not dev/) BEFORE sourcing javasetup.sh: javasetup.sh
	#reads the variable literally named DIR (its own "if [ -n "$DIR" ]; then . "$DIR/memdetect.sh""
	#check) -- a "../" prefix only on the source LINE below would still leave DIR itself pointing at
	#dev/, so javasetup.sh's internal memdetect.sh sourcing would still silently fail (caught non-
	#fatally, exits 0 anyway -- the same bug dev/testhbm.sh has; see the TODO left there, since fixing
	#that file is out of this change's scope).
	DIR="$DIR""../"
	. "$DIR""javasetup.sh"

	parseJavaArgs "--xmx=1g" "--xms=1g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prok.TrnaKmerIndexBoundTest $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
