#!/bin/bash

usage(){
echo "
Description:  Maps taxIDs to phylum names using the BBTools taxonomy tree.
Reads taxIDs from stdin (one per line), writes taxID<tab>phylum to stdout.

Usage:  cat taxids.txt | taxidtophylum.sh > phylum_map.tsv
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
	. "$DIR/../javasetup.sh"
	. "$DIR/../memdetect.sh"

	parseJavaArgs "--xmx=4g" "--xms=4g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prok.TaxIDToPhylum $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
