#!/bin/bash

usage(){
echo "
Description:  Tests HBM (BaseGraph) scoring for tRNA candidates.
Builds HBM models for one anticodon, scores known tRNAs and random
genome fragments to compare identity vs HBM score separation.

Usage:  testhbm.sh trnas=all_trnas.fa anticodon=GUC genome=ecoli.fna.gz
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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prok.TestHBM $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
