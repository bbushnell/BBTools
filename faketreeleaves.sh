#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 7, 2026

Description:  Assigns real taxonomy parents for permanent FakeTree leaf nodes.

Usage:  faketreeleaves.sh in=fake_leaf_genomes.tsv ssu=findssu_top_hits.tsv \\
        quickclade=quickclade_top_hits.tsv out=proposals.tsv \\
        outnodes=fake_leaf_nodes.tsv

Parameters:
in=<file>          Fake leaf genome map; required.
ssu=<file>         FindSSU top-hit table; required.
quickclade=<file>  QuickClade top-hit table; required.
out=<file>         Output proposal table; required.
outnodes=<file>    Output fake leaf node table; required.
outinspect=<file>  Optional subset needing manual inspection.
summary=<file>     Optional summary counts.
tree=auto          TaxTree path.
ow=t              Overwrite existing output files.

Java Parameters:
-Xmx               This will set Java's memory usage, overriding autodetection.
-eoom              Exit on out-of-memory.  Requires Java 8u92+.
-da                Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
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

	parseJavaArgs "--xmx=6g" "--xms=6g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP tax.FakeTreeLeafTaxonomy $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" tax.FakeTreeLeafTaxonomy "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
