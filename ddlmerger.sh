#!/bin/bash

usage(){
echo "
Written by Ady
Last modified April 18, 2026

Description:  Loads multiple DDL files, merges records sharing a TID
(e.g., combining mitochondrial/plastid/plasmid with host genomes),
sorts, renumbers, and writes a single combined DDL file.

Source categories are detected from filenames:
  *mito* = mitochondrial, *plastid* = plastid, *plasmid* = plasmid

Usage:  ddlmerger.sh in=a.ddl.gz,b.ddl.gz out=combined.ddl.gz

Parameters:
in=<file>           Input DDL file(s), comma-delimited.
                    Can also be given as bare arguments.
out=<file>          Output combined DDL file.
k=31                K-mer length (for validation).
merge=t             Merge all categories with matching TIDs.
mergemito=t         Merge mitochondrial records into host genome.
mergeplastid=t      Merge plastid records into host genome.
mergeplasmid=t      Merge plasmid records into host genome.
overwrite=f         Overwrite existing output.
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

	parseJavaArgs "--xmx=4g" "--xms=400m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ddl.DDLMerger $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
