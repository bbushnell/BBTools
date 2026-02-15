#!/bin/bash

usage(){
echo "
Written by Brian Bushnell.
Last modified December 19, 2019

Description:  Renames sequences to indicate their NCBI taxIDs.
The headers must be in NCBI or Silva format with gi numbers,
accessions, or organism names.  Only supports fasta and gff files.

Usage:  gi2taxid.sh in=<file> out=<file> server

Parameters:
in=<file>       Input sequences; required parameter.  Must be fasta.
                This can alternatively be a comma-delimited list,
                or just a bunch of space-delimited filenames, e.g.:
                gi2taxid.sh x.fa y.fa z.fa out=tid.fa tree=auto table=auto
out=<file>      Destination for renamed sequences.
invalid=<file>  Destination for headers with no taxid.
keepall=t       Keep sequences with no taxid in normal output.
prefix=t        Append the taxid as a prefix to the old header, but keep
                the old header.
title=tid       Set the title of the new number (e.g. ncbi, taxid, tid).
ziplevel=2      (zl) Compression level for gzip output.
pigz=t          Spawn a pigz (parallel gzip) process for faster 
                compression than Java.  Requires pigz to be installed.
silva=f         Parse headers in Silva format.
shrinknames=f   Replace multiple concatenated headers with the first.
deleteinvalid=f Delete the output file if there are any invalid headers.

Taxonomy File Parameters:
server=f        Use the taxonomy server instead of local files.
                Server mode only works for accessions (like RefSeq).
tree=           Specify a taxtree file.  On Genepool, use 'auto'.
gi=             Specify a gitable file.  On Genepool, use 'auto'.
accession=      Specify one or more comma-delimited NCBI accession to
                taxid files.  On Genepool, use 'auto'.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx800m will specify 800 megs.
                    The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

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
	CP="$DIR/current/"
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=7g" "--xms=7g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP tax.RenameGiToTaxid $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
