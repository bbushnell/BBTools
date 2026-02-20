#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified September 12, 2018

Description:  Sorts reads by name or other keys such as length,
quality, mapping position, flowcell coordinates, or taxonomy.
Intended to merge temp files produced by SortByName if the program
ran out of time during merging.

Usage:   mergesorted.sh sort_temp* out=<file>

Input may be fasta, fastq, or sam, compressed or uncompressed.

Parameters:

in=<file,file,...>  Input files.  Files may be specified without in=.
out=<file>          Output file.
delete=t            Delete input files after merging.
name=t              Sort reads by name.
length=f            Sort reads by length.
quality=f           Sort reads by quality.
position=f          Sort reads by position (for mapped reads).
taxa=f              Sort reads by taxonomy (for NCBI naming convention).
sequence=f          Sort reads by sequence, alphabetically.
flowcell=f          Sort reads by flowcell coordinates.
shuffle=f           Shuffle reads randomly (untested).
list=<file>         Sort reads according to this list of names.
ascending=t         Sort ascending.
memmult=.35         Write a temp file when used memory drops below this
                    fraction of total memory.

Taxonomy-sorting parameters:
tree=               Specify a taxtree file.  On Genepool, use 'auto'.
gi=                 Specify a gitable file.  On Genepool, use 'auto'.
accession=          Specify one or more comma-delimited NCBI accession to
                    taxid files.  On Genepool, use 'auto'.

Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.

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

	parseJavaArgs "--xmx=2g" "--xms=2g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP sort.MergeSorted $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
