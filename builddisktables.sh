#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Chloe
Last modified October 10, 2026

Description:   Builds binary disk-backed taxonomy tables for TaxServer.
The accession table is a memory-mapped hash table used by disk-backed
TaxServer mode.  The optional GI table is for legacy gi number support.

Usage:  builddisktables.sh accession=<files> pattern=<file> accout=<file> gi=<file> giout=<file>

Usage examples:
builddisktables.sh accession=auto pattern=auto accout=accession_disk.bin
builddisktables.sh accession=auto pattern=auto accout=accession_disk.bin gi=auto giout=gi_disk.bin

Parameters:
accession=       Comma-delimited NCBI accession-to-taxid files.
                 Use 'auto' on Dori when the default shrunk files are present.
pattern=         Accession pattern table produced by analyzeaccession.sh.
                 Use 'auto' on Dori when the default pattern table is present.
accout=          Output path for the binary disk accession table.
gi=              Optional input GI table, such as gitable.int2d.gz.
                 This may also be set with table= or gitable=.
giout=           Optional output path for the binary disk GI table.
tree=auto        Taxonomy tree.
taxpath=auto     Set the path to taxonomy files; auto only works at NERSC.

Java Parameters:
-Xmx             This will set Java's memory usage, overriding autodetection.
                 -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                 specify 200 megs.  The max is typically 85% of physical memory.
-eoom            This flag will cause the process to exit if an out-of-memory
                 exception occurs.  Requires Java 8u92+.
-da              Disable assertions.

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

	parseJavaArgs "--xmx=24g" "--xms=24g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP tax.BuildDiskTables $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" tax.BuildDiskTables "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
