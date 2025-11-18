#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified April 24, 2019

Description:  Reformats a fungal assembly for release.
Also creates contig and agp files.

Usage:  fungalrelease.sh in=<input file> out=<output file>

I/O parameters:
in=<file>           Input scaffolds.
out=<file>          Output scaffolds.
outc=<file>         Output contigs.
qfin=<file>         Optional quality scores input.
qfout=<file>        Optional quality scores output.
qfoutc=<file>       Optional contig quality scores output.
agp=<file>          Output AGP file.
legend=<file>       Output name legend file.
overwrite=f         (ow) Set to false to force the program to abort rather than
                    overwrite an existing file.

Processing parameters:
fastawrap=60        Wrap length for fasta lines.
tuc=t               Convert sequence to upper case.
baniupac=t          Crash on encountering a non-ACGTN base call.
mingap=10           Expand all gaps (Ns) to be at least this long.
mingapin=1          Only expand gaps that are at least this long.
sortcscaffolds=t    Sort scaffolds descending by length.
sortcontigs=f       Sort contigs descending by length.
renamescaffolds=t   Rename scaffolds to 'scaffold_#'.
scafnum=1           Number of first scaffold.
renamecontigs=f     Rename contigs to 'contig_#' instead of 'scafname_c#'.
contignum=1         Number of first contig; only used if renamecontigs=t.
minscaf=1           Only retain scaffolds at least this long.
mincontig=1         Only retain contigs at least this long.

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
	SCRIPT="$0"
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

	parseJavaArgs "--xmx=4g" "--xms=4g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.FungalRelease $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
