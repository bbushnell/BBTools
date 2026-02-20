#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 17, 2015

Description:  Separates paired reads into files of 'good' pairs and 'good' singletons by removing 'bad' reads that are shorter than a min length.
Designed to handle situations where reads become too short to be useful after trimming.  This program also optionally performs quality trimming.

Usage:        bbsplitpairs.sh in=<input file> out=<pair output file> outs=<singleton output file> minlen=<minimum read length, an integer>

Input may be fasta or fastq, compressed or uncompressed.

Optional parameters (and their defaults)

in=<file>           The 'in=' flag is needed if the input file is not the first parameter.  'in=stdin' will pipe from standard in.
in2=<file>          Use this if 2nd read of pairs are in a different file.
out=<file>          The 'out=' flag is needed if the output file is not the second parameter.  'out=stdout' will pipe to standard out.
out2=<file>         Use this to write 2nd read of pairs to a different file.
outsingle=<file>    (outs) Write singleton reads here.

overwrite=t         (ow) Set to false to force the program to abort rather than overwrite an existing file.
showspeed=t         (ss) Set to 'f' to suppress display of processing speed.
interleaved=auto    (int) If true, forces fastq input to be paired and interleaved.
qtrim=f             Trim read ends to remove bases with quality below trimq.
                    Values: rl (trim both ends), f (neither end), r (right end only), l (left end only).
trimq=6             Trim quality threshold.
minlen=20           (ml) Reads shorter than this after trimming will be discarded.
ziplevel=2          (zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster.
fixinterleaving=f   (fint) Fixes corrupted interleaved files by examining pair names.  Only use on files with broken interleaving.
repair=f            (rp) Fixes arbitrarily corrupted paired reads by examining read names.  High memory.
ain=f               (allowidenticalnames) When detecting pair names, allows identical names, instead of requiring /1 and /2 or 1: and 2:

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

	parseJavaArgs "--xmx=200m" "--xms=200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.SplitPairsAndSingles $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
