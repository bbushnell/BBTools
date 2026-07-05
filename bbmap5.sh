#!/bin/bash

usage(){
echo "
BBMap5
Written by Brian Bushnell, from Dec. 2010 - present
Last modified June 25, 2026

Description:  Fast and accurate splice-aware read aligner.
This is the BBMap5 variant, which packs reference coordinates as 32-bit
UNSIGNED values.  Its purpose is large chromosomes: it can index and map
references with single chromosomes up to ~1.07 Gbp (1,073,525,813 bp),
such as plant genomes like wheat, whose largest chromosomes exceed the
~536 Mbp single-chromosome limit of standard BBMap.

For most references, standard bbmap.sh and bbmap5.sh produce equivalent
results, and which is faster should be determined empirically.  But for a
reference with very long chromosomes, BBMap5 is currently the only option.

All parameters are identical to bbmap.sh.  For the full list of options and
documentation, run bbmap.sh (with no arguments) or read
bbmap/docs/guides/BBMapGuide.txt.

Usage:                      bbmap5.sh ref=<fasta> in=<reads> out=<sam>
Index only:                 bbmap5.sh ref=<fasta>
Map to existing index:      bbmap5.sh in=<reads> out=<sam>
Map without writing index:  bbmap5.sh ref=<fasta> in=<reads> out=<sam> nodisk

Java Parameters:
-Xmx                    This will set Java's memory usage, overriding
                        autodetection.  -Xmx20g will specify 20 gigs of RAM.
                        The index uses roughly 6 bytes per reference base, so
                        large genomes (e.g. wheat) need substantial memory.
-eoom                   Exit if an out-of-memory exception occurs.
-da                     Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter
any problems, or post at: http://seqanswers.com/forums/showthread.php?t=41057
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

	parseJavaArgs "--xmx=3200m" "--xms=3200m" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP align2.BBMap5 build=1 overwrite=true fastareadlen=500 $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
