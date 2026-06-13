#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and UMP45
Last modified June 13, 2026

Description:  Restores SEQ and QUAL on secondary (0x100) and supplementary (0x800)
alignments by copying them from the read's primary (non-supplementary) alignment,
matched by read name.  Aligners like minimap2 emit SEQ=* on secondary/supplementary
records, which makes them unusable for variant calling.  Restoring from the primary
recovers the TRUE bases and qualities with no reference bias (unlike MD-tag restore).

The input is partitioned by name-hash into temporary headerless sam.gz subfiles so all
alignments of a read share a subfile; each subfile is then name-sorted and restored
single-threaded.  The output header is rewritten with SO:unsorted.

Usage:  restorebases.sh in=<file> out=<file>

Parameters and their defaults:

in=<file>       Input sam or bam file (with secondary/supplementary alignments).
out=<file>      Output sam or bam file.
ways=31         Number of temporary subfiles.  Increase for very large inputs to
                reduce the memory needed to hold one subfile at a time.
                ways=1 skips partitioning and writes no temp files (loads the whole
                input into memory); best for bacterial-sized BAMs.
mode=fix        fix: write the bases/quals into SEQ/QUAL.
                tag: attach OS:Z: (bases) and OQ:Z: (quals) tags in original
                     sequencing orientation, leaving SEQ=*.
                'fix' and 'tag' may also be given as bare flags (no 'mode=').
hardclip=convert  How to handle supplementary hard-clips in fix mode:
                convert: change H->S and copy the full read (default).
                truncate: keep the H cigar and copy only the supplementary's own
                     segment.  Smaller output; gives identical variant calls.
tmpdir=<dir>    Directory for the temporary subfiles (default: out + '_rbtmp').
                Temp files and the directory are always deleted on completion.
reads=-1        Process only this many input alignments (-1 = all); for testing.
ow=t            (overwrite) Overwrite existing output files.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
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
	SCRIPT="$0"
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

	parseJavaArgs "--xmx=4000m" "--xms=4000m" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP var2.RestoreBases $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
