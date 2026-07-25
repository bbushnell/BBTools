#!/bin/bash

usage(){
echo "
Written by Eru with Brian Bushnell
Last modified July 24, 2026

Description:  Renames CAMI benchmark contigs to include tid_ labels.
Reads a CAMI binning_gs.tsv key file and appends _tid_TAXID to contig
names in FASTA and/or SAM/BAM files, making them compatible with
BBTools tools that use tid_ labels (QuickBin, GradeBins, etc.).

Usage:  renamecami.sh in=contigs.fa out=renamed.fa key=binning_gs.tsv
or
renamecami.sh in=contigs.fa sam=mapped.bam key=binning_gs.tsv out=renamed.fa outsam=renamed.bam

Parameters:
in=<file>       Input FASTA assembly.
out=<file>      Output renamed FASTA.
sam=<file>      Input SAM/BAM alignment file.
outsam=<file>   Output renamed SAM/BAM.
key=<file>      CAMI binning_gs.tsv file mapping contig names to taxIDs.
                Format: tab-separated, columns SEQUENCEID BINID TAXID.
                Lines starting with @ or # are skipped.
drop=f          Drop contigs/alignments not found in the key file.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
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
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=4g" "--xms=4g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP bin.RenameCAMI $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
