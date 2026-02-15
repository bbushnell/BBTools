#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified September 18, 2019

Description:  Runs TadpoleWrapper after some preprocessing,
to allow optimal assemblies using long kmers.
Only paired reads are supported.

Usage:
tadpipe.sh in=reads.fq out=contigs.fa


Parameters:
in=<file>           Input reads.
in2=<file>          Optional read 2, if reads are in two files.
out=contigs.fa      Output file name.
temp=$TMPDIR        Path to a directory for temp files.
delete=t            Delete intermediate files.
gz=f                Gzip intermediate files.

Other parameters can be passed to individual phases like this:

assemble_k=200,250  Set kmer lengths for assembly phase.
merge_strict        Set the strict flag in merge phase.
extend_el=120       Set the left-extension distance in the extension phase.

Valid prefixes:

filter_             PhiX and contaminant filtering.
trim_               Adapter trimmming.
merge_              Paired-read merging.
correct_            Error correction.
extend_             Read extension.
assemble_           Final assembly.
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

	parseJavaArgs "--xmx=14g" "--xms=14g" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP assemble.TadPipe $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
