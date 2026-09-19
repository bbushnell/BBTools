#!/bin/bash

usage(){
echo "
BBMapS — BBMap alignment with the Streamer/Writer interface.

Single-ended:  bbmapS.sh ref=reference.fa in=reads.fq out=mapped.sam
Paired-end:    bbmapS.sh ref=reference.fa in=R1.fq in2=R2.fq out=mapped.sam
Index only:    bbmapS.sh ref=reference.fa path=index
Reuse index:   bbmapS.sh in=reads.fq out=mapped.sam path=index
Split reads:   bbsplitS.sh ref_a=a.fa ref_b=b.fa in=reads.fq basename=out_%.fq

Same flag surface as bbmap.sh (build=, in=, in2=, ref=, t=, out=, etc.).
Optional short-indel acceleration: quantumonebase=t (default f).
Optional no-MSA speed mode: quantumonly=t (default f; changes scoring/mapping).
Optional k-mer pseudoalignment: pseudoalign=t (default f; mapped SAM uses CIGAR=*;
  no base alignment or identity/edit filters; intended for coverage/counting).
Run bbmap.sh -h for the full flag reference.
Java SIMD is detected by the standard BBTools launcher setup.
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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP align2.BBMapS build=1 overwrite=true fastareadlen=500 $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" align2.BBMapS build=1 overwrite=true fastareadlen=500 "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
