#!/bin/bash

usage(){
echo "
Last modified February 13, 2020

Description:  Wrapper for BBMap to allow multiple input and output files for the same reference without reloading the index each time.

Usage:  bbwrap.sh ref=<reference fasta> in=<file,file,...> out=<file,file,...> nodisk
To index only:                bbwrap.sh ref=<reference fasta>
To map to an existing index:  bbwrap.sh in=<file,file,...> out=<file,file,...>
To map pairs and singletons and output them into the same file:
bbwrap.sh in1=read1.fq,singleton.fq in2=read2.fq,null out=mapped.sam append

BBWrap will not work with stdin and stdout, or histogram output.

Other Parameters:

in=<file,file>  Input sequences to map.
inlist=<fofn>   Alternately, input and output can be a file of filenames,
                one line per file, using the flag inlist, outlist, outmlist,
                in2list, etc.
mapper=bbmap    Select mapper.  May be BBMap, BBMapPacBio,
                or BBMapPacBioSkimmer.
append=f        Append to files rather than overwriting them.
                If append is enabled, and there is exactly one output file,
                all output will be written to that file.

***** All BBMap parameters can be used; see bbmap.sh for more details. *****
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
	JNI=""
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=3200m" "--xms=3200m" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS $JNI -cp $CP align2.BBWrap build=1 overwrite=true fastareadlen=500 $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
