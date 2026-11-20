#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 16, 2026

***DEPRECATED***

Description:  Generates adapter-contaminated reads without using insert size.
The optional grade mode compares surviving read lengths with header annotations.
The input is a set of reads, paired or unpaired.
The output is those same reads with adapter sequence replacing some of the bases in some reads.
By default, paired reads use the same replacement position, independently of insert size.
The generated data can be used to benchmark adapter-trimming software.
Use randomreads.sh when paired adapter placement must respect insert size,
allowing overlap-based adapter detection.

Usage:  addadapters.sh in=<file> in2=<file2> out=<outfile> out2=<outfile2> adapters=<file>

in2 and out2 are for paired reads and are optional.
If input is paired and there is only one output file, it will be written interleaved.

Parameters:
ow=t                (overwrite) Overwrites files that already exist.
int=f               (interleaved) Determines whether INPUT file is considered interleaved.
qin=auto            ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.
qout=auto           ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input).
mode=add            Add adapters to input files.  Default mode.
mode=grade          Compare trimmed-read lengths with their recorded expected lengths.
adapters=<file>     Fasta file of adapter sequences.
literal=<sequence>  Comma-delimited list of adapter sequences.
left=f              Put adapters on the left (5') end of the read.
right=t             Put adapters on the right (3') end of the read.  Default mode.
end=right           Select right/3prime or left/5prime.
adderrors=t         Add errors to adapters based on the quality scores.
addpaired=t         Use the same replacement position for both mates, ignoring insert size.
arc=f               Add reverse-complemented adapters as well as forward.
rate=0.5            Add adapters to this fraction of reads.
minlen=1            Minimum expected remainder for validity statistics, not an output filter.
parsecustom=f       Enable generic SYN_ input parsing; explicit t requires valid SYN_ names.
                    Grade mode decodes length annotations locally with the default f.

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

	parseJavaArgs "--xmx=200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.AddAdapters $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" jgi.AddAdapters "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
