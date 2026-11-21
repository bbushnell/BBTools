#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 19, 2026

Description:  In progress.
Generates a fixed nine-column, tab-separated assembly statistics table for multiple files.

Usage:        stats3.sh in=file
Or:           stats3.sh in=file,file
Or:           stats3.sh file file file

Parameters:
in=file         Specify the input fasta file(s), or stdin.
                Aliases: input, in1, input1.
                Multiple files can be listed without a 'in=' flag.
out=stdout      Destination of primary output; may be directed to a file.
extout=null     Output-format overrides are unsupported; nonempty values are rejected.
                The last extout setting applies; null or an empty value clears it.
reads=-1        Maximum nonempty sequence records per input, in input order.
                Alias: maxreads. Negative is unlimited; zero reports empty statistics.
                Empty header records do not count; preheader sequence counts as a record.
ow=t            Allow replacing existing output files when app=f.
app=f           Append a complete statistics table, including its header.
                When enabled, append preserves existing contents and takes precedence over ow.

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

	parseJavaArgs "--xmx=120m" "--xms=120m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.AssemblyStats3 $@"
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" jgi.AssemblyStats3 "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
