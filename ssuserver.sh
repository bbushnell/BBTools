#!/bin/bash

usage(){
echo "
Written by Chloe and Brian Bushnell
Last modified May 20, 2026

Description:  Starts a persistent HTTP server for SSU (16S/18S) ribosomal
sequence classification using DDL sketching against a pre-built reference
database of 276k organisms.

The server preloads all reference data at startup and serves queries over
HTTP with sub-second response times.  Use findssu.sh with address= to
send queries to this server, or POST FASTA sequences directly.

Usage:  ssuserver.sh [port=3070] [ref=<ssu_ddls.tsv.gz>]

Server Parameters:
port=3070       HTTP listen port.
kill=<code>     Kill code for graceful remote shutdown via /kill/<code>.
prefix=<addr>   Restrict access to addresses starting with this prefix.
domain=<url>    CORS allowed origin (default: * for any origin).
ref=<file>      SSU DDL reference file (default: resources/ssuSketchDDL.tsv.gz).
ref16s=<file>   Separate 16S reference file.
ref18s=<file>   Separate 18S reference file.
k=19            K-mer length for hashing.
buckets=128     Number of DDL buckets.
exponent=4      Exponent bits.
records=5       Max hits to return per query.
minhits=8       Minimum shared index keys to compare a ref.
buffer=0        Alignment buffer size.
maxsize=100000000  Max request body size in bytes (default 100MB).
t=4             Number of handler threads.
verbose         Enable request logging.

Request Format:
  POST raw FASTA to / for SSU classification.
  Prefix body with //Call to enable gene-calling mode.
  Prefix body with //JSON for JSON output.
  GET / returns usage information.
  POST //Status returns server health check.

Java Parameters:
-Xmx            Set Java memory.  Default is 8g for the server.
-eoom           Exit on out-of-memory exception.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
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

	parseJavaArgs "--xmx=8g" "--xms=8g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ddl.SSUServer $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
