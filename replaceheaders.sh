#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 23, 2016

Description:  Replaces read names with names from another file.
The other file can either be sequences or simply names, with
one name per line (and no > or @ symbols).  If you use one name
per line, please give the file a .header extension.

Usage:  replaceheaders.sh in=<file> hin=<headers file> out=<out file>

Parameters:
in=                 Input sequences.  Use in2 for a second paired file.
in=                 Header input sequences.  Use hin2 for a second paired file.
out=                Output sequences.  Use out2 for a second paired file.
ow=f                (overwrite) Overwrites files that already exist.
zl=4                (ziplevel) Set compression level, 1 (low) to 9 (max).
int=f               (interleaved) Determines whether INPUT file is considered interleaved.
fastawrap=70        Length of lines in fasta output.
qin=auto            ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.
qout=auto           ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input).

Renaming mode parameters (if not default):
addprefix=f         Rename the read by prepending the new name to the existing name.

Sampling parameters:
reads=-1            Set to a positive number to only process this many INPUT reads (or pairs), then quit.

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

#This block allows symlinked shellscripts to correctly set classpath.
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx1g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
}
calcXmx "$@"

function reheader() {
	local CMD="java $EA $SIMD $EOOM $z -cp $CP jgi.ReplaceHeaders $@"
	echo $CMD >&2
	eval $CMD
}

reheader "$@"
