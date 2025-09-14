#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 17, 2015

Description:  Generates fake read pairs from ends of contigs or single reads.
Specifically for simulating a fake LMP library from long reads or an assembly;
for synthetic read generation from a reference see randomreads.sh or randomreadsmg.sh.

Usage:        bbfakereads.sh in=<file> out=<outfile> out2=<outfile2>

Out2 is optional; if there is only one output file, it will be written interleaved.

Standard parameters:
ow=f                (overwrite) Overwrites files that already exist.
zl=4                (ziplevel) Set compression level, 1 (low) to 9 (max).
fastawrap=100       Length of lines in fasta output.
tuc=f               (touppercase) Change lowercase letters in reads to uppercase.
qin=auto            ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.
qout=auto           ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input).
qfin=<.qual file>   Read qualities from this qual file, for the reads coming from 'in=<fasta file>'
qfout=<.qual file>  Write qualities from this qual file, for the reads going to 'out=<fasta file>'
qfout2=<.qual file> Write qualities from this qual file, for the reads coming from 'out2=<fasta file>'
verifyinterleaved=f (vint) When true, checks a file to see if the names look paired.  Prints an error message if not.
tossbrokenreads=f   (tbr) Discard reads that have different numbers of bases and qualities.  By default this will be detected and cause a crash.

Faking parameters:
length=250          Generate reads of this length.
minlength=1         Don't generate reads shorter than this.
overlap=0           If you set overlap, then reads will by variable length, overlapping by 'overlap' in the middle.
identifier=null     (id) Output read names are prefixed with this.
addspace=t          Set to false to omit the  space before /1 and /2 of paired reads.

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

z="-Xmx600m"
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

function fakereads() {
	local CMD="java $EA $EOOM $z -cp $CP synth.FakeReads $@"
	echo $CMD >&2
	eval $CMD
}

fakereads "$@"
