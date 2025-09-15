#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 17, 2017

Description:  Masks sequences of low-complexity, or containing repeat kmers, or covered by mapped reads.
By default this program will mask using entropy with a window=80 and entropy=0.75
Please read bbmap/docs/guides/BBMaskGuide.txt for more information.

Usage:   bbmask.sh in=<file> out=<file> sam=<file,file,...file>

Input may be stdin or a fasta or fastq file, raw or gzipped.
sam is optional, but may be a comma-delimited list of sam files to mask.
Sam files may also be used as arguments without sam=, so you can use *.sam for example.
If you pipe via stdin/stdout, please include the file type; e.g. for gzipped fasta input, set in=stdin.fa.gz

Input parameters:
in=<file>           Input sequences to mask. 'in=stdin.fa' will pipe from standard in.
sam=<file,file>     Comma-delimited list of sam files.  Optional.  Their mapped coordinates will be masked.
touppercase=f       (tuc) Change all letters to upper-case.
interleaved=auto    (int) If true, forces fastq input to be paired and interleaved.
qin=auto            ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.

Output parameters:
out=<file>          Write masked sequences here.  'out=stdout.fa' will pipe to standard out.
overwrite=t         (ow) Set to false to force the program to abort rather than overwrite an existing file.
ziplevel=2          (zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster.
fastawrap=70        Length of lines in fasta output.
qout=auto           ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input).

Processing parameters:
threads=auto        (t) Set number of threads to use; default is number of logical processors.
maskrepeats=f       (mr) Mask areas covered by exact repeat kmers.
kr=5                Kmer size to use for repeat detection (1-15).  Use minkr and maxkr to sweep a range of kmers.
minlen=40           Minimum length of repeat area to mask.
mincount=4          Minimum number of repeats to mask.
masklowentropy=t    (mle) Mask areas with low complexity by calculating entropy over a window for a fixed kmer size.
ke=5                Kmer size to use for entropy calculation (1-15).  Use minke and maxke to sweep a range.  Large ke uses more memory.
window=80           (w) Window size for entropy calculation.
entropy=0.70        (e) Mask windows with entropy under this value (0-1).  0.0001 will mask only homopolymers and 1 will mask everything.
lowercase=f         (lc) Convert masked bases to lower case.  Default is to convert them to N.
split=f             Split into unmasked pieces and discard masked pieces.

Coverage parameters (only relevant if sam files are specified):
mincov=-1           If nonnegative, mask bases with coverage outside this range.
maxcov=-1           If nonnegative, mask bases with coverage outside this range.
delcov=t            Include deletions when calculating coverage.
NOTE: If neither mincov nor maxcov are set, all covered bases will be masked.

Other parameters:
pigz=t              Use pigz to compress.  If argument is a number, that will set the number of pigz threads.
unpigz=t            Use pigz to decompress.

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
z2="-Xms1g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

bbmask() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP jgi.BBMask $@"
	echo $CMD >&2
	eval $CMD
}

bbmask "$@"
