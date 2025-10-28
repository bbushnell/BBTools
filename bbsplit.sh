#!/bin/bash

usage(){
echo "
BBSplit
Written by Brian Bushnell, from Dec. 2010 - present
Last modified June 11, 2018

Description:  Maps reads to multiple references simultaneously.
Outputs reads to a file for the reference they best match, with multiple options for dealing with ambiguous mappings.

Usage: bbsplit.sh ref=x.fa,y.fa in=reads.fq basename=o%.fq
that is equivalent to
bbsplit.sh build=1 in=reads.fq ref_x=x.fa ref_y=y.fa out_x=ox.fq out_y=oy.fq

To index:     bbsplit.sh build=<1> ref_x=<reference fasta> ref_y=<another reference fasta>
To map:       bbsplit.sh build=<1> in=<reads> out_x=<output file> out_y=<another output file>

By default paired reads will yield interleaved output, but you can use the # symbol to produce twin output files.
For example, basename=o%_#.fq will produce ox_1.fq, ox_2.fq, oy_1.fq, and oy_2.fq.


Indexing Parameters (required when building the index):
ref=<file,file>     A list of references, or directories containing fasta files.
ref_<name>=<ref.fa> Alternate, longer way to specify references. e.g., ref_ecoli=ecoli.fa
                    These can also be comma-delimited lists of files; e.g., ref_a=a1.fa,a2.fa,a3.fa
build=<1>           Designate index to use.  Corresponds to the number specified when building the index.
path=<.>            Specify the location to write the index, if you don't want it in the current working directory.

Input Parameters:
in=<reads.fq>       Primary reads input; required parameter.
in2=<reads2.fq>     For paired reads in two files.
qin=<auto>          Set to 33 or 64 to specify input quality value ASCII offset.
interleaved=<auto>  True forces paired/interleaved input; false forces single-ended mapping.
                    If not specified, interleaved status will be autodetected from read names.

Mapping Parameters:
maxindel=<20>       Don't look for indels longer than this.  Lower is faster.  Set to >=100k for RNA-seq.
minratio=<0.56>     Fraction of max alignment score required to keep a site.  Higher is faster.
minhits=<1>         Minimum number of seed hits required for candidate sites.  Higher is faster.
ambiguous=<best>    Set behavior on ambiguously-mapped reads (with multiple top-scoring mapping locations).
                       best   (use the first best site)
                       toss   (consider unmapped)
                       random   (select one top-scoring site randomly)
                       all   (retain all top-scoring sites.  Does not work yet with SAM output)
ambiguous2=<best>   Set behavior only for reads that map ambiguously to multiple different references.
                    Normal 'ambiguous=' controls behavior on all ambiguous reads;
                    Ambiguous2 excludes reads that map ambiguously within a single reference.
                       best   (use the first best site)
                       toss   (consider unmapped)
                       all   (write a copy to the output for each reference to which it maps)
                       split   (write a copy to the AMBIGUOUS_ output for each reference to which it maps)
qtrim=<true>        Quality-trim ends to Q5 before mapping.  Options are 'l' (left), 'r' (right), and 'lr' (both).
untrim=<true>       Undo trimming after mapping.  Untrimmed bases will be soft-clipped in cigar strings.

Output Parameters:
out_<name>=<file>   Output reads that map to the reference <name> to <file>.
basename=prefix%suffix     Equivalent to multiple out_%=prefix%suffix expressions, in which each % is replaced by the name of a reference file.
bs=<file>           Write a shell script to 'file' that will turn the sam output into a sorted, indexed bam file.
scafstats=<file>    Write statistics on how many reads mapped to which scaffold to this file.
refstats=<file>     Write statistics on how many reads were assigned to which reference to this file.
                    Unmapped reads whose mate mapped to a reference are considered assigned and will be counted.
nzo=t               Only print lines with nonzero coverage.

***** Notes *****
Almost all BBMap parameters can be used; run bbmap.sh for more details.
Exceptions include the 'nodisk' flag, which BBSplit does not support.
BBSplit is recommended for fastq and fasta output, not for sam/bam output.
When the reference sequences are shorter than read length, use Seal instead of BBSplit.

Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.

This list is not complete.  For more information, please consult $DIRdocs/readme.txt
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
JNI="-Djava.library.path=""$DIR""jni/"
JNI=""

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

function bbsplit() {
	local CMD="java $EA $SIMD $EOOM $z $z2 $JNI -cp $CP align2.BBSplitter ow=t fastareadlen=500 minhits=1 minratio=0.56 maxindel=20 qtrim=rl untrim=t trimq=6 $@"
	echo $CMD >&2
	eval $CMD
}

bbsplit "$@"
