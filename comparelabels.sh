#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 9, 2024

Description:  Compares delimited labels in read headers to count 
how many match.  The 'unknown' label is a special case.  The original 
goal was to measure the differences between demultiplexing methods.
Labels can be added with the rename.sh suffix flag, or the 
novademux.sh rename+nosplit flags, or seal.sh with rename, addcount=f,
and tophitonly.  The assumption is that a header will look like:
@VP2:12:H7:2:1101:8:2 1:N:0:CAAC (tab) CAAC (tab) CAAC
...in which case the labels CAAC would be compared and found equal.


Usage:  comparelabels.sh in=<input file> out=<output file>

Input may be fasta or fastq, compressed or uncompressed.

Standard parameters:
in=<file>       Primary input, or read 1 input.
out=stdout      Print the results to this destination.  Default is stdout
                but a file may be specified.
labelstats=     Optional destination for per-label stats.
quantset=<file> If set, ignore reads with labels not contained in this file;
                one label per line.  'unknown' is automatically included.
swap=f          Swap the order of label 1 and label 2.
delimiter=tab   Compare the last two terms in the header, using this 
                single-character delimiter.  Most symbols can be expressed
                as literals (e.g. 'delimiter=_' for underscore) but you can
                also spell out some of the problematic ones:
                   space, tab, pound, greaterthan, lessthan, equals,
                   colon, semicolon, bang, and, quote, singlequote,
                   backslash, hat, dollar, dot, pipe, questionmark, star,
                   plus, openparen, closeparen, opensquare, opencurly

Output Terminology:
aa             Both labels were equal.
uu             Both labels were unknown.
au             Label 1 was assigned, label 2 was unknown.
ua             Label 1 was unknown, label 2 was assigned.
ab             Both labels were assigned, but not equal.  For per-label
               stats, indicates label 1 was assigned to this, and label 2
               was assigned to something else.
ba             In per-label stats, indicates label 2 was assigned to this and
               label 1 was assigned to something else. 
yield          Fraction of reads assigned to the same label.  E.g. if aa=10,
               au=1, ab=2, then yield2 = aa/(aa+au+ab) = 10/13 = 0.77.
contam         Fraction of reads assigned to a different label, using the
               other as ground truth.  For example, if aa=10, au=1, ab=2,
               then contam1=ab/(aa+au+ab) = 2/13 = 153846 PPM.

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

z="-Xmx300m"
z2="-Xms300m"
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

function comparelabels() {
	local CMD="java $EA $SIMD $EOOM $z $z2 -cp $CP barcode.CompareLabels $@"
	echo $CMD >&2
	eval $CMD
}

comparelabels "$@"
