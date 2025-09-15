#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified June 22, 2016

Description:  Summarizes the stats output of Seal for evaluation of 
cross-contamination.  The intended use is to map multiple libraries or 
assemblies, of different multiplexed organisms, to a concatenated reference 
containing one fused scaffold per organism.  This will convert all of the 
resulting stats files (one per library) to a single text file, with multiple 
columns, indicating how much of the input hit the primary versus nonprimary 
scaffolds.

If ingoresametaxa or ignoresamebarcode are used, ref names must be 
in this format:
barcode,library,tax,location
For example:
6-G,N0296,gammaproteobacteria_bacterium,deep_ocean

Usage:  summarizeseal.sh in=<file,file...> out=<file>

You can alternately run 'summarizeseal.sh *.txt out=out.txt'

Parameters:
in=<file>              A list of stats files, or a text file containing one stats file name per line.
out=<file>             Destination for summary.
printtotal=t           (pt) Print a line summarizing the total contamination rate of all assemblies.
ignoresametaxa=f       Ignore secondary hits sharing taxonomy. 
ignoresamebarcode=f    Ignore secondary hits sharing a barcode.
ignoresamelocation=f   Ignore secondary hits sharing a sampling site.
totaldenominator=f     (td) Use all bases as denominator rather than mapped bases.

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

z="-Xmx120m"
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

summarizeseal() {
	local CMD="java $EA $EOOM $z -cp $CP driver.SummarizeSealStats $@"
#	echo $CMD >&2
	eval $CMD
}

summarizeseal "$@"
