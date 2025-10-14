#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 12, 2025

Description:  Calculates some scalars from nucleotide sequence data.
Prints the averages for each input file.
Also prints standard deviation of each file if windowed.

Usage:  scalars.sh in=<input file> out=<output file>


Standard parameters:
in=<file>       Primary input; fasta or fastq.
                This can also be a directory or comma-delimited list.
		Filenames can also be used without in=
out=stdout      Set to a file to redirect output.

Processing parameters:
header=f        Print a header line.
rowheader=f     Print a row header.
window=0        If nonzero, calculate and average over windows.
break=f         Set to true to break data at contig bounds,
                in windowed mode.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
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

calcXmx () {
    # Source the new scripts
    source "$DIR""/memdetect.sh"
    source "$DIR""/javasetup.sh"

    parseJavaArgs "--mem=800m" "--mode=fixed" "$@"

    # Set environment paths
    setEnvironment
}
calcXmx "$@"

scalars() {
	if [[ $# -eq 0 ]]; then
		usage
		return
	fi
	local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP scalar.Scalars $@"
	#echo $CMD >&2
	eval $CMD
}

scalars "$@"
