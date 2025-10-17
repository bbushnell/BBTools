#!/bin/bash

usage(){
echo "
Written by Chloe
Last modified October 15, 2025

Description:  Converts Illumina CBCL (Compressed Base Call) files to text format.
Extracts base calls, quality scores, and flowcell coordinates from binary CBCL files.

Usage:  cbcl2text.sh runfolder=<path> out=<file> lane=<int>


Standard parameters:
runfolder=<dir>  Path to Illumina run folder containing Data/Intensities.
out=<file>       Output file (tab-delimited text).
lane=<int>       Lane number to process (default 1).

Optional parameters:
tiles=<list>     Comma-separated tile numbers (e.g., tiles=1101,1102).
                 Default: process all tiles found in lane directory.
length=<mode>    Read splitting mode:
                   (none)          - Concatenate all cycles (default)
                   auto            - Parse RunInfo.xml for read structure
                   151,19,10,151   - Manual read lengths (comma-delimited)

Output format (default):
tile    X       Y       PF      bases(all_cycles)       quals(all_cycles)

Output format (with length):
tile    X       Y       PF      R1,I1,I2,R2             Q1,QI1,QI2,Q2

Coordinates:
X and Y are transformed to Illumina FASTQ format: round(10*raw + 1000)

Quality scores:
Illumina bins qualities to 2 bits (values 0-3) in CBCL files.

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

    parseJavaArgs "--mem=8g" "--mode=fixed" "$@"

    # Set environment paths
    setEnvironment
}
calcXmx "$@"

cbcl2text() {
	if [[ $# -eq 0 ]]; then
		usage
		return
	fi
	local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP illumina.Cbcl2Text $@"
	#echo $CMD >&2
	eval $CMD
}

cbcl2text "$@"
