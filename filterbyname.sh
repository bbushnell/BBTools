#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified September 1, 2016

Description:  Filters reads by name.

Usage:  filterbyname.sh in=<file> in2=<file2> out=<outfile> out2=<outfile2> names=<string,string,string> include=<t/f>

in2 and out2 are for paired reads and are optional.
If input is paired and there is only one output file, it will be written interleaved.
Important!  Leading > and @ symbols are NOT part of sequence names;  they are part of
the fasta, fastq, and sam specifications.  Therefore, this is correct:
names=e.coli_K12
And these are incorrect:
names=>e.coli_K12
names=@e.coli_K12

Parameters:
include=f       Set to 'true' to include the filtered names rather than excluding them.
substring=f     Allow one name to be a substring of the other, rather than a full match.
                   f: No substring matching.
                   t: Bidirectional substring matching.
                   header: Allow input read headers to be substrings of names in list.
                   name: Allow names in list to be substrings of input read headers.
prefix=f        Allow names to match read header prefixes.
case=t          (casesensitive) Match case also.
ow=t            (overwrite) Overwrites files that already exist.
app=f           (append) Append to files that already exist.
zl=4            (ziplevel) Set compression level, 1 (low) to 9 (max).
int=f           (interleaved) Determines whether INPUT file is considered interleaved.
names=          A list of strings or files.  The files can have one name per line, or
                be a standard read file (fasta, fastq, or sam).
minlen=0        Do not output reads shorter than this.
ths=f           (truncateheadersymbol) Ignore a leading @ or > symbol in the names file.
tws=f           (truncatewhitespace) Ignore leading or trailing whitespace in the names file.
truncate=f      Set both ths and tws at the same time.

Positional parameters:
These optionally allow you to output only a portion of a sequence.  Zero-based, inclusive.
Intended for a single sequence and include=t mode.
from=-1         Only print bases starting at this position.
to=-1           Only print bases up to this position.
range=          Set from and to with a single flag.


Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

To read from stdin, set 'in=stdin'.  The format should be specified with an extension, like 'in=stdin.fq.gz'
To write to stdout, set 'out=stdout'.  The format should be specified with an extension, like 'out=stdout.fasta'

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

z="-Xmx800m"
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
	freeRam 800m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

function filterbyname() {
	local CMD="java $EA $EOOM $z -cp $CP driver.FilterReadsByName $@"
	echo $CMD >&2
	eval $CMD
}

filterbyname "$@"
