#!/bin/bash

usage(){
echo "
Written by UMP45
Last modified June 22, 2026
Description:  Shreds long single-molecule reads (e.g. PacBio HiFi) into fake
              paired-end reads, walking each molecule front-to-back and using
              every base exactly once.  Each 2*readlen window becomes one pair:
              first readlen bases = read 1, next readlen bases = read 2
              (reverse-complemented by default so pairs are FR-oriented).
              Rationale: a 2*readlen span can be placed with two readlen-bp
              alignments (cheaper than one long alignment, helps in repeats),
              and gives CallVariants a bbmap-style mapping.

Usage: shredpacbio.sh in=<file> out=<pairs> outs=<singletons> outd=<discards>

File Parameters:
in=<file>       Input long reads (fastq), single-ended.
out=<file>      Interleaved proper pairs (fastq).
outs=<file>     Singletons: a leftover read with no mate >= minlen.  Optional.
outd=<file>     Discards: leftover fragments shorter than minlen.  Optional.

Processing Parameters:
readlen=500     Length of each emitted read; a pair consumes 2*readlen bases.
                Aliases: rl, length, len.
minlen=200      Minimum length for any emitted read; shorter residuals go to
                discards.  Must be <= readlen.  Aliases: ml.
rcompmate=t     Reverse-complement read 2 (and reverse its quality) so each
                pair is FR-oriented like real Illumina data.  Aliases: rcomp, rc.
reads=-1        If nonnegative, stop after this many input molecules.
overwrite=t     (ow) Permit overwriting existing output files.

Residual handling (per molecule, after full 2*readlen pairs are taken):
  residual < minlen          -> discards (outd)
  [minlen, readlen)          -> one singleton (no room for a mate)
  [readlen, 2*readlen)       -> r1 is a full readlen read; the rest is its mate,
                                emitted as a (length-uneven) pair if the mate is
                                >= minlen, else r1 -> singleton and the sub-minlen
                                tail -> discards.

Output read headers (tab-delimited provenance):
  Pairs:      shredpb_<n> <1|2>:N:0 <tab> originalHeader <tab> pairStartInMolecule
              The '<1|2>:' token lets BBTools auto-detect interleaving; <n> is a
              global pair index shared by both reads of a pair.
  Singletons: shredpb_s<n> <tab> originalHeader <tab> startInMolecule
  Discards:   shredpb_d<n> <tab> originalHeader <tab> startInMolecule

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

if [ -z "$1" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
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

	parseJavaArgs "--xmx=4000m" "--xms=4000m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP synth.ShredPacBio $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
