#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and G11
Last modified September 9, 2026

Description:  Emits fungal ITS regions (ITS1, ITS2, and the full contiguous
ITS1+5.8S+ITS2 span) from a caller GFF containing 18S, 5.8S and LSU/28S rRNA
rows — e.g. callgenes.sh output with 18s=t plus the r58/lsu ncRNA families.
The primary output is the full ITS span bounded by 18S and LSU; a 5.8S call
is not required for it when the outer pairing is unambiguous.  Individual
ITS1/ITS2 also emit when a uniquely paired 5.8S pivot exists.  Ambiguous,
overlapping, empty or inverted pairings are rejected with diagnosed statuses
rather than guessed; the report file records every pairing decision.

Usage:  callits.sh gff=<caller gff> out=<its gff>

File Parameters:
gff=<file>          Input GFF with the rRNA flank rows.
out=<file>          Output GFF of accepted ITS regions (types ITS, ITS1, ITS2).
report=<file>       Optional TSV recording every pairing decision, including
                    rejections (AMBIGUOUS_PAIRING, OVERLAPPING_FLANK,
                    REJECT_EMPTY, MISSING_18S, MISSING_LSU).
in=<file>           Genome fasta; required for outfasta.
outfasta=<file>     Optional fasta of region sequences, reverse-complemented
                    on the minus strand (coordinates in the GFF always stay
                    ascending-genomic).  NOTE: this path loads the genome
                    fully; it is intended for single genomes.  For huge or
                    bulk inputs, run without outfasta and extract sequences
                    from the emitted GFF with cutgff2.sh (streaming).

Other Parameters:
full=t              Emit the full-span ITS records.
individual=t        Emit individual ITS1/ITS2 records.
pattern18s=18S,     Attribute prefix identifying 18S rows.
pattern58s=model:r58_   Attribute substring identifying 5.8S rows.
patternlsu=model:lsu_   Attribute substring identifying LSU rows.

Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
-eoom               This flag will cause the process to exit if an
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

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

	parseJavaArgs "--xmx=1g" "--xms=1g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	echo "java $EA $EOOM $SIMD $XMX $XMS -cp $CP prok.CallITS $@" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" prok.CallITS "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
