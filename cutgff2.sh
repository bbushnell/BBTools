#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and G11
Last modified September 9, 2026

Description:  Streaming version of cutgff.sh: holds the GFF in memory and
STREAMS the fasta, so memory is proportional to the annotation plus one
batch of sequence rather than the whole genome.  Use this for huge inputs
(whole-clade fastas) where cutgff.sh's load-everything approach needs
hundreds of gigabytes.  Features are output in their sense strand.
Per-feature output is identical to cutgff.sh; output ORDER follows the
fasta (contig-major), not the gff.

Usage:  cutgff2.sh in=<fna file> gff=<gff file> out=<fna file>

Exactly one fna/gff pair per invocation.

File Parameters:
in=<file>           Input FNA (fasta) file.
gff=<file>          Input GFF file (optional; assumed from fasta name).
out=<file>          Output FNA file.

Other Parameters:
types=CDS           Types of features to cut.
invert=false        Invert selection: rather than outputting the features,
                    mask them with Ns and output the original sequences.
attributes=         A comma-delimited list of strings.  If present, one of
                    these strings must be in the gff line attributes.
bannedattributes=   A comma-delimited list of banned strings.
banpartial=t        Ignore lines with 'partial=true' in attributes.
minlen=1            Ignore lines shorter than this.
maxlen=2147483647   Ignore lines longer than this.
flank=0             Add this many bases of genomic flank to each side of
                    every extracted feature; same semantics as cutgff.sh.
gccontig=f          Append contig_gc= (GC of the full source contig).
filename=f          Append source=<fasta basename> to headers.
allowmissingseqids=f   By default, gff features whose seqid never appears
                    in the fasta are FATAL at end of input (the counts are
                    printed).  Set to t when the annotation legitimately
                    covers records absent from the fasta; the counts are
                    reported to stderr instead (record them).
workers=1           Worker threads.  1 (default) keeps output deterministic
                    and resources minimal; the pipeline is I/O-bound.

Notes:
  - Sequence headers like tid|taxid|ACC automatically fall back to bare ACC
    for gff-seqid matching.
  - If the same contig id appears twice, the FIRST occurrence gets the
    features; later ones are counted and reported (cutgff.sh was last-wins).
  - NOT supported here (use cutgff.sh): pickbest, oneperfile, multiple
    fna/gff pairs, renamebytaxid with taxmode accession/gi/header.

Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM.  The default of 1g is
                    enough for most annotations; raise it if the gff is huge.
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
	echo "java $EA $EOOM $SIMD $XMX $XMS -cp $CP gff.CutGff2 $@" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" gff.CutGff2 "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
