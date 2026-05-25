#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Noire
Last modified May 24, 2026

Description:  Identifies and classifies ribosomal SSU (16S/18S) and ITS
sequences using DynamicDemiLog (DDL) sketching against pre-built reference
databases.  Query type is determined automatically:
  Sequences aligning >64% to 16S or 18S consensus are classified as SSU.
  Sequences aligning <56% to all SSU consensuses are classified as ITS.
  Others are classified as unknown and compared to all reference types.

Usage:  findssu.sh ssu1.fa [ssu2.fa ...] [records=5]
    or: findssu.sh genome.fa call
    or: findssu.sh literal=ACGTACGT...
    or: findssu.sh name=Escherichia_coli
    or: findssu.sh name=Saccharomyces_cerevisiae its
    or: findssu.sh tid=562
    or: findssu.sh ssu.fa ref16s=<16S.tsv> ref18s=<18S.tsv>

Required resource files are loaded automatically from BBTools/resources/:
  ssuSketchDDL.tsv.gz          SSU DDL reference sketches (276k organisms)
  itsSketchDDL.tsv.gz          ITS DDL reference sketches (35k organisms)
  all_prok_16S_best_taxsorted.fa.gz   16S rRNA sequences (for alignment ANI)
  all_euk_18S_best_taxsorted.fa.gz    18S rRNA sequences (for alignment ANI)
  all_ITS_best_taxsorted.fa.gz        ITS sequences (for alignment ANI)
  16S_consensus_sequence.fa    16S consensus (for type classification)
  18S_consensus_sequence.fa    18S consensus (for type classification)
  ITS_*_consensus_sequence.fq  ITS consensuses (fungi, plant, animal, other)
If missing, download from:
  https://sourceforge.net/projects/bbmap/files/Resources/

SSU Mode (default):
Each input sequence is classified by alignment to consensus sequences
(16S, 18S, or ITS), then sketched with DDL and compared to the reference.

Call Mode:
Gene-calling identifies all SSU sequences in the input genome(s).
Each SSU found is individually sketched and classified.

Parameters:
ref=<file>      Pre-built SSU DDL reference file (TSV format).
                Default: resources/ssuSketchDDL.tsv.gz
ref16s=<file>   Separate 16S reference file.
ref18s=<file>   Separate 18S reference file.
refits=<file>   Separate ITS reference file.
                Default: resources/itsSketchDDL.tsv.gz (if present)
qf=<file>       Pre-built DDL query file for batch comparison.
call            Enable gene-calling mode for genomic input.
literal=<seq>   Provide a query sequence directly on the command line
                instead of from a file.
name=<name>     Look up a reference by organism name.  Accepts full names
                (name=Escherichia_coli), abbreviated (name=E.coli), or
                partial prefix matches.  Outputs TID, Type, Name, Sequence.
                Works via server (default) or locally.
tid=<int>       Look up a reference by NCBI TaxID (e.g. tid=562).
                Outputs TID, Type, Name, and Sequence.
                Works via server (default) or locally.
its             In lookup mode, return only ITS records.
16s             In lookup mode, return only 16S records.
18s             In lookup mode, return only 18S records.
ssu             In lookup mode, return only SSU records (16S + 18S).
                Flags are combinable: 'its 16s' returns both ITS and 16S.
records=5       Max hits to display per query.
minhits=8       Minimum shared index keys to compare a ref.
buffer=0        Alignment buffer size.  After index filtering, the top
                max(buffer, 20+2*records) candidates are aligned, then
                re-sorted by alignment ANI.  Bounds alignment cost while
                ensuring the best match is captured.
index=t         Use inverted index for query acceleration.
align=t         Perform SSU alignment for ANI calculation.
banself=f       Skip self-comparisons (when query and ref share a TaxID).
sequence=f      Print SSU/ITS sequence as last output column.
rank=f          Print rank column in output.
lineage=f       Print lineage column in output.
printname=t     Show Name column in output.
printtid=t      Show TID column in output.
loud=f          Print detailed timing and configuration info.
local=f         Force local processing (skip server, load refs locally).
server=t        Use the JGI SSU server (default).  Equivalent to local=f.
k=19            K-mer length for hashing.
buckets=128     Number of DDL buckets.
exponent=4      Exponent bits.
t=auto          Number of threads (default: all available cores).

Output columns:
ANI, WKID, Rank, Matches, Type, qLen, rLen, TID, Query, Name,
File, Contig, Start, Strand, Lineage

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

	parseJavaArgs "--xmx=3200m" "--xms=3200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ddl.SSUCompare $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
