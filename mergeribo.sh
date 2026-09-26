#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified September 26, 2026

Description:  Merges ribosomal sequence files to keep ranked members per TaxID.
By default, a consensus is generated per TaxID, then the sequence
best matching that consensus is used:
First, all sequences per TaxID are aligned to a reference consensus.
Second, the best-matching sequence is used as a seed, and all
sequences for that TaxID are aligned to the seed to generate a new consensus.
Third, in 'consensus' mode, that consensus is simply output.
In 'best' mode (default), all sequences are aligned again to the new consensus,
and up to maxpertaxid members are output (one by default).
All copies, including duplicates, contribute to the per-TaxID consensus.
Optional exact deduplication occurs after ranking and before selecting outputs.

Usage:  mergeribo.sh in=<file,file> out=<file>
        mergeribo.sh in=5s.fa out=best5s.fa 5S=t maxpertaxid=10 dedupe=t

Standard parameters:
in=<file,file>  Comma-delimited list of files.
out=<file>      Output file.
overwrite=t     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.
t=auto          Number of worker threads.
reads=-1        Maximum reads from each input file; -1 means unlimited.
verbose=f       Print additional loading, alignment and selection diagnostics.
ordered=f       Use the ordered output writer. Taxon order still depends on
                worker completion, rather than input or taxonomic order.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression
                level; lower compression is faster.
fastawrap=70    4000 is recommended to minimize filesize.

Processing parameters:
alt=<file>      Lower priority data, used only for TaxIDs with no accepted
                representative in the primary input.
best=t          Output ranked real representatives per TaxID.
                Rank is identity times min(length,ideal)/max(length,ideal),
                with longer sequences first on ties. The FIRST reference sets
                the ideal length, even when another reference matches better.
maxpertaxid=1   Maximum representatives per TaxID; must be at least 1.
dedupe=f        Remove exact duplicate sequences within each TaxID AFTER
                consensus building and ranking, keeping the highest-ranked copy.
                Use maxpertaxid=10 dedupe=t for up to 10 distinct copies.
                With dedupe=f, repeated sequences may occupy multiple outputs.
consensus=f     Output a consensus per taxID instead of the best input
                sequence.  Mutually exclusive with best; requires maxpertaxid=1.
fast=f          Rank against the reference set without building a per-TaxID
                consensus; still honors maxpertaxid and dedupe.
                Fast mode and groups with fewer than 3 members use global-seed
                ranking; larger non-fast groups use their own consensus.
minid=0.62      Ignore sequences with identity lower than this to the global
                consensus.
maxns=-1        Ignore sequences with more than this many Ns, if non-negative.
minlen=1        Ignore sequences shorter than this; mode defaults listed below.
maxlen=4000     Ignore sequences longer than this; mode defaults listed below.
                Explicit minlen/maxlen override mode defaults in any arg order.
16S=t           Use the default 16S reference, or ref= if supplied.
18S=f           Use the default 18S reference, or ref= if supplied.
ITS=f           Use the default ITS lineage references, or ref= if supplied.
LSU=f           Large-subunit mode (23S/25S/26S/28S name one molecule); aligns
                to the prokaryotic 23S consensus, or to ref= if given (required
                for eukaryotic LSU, which is too divergent from 23S).
                Sets maxlen=6000 unless overridden. Mutually exclusive.
5S=f            Named 5S mode; use ALL shipped 5S reference records, or ALL
                records in ref=. Defaults to minlen=80 maxlen=200.
5.8S=f          5.8S mode; requires ref= (no configured MergeRibo default).
                Sets minlen=50 maxlen=300 unless overridden. Mutually exclusive.
                Modes are mutually exclusive; the last enabled mode wins.
                Aliases: process16S, process18S, processITS, processLSU,
                process5S, process5.8S; 58S also selects 5.8S.
ref=<file>      Override the reference set in EVERY mode, including 16S/18S/ITS
                (older versions ignored ref= in those modes). ALL records are
                used; each input receives its maximum identity across them.
                Older 5.8S used only the first record; now it also uses all.
                The first record still sets the ideal length for ranking.
                Empty or unreadable reference sets are errors.
level=          (taxlevel) If set to a term like 'species' or 'genus', nodes
                will be promoted to that level, minimum, before consensus.
dada2=f         Output headers in dada2 format.
tree=f          Load a TaxTree (default f); tree=t, tree=auto, or tree=<path> enables it.
usetree=f       Alias for tree=f.  dada2=t or level=<level> requires a tree; explicit f fails.
treefile=auto   Select the TaxTree source without changing the local gate.

TaxIDs: Headers may contain tid|N, tid_N, ncbi|N or ncbi_N anywhere in the
header. Invalid or missing positive TaxIDs are skipped. Inputs are unpaired.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
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
		LINK_TARGET="$(readlink "$SCRIPT")"
		SCRIPT="$LINK_TARGET"
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

	parseJavaArgs "--xmx=4g" "--xms=4g" "--percent=42" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prok.MergeRibo $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" prok.MergeRibo "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
