#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified July 26, 2026

Description:  Loads sequence or clade files and writes clade (spectra) files.
              Input may be fasta (tid in headers) OR existing clade/spectra
              files, which are re-read and merged -- so this tool also rewrites
              or augments a clade/spectra database (e.g. adddomain=t adds a field).

Usage: cladeloader.sh in=contigs.fa out=clades.clade

Parameters:
in=<file,file>  Input files: fasta with tid in headers, or clade/spectra files
                (re-read and merged).  Mixed input is allowed.
out=<file>      Output clade/spectra file.
adddomain=f     Embed a domain category (0-6: bacteria,archaea,virus,animal,
                plant,fungi,other-euk) in each output record, computed from the
                tree at build time, so downstream tools derive the domain without
                loading the tree.  Requires the tree (usetree=t, the default).
relineage=f     Regenerate each record's lineage string from the tree (up to the
                domain rank), replacing any stored lineage.  Enriches an existing
                clade/spectra DB's lineages with d__/sk__.  Requires the tree.
maxk=5          Limit max kmer length (range 3-5).
a48             Output counts in ASCII-48 instead of decimal (smaller).
a48o            Output counts in offset/delta ASCII-48 (smallest).
deco            Output counts in offset/delta decimal.
concise=t       Omit higher-k count arrays for records too short to use them.
16s=<file,file> Optional tax-labeled file of 16S sequences to attach.
18s=<file,file> Optional tax-labeled file of 18S sequences to attach.
replaceribo=f   Replace existing ssu with the ones supplied above.
callssu=t       Gene-call sequences over 900bp to find 16S/18S.  Slow on large
                inputs; set f when ssu is supplied via 16s= and 18s= instead.
usetree=t       Use a taxonomic tree, for names, levels, and lineage strings.
tree=auto       Tree location; auto finds it in resources, or give a path.
ddls=t          Build a DDL sketch per clade.  Set f when sketches are kept in
                a separate file.
ddl=<file>      Load DDL sketches from this file and attach them to clades by
                taxID (loaded in parallel with the main input).  Implies ddls=t.
ddlk=25         Kmer length for DDL sketches.
ddlbuckets=4096 Bucket count for DDL sketches.
ordered=f       Preserve input order in the output.
percontig=f     Emit one clade per sequence rather than per taxID.
keeptid=f       With percontig: label each record with the taxID parsed from its
                header (tid|NNN| or tid_NNN) instead of the sequence index, while
                still emitting one record per sequence -- sequences sharing a
                taxID stay separate records.  Headers with no parsable taxID
                keep the sequence-index label.  16s=/18s= attachment is by
                taxID and thus ambiguous here; attach ribo in a later pass.
mergedupes=f    Merge records sharing a taxID rather than keeping the largest.
addnovelonly=f  Add a record only if its taxID is not already present; on a
                duplicate, keep the EXISTING record (never merge or replace).
                For unioning a trusted base DB with best-effort additions, list
                the base file(s) first (files load in order).  Alias: keepexisting.
whitelist=<file>  Keep only records whose taxID appears in this file (one
                integer per line; '#' comments and blanks ignored).  Everything
                else is dropped at load, so it never enters the output.
                Designed to consume ddlmerger.sh tidsout= directly, making a
                spectra database congruent with a filtered sketch database.
                A whitelist subsumes a blacklist: keeping only the taxa the
                sketch database kept removes both the size-filtered taxa and any
                taxa this database holds alone, in one pass.
                Reports kept/dropped counts, and how many listed taxIDs had no
                record here (not an error, but a large number means the two
                databases disagree more than expected).
aligner=quantum Options include ssa2, glocal, drifting, banded, crosscut.

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

	parseJavaArgs "--xmx=4g" "--xms=4g" "--percent=42" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP clade.CladeLoader $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
