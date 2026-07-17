#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Noire
Last modified May 25, 2026

Description:  Builds a DDL kmer blacklist from pre-built DDL sketch files.
Reads sketches with kmer arrays (built with ddlwriter.sh kmers=t), promotes
each TaxID to genus level, counts distinct genera per kmer across all input
files, and outputs overrepresented kmers as FASTA with multi-level taxonomic
counts in the header.  Accepts comma-separated inputs; processes each file
sequentially and discards sketch data after kmer extraction.

Also supports condense mode (condense=t blacklist=X): condenses double-sized
sketches to half size while preferring non-blacklisted kmers, then validates.
Approximates full re-sketching without re-reading raw genome data.

Usage:  ddlblacklist.sh in=sketches.tsv out=dump.fa [mincount=5] [k=31]
        ddlblacklist.sh in=a.tsv,b.tsv out=dump.fa k=31 exponent=5 mincount=5
        ddlblacklist.sh in=sketches.tsv condense=t validate=t blacklist=bl.fa

Parameters:
in=<file>       Input DDL sketch file(s), comma-delimited.
                Must have been built with kmers=t.
out=<file>      Output FASTA of overrepresented kmers.
                Headers: >kmer_HEX raw=N g=N f=N o=N c=N p=N k=N sk=N e=0.XXX
mincount=5      Minimum genus count to include in output.
k=19            K-mer length (must match input sketches).
exponent=5      Exponent bits (must match input sketches).
blacklist=<file> FASTA blacklist for condense mode.
condense=f      Condense 4096-bucket sketches to 2048, preferring
                non-blacklisted kmers.  Use with validate=t.
validate=f      Validate shared keys at each taxonomic distance.
samples=200000  Number of random pairs for validation.

Output header format:
  raw = distinct genus-level TaxIDs containing this kmer
  g/f/o/c/p/k/sk = distinct taxa at genus/family/order/class/phylum/kingdom/
                    superkingdom level (promoted from genus TaxIDs)
  e = Shannon entropy of the kmer's trimer distribution

--- Blacklist Generation Methodology ---

=== SSU/ITS ribosomal blacklists (k=19, buckets=128, exponent=4) ===

Build double-sized sketches (buckets=256) with kmers=t from ribo databases:
  ddlwriter.sh in=16S.fa.gz out=16S_sketches.tsv k=19 buckets=256 exponent=4 \\
    mode=pertid kmers=t lineage=t

Generate blacklists at each level (run once per type per level):
  ddlblacklist.sh in=16S_sketches.tsv out=bl_family.fa k=19 exponent=4 \\
    mincount=160 validate=t samples=200000

Tuned cutoffs (each level run independently, merged with MergeDDLBlacklists):
  16S:  family/160, order/80,  class/32, phylum/18  -> 1,219 unique kmers
  18S:  family/190, order/70,  class/25, phylum/12  -> 2,758 unique kmers
  ITS:  family/200, order/70,  class/24, phylum/8   ->   217 unique kmers
  Combined: 4,230 unique kmers -> riboDDLBlacklist.fa.gz
  Noise reduction: 10-29x.  Speed improvement: 2.6x on all-to-all.

=== Whole-genome blacklists (k=25, exponent=5; production buckets=4096) ===

IMPORTANT - use SOFT (high) cutoffs for genome blacklists.  A genome blacklist
should remove only kingdom/phylum-level "universal" noise (the plant-vs-octopus
kmers) while PRESERVING the genus/family signal needed for distant-relative
detection.  Aggressive low cutoffs (e.g. genus/250) strip phylogenetic signal
and must NOT be used here.  Blacklisted kmers have higher trimer entropy (0.919)
than the bulk of genus>=5 kmers (0.890) - they are compositionally generic.

Step 1: Build sketches with kmers=t per clade at the target bucket count:
  ddlwriter.sh in=refseq.CLADE.fna.gz out=CLADE.tsv.gz k=25 buckets=4096 \\
    exponent=5 mode=pertid kmers=t lineage=t tossjunk t=32
  (Denser blacklists use buckets=32768 or 65536 - one dump per bucket count.)

Step 2: Generate combined genus-promoted kmer dump across all clades:
  ddlblacklist.sh in=archaea.tsv.gz,bacteria.tsv.gz,...,viral.tsv.gz \\
    out=combined_dump.fa.gz k=25 exponent=5 mincount=5
  TaxIDs promoted to genus level during collection; sketch data discarded
  after kmer extraction.

Step 3: Filter at SOFT thresholds (union) using awk on the dump headers.
  A kmer is blacklisted if it passes ANY of:
    genus g>=1600 | family f>=270 | order o>=110 | class c>=40
  NO phylum/kingdom filter (those levels are already covered by genus+family).

  A single genomic blacklist ships: refseqGenomeDDLBlacklist_k25e5b65536.fa.gz (5,944 kmers, built
  from the 65536-bucket dump).  It is bucket-count-independent (just a kmer set), so it masks 4k/32k/64k
  sketches alike and is the default for both DDLWriter (build time) and QuickClade (query time).
  Counts scale gently with the dump's bucket count (4096->1,651; 32768->3,607; 65536->5,944).
  Noise floor (class+ distance, 2M pairs): 5.06 -> 2.64; genus signal kept.

Step 4 (optional): kcompress fuse=160 shrinks the blacklist for distribution;
  verify lossless with kmercountexact (symmetric difference MUST be 0).  The
  *_fused.fa.gz forms require the loadBlacklist sliding-window reader.

Loaded by DDLWriter (blacklist= at sketch build time) and by CladeSearcher /
QuickClade (blacklist= applied to the query).  DDLCompare autoloads the default.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM.  The max is typically
                85% of physical memory.
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

	parseJavaArgs "--xmx=4g" "--xms=4g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ddl.DDLBlacklistMaker $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
