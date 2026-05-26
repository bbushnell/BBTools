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

=== Whole-genome blacklists (k=31, buckets=2048, exponent=5) ===

Step 1: Build double-sized sketches (buckets=4096) with kmers=t per clade:
  ddlwriter.sh in=refseq.CLADE.fna.gz out=CLADE.tsv.gz k=31 buckets=4096 \\
    exponent=5 mode=pertid kmers=t lineage=t tossjunk t=32

Step 2: Generate combined genus-promoted kmer dump across all clades:
  ddlblacklist.sh in=archaea.tsv.gz,bacteria.tsv.gz,...,viral.tsv.gz \\
    out=combined_dump.fa.gz k=31 exponent=5 mincount=5

  Files processed sequentially; each clade's sketch data is discarded
  after kmer extraction.  TaxIDs promoted to genus level during collection.
  304M distinct kmers across 200k records from 12 RefSeq clades.

Step 3: Filter at chosen thresholds using awk on combined_dump.fa.gz headers,
  merge with deduplication.  Cutoffs chosen from cumulative distribution plot
  (Y = kmers with count >= X):
  genus/250, family/50, order/20, class/10, phylum/5

  Individual level counts: g/400=3853, f/75=6309, o/30=6785, c/10=6199, p/5=8166
  Merged: 15,726 unique kmers

Step 4: Validate via condense mode (no re-sketching needed):
  ddlblacklist.sh in=archaea.tsv.gz,...,viral.tsv.gz condense=t validate=t \\
    blacklist=bl_merged.fa k=31 exponent=5 samples=200000

  Condenses 4096->2048 buckets while preferring non-blacklisted kmers.
  Approximates real blacklist effect without reading terabytes of genomes.

Step 5: Iterate (optional, marginal benefit for genomes):
  Re-sketch with blacklist active to expose second-tier masked kmers:
    ddlwriter.sh ... blacklist=bl_merged.fa
  Re-run ddlblacklist.sh on new sketches with same thresholds.
  Merge iter1+iter2 blacklists.  For RefSeq genomes, iteration added
  2,022 new kmers with negligible noise improvement
  (0.70 -> 0.71), confirming convergence after one iteration.

Final result: genomeDDLBlacklist.fa.gz (23,020 kmers, 571KB)
  Noise floor: 4.66 -> 0.71 avg shared keys at class+ distance (6.6x)
  Genus signal preserved: 145.4 -> 143.7 (99%)
  Autoloaded by DDLCompare.  Also loaded by DDLWriter via blacklist= flag.

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
