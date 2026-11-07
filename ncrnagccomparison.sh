#!/bin/bash
# Diagnostic measurement only; never trains or changes resources.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden

Usage: ncrnagccomparison.sh in=contexts.fa.gz out=fresh_directory family=tpp \
  graph=original.support.tsv original=original_consensus.fa roster=roster.tsv \
  loci=loci.tsv ref=accepted_consensus.fa hbm=accepted_models.hbmc \
  seeds=seeds.fa.gz splitseed=1

Measures GC differences through the actual boundary vector writer, classifier
reader, regression reader and Scavenger network input. No network is trained.
Uses the full existing original population and contexts; only selected source
genomes are reread. Inputs must be the already admitted training population.

in=          Complete sense-oriented contexts with contig/source/GC/flank fields.
out=         Fresh directory beneath an existing parent; never overwrites.
graph=       Complete ORIGINAL support graph, not reduced membership.
original=    Consensus FASTA matching that graph's model axis.
roster=      Native eleven-column source roster with literal FASTA paths.
loci=        Native 21-column corrected locus ledger with contiguous intervals.
ref=         Accepted consensus library used for actual feature generation.
hbm=         Index-aligned accepted HBM models.
seeds=       Accepted forward seed FASTA.
family=      Registered shared-boundary family name.
splitseed=   Required nonnegative whole-TaxID split seed. Exclusions do not change
             this split. Diagnostic refit includes selected training loci only.
requireexact=f  If true, fail on any bit difference between raw vector features,
             both native readers and the scorer using the same leave-one-out
             counts. Also require sample GC equality with actual Scavenger input.

One lexical source/contig/locus is selected per GC decile and contig-length bin
(<=1k, <=10k, <=100k, <=1M, <=16777216, >16777216), at most 60 real contigs.
Empty bins and missing-model loci remain explicit. A single worker and shard
preserve row order. Radius-one vectors contain each true end plus valid candidates
one base either side. Both labels are required by the native classifier reader.
Sampled endpoint counts are diagnostic inputs, not release resources.

Outputs: bins.tsv, exclusions.txt, vectors.train.0.tsv.gz, loci.tsv,
differences.tsv, accounting.tsv and features.tsv. Differences report exact float-bit counts,
absolute-error and float-step-distance nearest-rank median/95th/max values.
The GC report measures actual complete-contig inference. The ten-feature report
compares arithmetic with MATCHING leave-one-out counts; it does not equate training
counts with ordinary full-count inference or measure neural output sensitivity.
Zero observed mismatches establish sample agreement only. No tolerance is chosen.
HELP
}
if [ $# -eq 0 ] || [ "$1" = -h ] || [ "$1" = --help ]; then
	usage
	exit 0
fi
SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
while [ -h "$SCRIPT" ]; do
	DIR="$(dirname "$SCRIPT")"
	RESOLVED_SCRIPT="$(readlink "$SCRIPT")"
	SCRIPT="$RESOLVED_SCRIPT"
	if [ "${SCRIPT#/}" = "$SCRIPT" ]; then
		RESOLVED_SCRIPT="$DIR/$SCRIPT"
		SCRIPT="$RESOLVED_SCRIPT"
	fi
done
DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
if [ -f "$DIR/bbtools.jar" ]; then
	CP="$DIR/bbtools.jar"
else
	CP="$DIR/current/"
fi
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=4g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$CP" prok.NcrnaGcComparison "$@"
