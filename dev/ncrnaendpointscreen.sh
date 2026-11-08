#!/bin/bash
# Screens native endpoint-count geometry; no network training or gene calling.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden

Usage: dev/ncrnaendpointscreen.sh in=contexts.fa.gz out=fresh_directory \
  graph=original.support.tsv original=original_consensus.fa \
  roster=roster.tsv loci=loci.tsv splitseed=1 ks=7,9 radius=10

Screens all inside/outside compositions for each K using the native whole-TaxID
split and endpoint tables. Inputs are the complete original admitted population,
not a small benchmark or reduced model membership. No exclusions or row sampling.

in=          Complete sense-oriented contexts with source, actual flanks and GC.
out=         Fresh output directory beneath an existing parent.
graph=       Complete ORIGINAL support graph.
original=    Consensus FASTA matching the graph's original model-name axis.
roster=      Native eleven-column original source roster.
loci=        Native 21-column original locus metadata.
splitseed=   Required nonnegative whole-TaxID split seed.
ks=7,9       Distinct K values, each 1..15. All outside counts 0..K are tested.
radius=10    Positive symmetric candidate radius in bases, with offset zero included.
ratios=1.5,2,3,5,10
             Increasing finite count-decisiveness ratios, all greater than one.
t=           Workers for native table construction. Scoring uses one worker.

Tables learn only eligible internal TRAIN loci. TRAIN_LOO removes that locus's
true word at EVERY matching occurrence through the existing native LOO view;
VALIDATION uses unchanged training counts. START and STOP use their distinct
native inside/outside formulas. Reports keep endpoints and partitions separate;
validation never contributes to the learned counts.

For each feature-complete eligible locus, score every candidate whose gene span
is at least 16 bases and fits the context, and whose entire tip word is defined
and available. Missing/ambiguous words and invalid spans have separate counters.
A curve requires a valid true word and at least one valid nonzero-offset word;
incomplete competitor sets remain eligible, with their missing counts visible.

Compare true+1 against max(other)+1. At ratio M, decisive_true means true+1 is
at least M times every valid other count+1; decisive_other means a valid other
count+1 is at least M times true+1. Otherwise abstain. Ties abstain for every M>1.
Fractions use eligible_loci, not all original loci; empty denominators are NA.
These are count-screening curves, NOT caller FPR/recall or neural accuracy.

Outputs: population.tsv (exact split/eligibility), summary.tsv (original/context,
table and candidate accounting), curves.tsv (counts and fractions by ratio).
Missing GC prevents scoring, as in native vectors, but a geometrically valid
training record can still contribute to the native endpoint table.
Only one geometry's tables are retained. Two passes over the existing small
context file per geometry; no original source-genome reread or HBM alignment.
Final K/placement choices still need native network and whole-genome evaluation.
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
DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)"
if [ -f "$DIR/bbtools.jar" ]; then
	CP="$DIR/bbtools.jar"
else
	CP="$DIR/current/"
fi
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=4g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$CP" prok.NcrnaEndpointScreen "$@"
