#!/bin/bash
# Exact source-qualified exclusion lists; no sequence analysis or neural fitting.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Export genome-level panel exclusions for ncRNA boundary testing.
Usage: dev/ncrnalearningexclusions.sh in=corpus_MANIFEST.tsv panel=panel1000.list \
 optionalpanel=old300.list original=original_consensus.fa graph=original.support.tsv \
 roster=genomes.tsv loci=loci.tsv out=fresh_directory t=1

panel= lists exact FASTA paths that MUST all map to the working corpus manifest.
optionalpanel= may have unmatched path aliases; every one is reported and must
be reconciled against authoritative assembly/source bindings before any fit.
An unmatched path does NOT prove that its genome is outside the working corpus.
Only the manifest's assembly and source_fasta paths are accepted aliases.
Basenames and TaxIDs never identify a held-out assembly. The existing population
loader validates source-qualified graph/roster/locus joins before export.

Outputs: exclude.txt (exact original IDs, no header), panel_sources.tsv and
summary.tsv (original, excluded, retained and eligible counts). This prepares
HELD OUT (panel genomes excluded from tables and training) bindings; it does not
fit anything or certify subsequent training. Compare remaining minor-member
populations before fitting. Relatives under different assembly IDs can remain.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.NcrnaLearningExclusions "$@"
