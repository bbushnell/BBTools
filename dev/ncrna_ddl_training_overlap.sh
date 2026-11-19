#!/bin/bash
# Count actual kept-member assembly contributions through accepted aliases.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_ddl_training_overlap.sh in=assignments.tsv genomes=corpus.genomes.tsv groups=groups.inputs.tsv.gz sample=rank.inputs.tsv census=membership.models.tsv out=fresh_prefix t=1 -Xmx512m -eoom
Reads five saved metadata tables only; no sequences, caller or alignment.
Counts actual kept votes, excluding orphans and discarded small clusters, and
joins donor sources to sample inputs through accepted assembly-group aliases.
Reports sample paths and distinct resolved assembly groups separately. Missing
bindings and unresolved assembly pairs remain explicit. It does not determine
exact-sequence overlap, same-source-locus overlap, or independent holdout status.
Writes .inputs.tsv and .summary.tsv. Expected complete counts default to
expectedinputs=1000 expectedgenomes=39979 expectedmodels=508
expectedmembers=85741 expectedorphans=15 expectedsmall=0.
Fresh outputs and complete populations are required; no overwrite or sampling.
sample= also accepts the actual eight-column SelectNcrnaBenchmarkInputs selection
table. To select additional noncontributing assembly groups, supply the COMPLETE
evaluation selection with expectedinputs=11769 and captured=prior.rank.inputs.tsv
expectedcaptured=1000. Writes .reuse.tsv, .uncaptured.selection.tsv,
.uncaptured.config and .expansion.tsv. Already captured aliases take precedence;
one representative per assembly, explicit unresolved groups excluded. Sequence
overlap is still unassessed. No genome is opened and no work is submitted.
HELP
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ]; do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
if [ -f "${DIR}bbtools.jar" ]; then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=512m --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.NcrnaDdlTrainingOverlap "$@"
