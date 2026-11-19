#!/bin/bash
# Summarize saved locus ranks and preserve named unresolved/zero-support cases.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_saved_rank_summary.sh in=saved_rank_prefix out=fresh_prefix methods=0,54,57,84 expectedmethods=148 expectedmodels=508 expectedloci=2156 expectedcommon=2007 expectedzero=3 expectedgaps=110 t=1 -Xmx512m -eoom
Reads .loci.tsv.gz and .parents.tsv.gz only; no sequences, caller or alignment.
Reports sums, means, within-10/50/100 and quantiles on the common native-rankable
loci. These are minimum ranks over saved events, not actual caller work or recall.
Optional strata=labels.tsv supplies every target in saved order with four fields:
input_slot, input_id, parent, seen_stratum (SEEN, NOT_SEEN or UNKNOWN).
Without that file all seen status is UNKNOWN, never assumed independent holdout.
Writes .summary.tsv, .zero_support.tsv and .no_compatible.tsv. The last list names
cases with no accepted compatible model under the captured windows/current verifier;
the cause remains unresolved between window geometry, verifier and library coverage.
No hybrid reconstruction is performed. Complete expected counts and fresh outputs
are required; sampling, append and overwrite are rejected.

To combine disjoint assembly groups, use panels=manifest.tsv instead of in=.
Supply expectedgroups= and the same methods/expectedmethods/expectedmodels.
Manifest fields (tab-separated): panel, rank_prefix, strata, assemblies, include,
expected_inputs, expected_selected_inputs, expected_loci, expected_common,
expected_zero, expected_gaps. Paths must be absolute. assemblies= the complete
seven-column assembly-contribution census; include is the four-column reused-group
selection (resolved_assembly_group,captured_slot,input_id,fasta), or . for all rows.
Every selected input must be a resolved noncontributor, with groups/IDs/paths
disjoint across panels. Every original input, parent, method and label is checked,
including excluded inputs and selected groups without target loci. Combined
summaries report NOT_SEEN, SEEN and UNKNOWN separately, with no pooled rows.
Named cases gain a panel column. Additional .inputs.tsv retains all selected
inputs, including zero-target groups, and their complete target dispositions.

fullcapture=t instead keeps contributing assemblies too. Every panel must use
include=.; all original inputs are candidates. The first captured representative
of each resolved assembly group is retained, with aliases named in .aliases.tsv.
expectedgroups=0 derives the unique census; a positive value checks that census.
Additional .panels.tsv binds the effective selected counts for downstream readers;
.population.tsv conserves captured inputs = distinct groups + discarded aliases.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.NcrnaSavedRankSummary "$@"
