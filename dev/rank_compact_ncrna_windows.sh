#!/bin/bash
# Compact consensus/member construction on already accepted saved windows.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/rank_compact_ncrna_windows.sh panels=panels.tsv selection=selection.tsv \
       baseline=corrected_loci.tsv.gz identities=accepted_rank_prefix \
       members=members.fa ledger=ledger.tsv assignments=assignments.tsv \
       census=models.tsv summary=members.summary.tsv ref=consensus.fa \
       out=fresh_prefix expectedmodels=508 expectedgroups=3153 t=1 -Xmx2g -eoom

Tests the simplest unpartitioned forward masked10 profile at512/1024 bins, with
one2048-bin reference. Each distinct raw consensus bucket/key gets initial weight
W once, plus actual member support; highest total wins, tie to highest key.
Weights1000000/1000/100/10 and controls1/0 run at512/1024. The2048 reference
uses1000000/0. Original member counts exclude pseudo-votes. Members hash once per
width. Weight one million is compared bin-for-bin with the completed member-count
plus-one rule without re-ranking that old profile. Differences include bucket
and old/new key; collision bins count distinct masked words BEFORE key compression.
Controls: highest-key member-only, consensus-only DDL, native-first hybrid at the
same width, original native and the accepted historical2048 lowest-key hybrid.

partitions=t runs only the separate compact partition designs: three-quarter
masked10 members plus one-quarter consensus, or one-quarter plain7/8 members,
one-half masked10 members and one shared consensus quarter. Both tagged lengths
contribute to that consensus sketch; queries probe their own member area plus
the shared quarter. Each complete bucket/key pair scores once. Widths512/1024
are decision sizes;2048 is one reference per design. Reuse accepted same-size S1
member/consensus/hybrid control tables; this mode does not rebuild them. Reserved
consensus-quarter occupancy is labelled separately from overall matching bins.

fineweights=t runs only NEW masked10 weights2/3/5/8/20/50 at1024bins.
Combine with the accepted weight10/100 rows externally; those finished comparisons
are not rerun. No native model-count cache or hybrid order is needed in this mode.
It cannot be combined with partitions=t. All other input bindings stay unchanged.

Reads accepted identity/bins/groups and original saved cache/reference/member
bindings; cannot compute new identities or read genomes. No caller or alignment.
Ranking includes zero scores in original-model tie order. Headline ranks use
the same original native-common loci and frozen low-identity bins; unknown and
unrankable targets and zero-target inputs remain in the outputs.

Writes .methods.tsv (stored-key bytes and consensus occupancy/length mean/range),
.profile_models.tsv (each model's counts/occupancy/length, raw consensus words,
collision bins and million-versus-always selected-key differences), .loci.tsv.gz,
.summary.tsv (ranks,within10/50,last-input losses versus native/historical hybrid)
and .last_losses.tsv plus .inputs.tsv. UInt16 DDL key payload is2*total_bins bytes
per model; total multiplies actual model count. Native keys are implicit dense
array addresses and add zero explicit key bytes; postings,headers,text size and
heap are excluded from these KEY-ONLY figures. They are not total index memory.

Single owner, reused native byte I/O/primitive scratch. Fresh outputs only;
ow=f,append=f,no sampling. Recipe must independently bind all accepted input
files. This is saved ranking evidence, not caller detection or adoption. Other
consensus tie-only/hash policies are not run here. Losses at50 are reporting rows,
not a rejection gate. Compact caller cap/floor/patience are tuned separately.
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
CP="${DIR}current/"
if [ -f "${DIR}bbtools.jar" ]; then CP="${DIR}bbtools.jar"; fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=2g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$CP" prok.RankCompactNcrnaWindows "$@"
