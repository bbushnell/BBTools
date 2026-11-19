#!/bin/bash
# Raw source-locus projection and independent native CutGff output validation.
set -eo pipefail
if [[ $# == 0 || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/cobalaminconstructioninputs.sh mode=prepare in=selected_MANIFEST.tsv loci=native_loci.tsv.gz out=fresh_directory expectedgenomes=40686 expected=91072,2270,1786 chunkgenomes=1000 seed=20261115 t=1 -Xmx4g -eoom
mode=prepare writes per-assembly extraction GFFs, source-ID provenance, and
disjoint fixed-seed shuffled chunk configs. All selected genomes participate in
chunk assignment; a zero-locus assembly requires no extraction read. Native
truncation alone excludes construction. Equal DNA never merges source records.
Run native cutgff.sh with each chunk config and the explicit extraction recipe.
mode=validate uses the SAME manifest/loci/expected/chunkgenomes/seed arguments,
plus sequences=chunk.fa chunk=N and another fresh out= directory. It checks the
exact admitted source-ID set and lengths, rejects unexpected/duplicate records,
and splits accepted originals into the three member FASTAs. On failure only,
affected assemblies are read for contig-length diagnostics. complete.tsv is
written last on success; a directory without it is not an accepted input set.
flank=N (default0) prepares native boundary contexts, retaining original IDs and
adding source identity. Extract with CutGff flank=N gccontig=caller filename=f.
Validate with the same flank=N and cores=RF00174.fa,RF01482.fa,RF01689.fa, naming
the three accepted exact original FASTAs. Every de-flanked core must agree byte
for byte. Per-member .genomes.tsv/.loci.tsv then bind all original loci for refit;
missing TaxIDs stay missing, and unmeasured taxonomy/contig census stays UNKNOWN/.
No training, validation partition or endpoint acceptance is performed here.
member=RF00174 optionally selects exactly one producing-CM member (RF01482 or
RF01689 also supported). The complete ledger's IDs and ALL THREE expected counts
are validated first, so selection never renumbers originals or hides bad rows.
Only that member is projected/extracted/validated; no attribute substring filter.
For flanked selected validation, cores= names ONE accepted original-core FASTA.
Only the selected member's FASTA/learning files and context_geometry.tsv are
written. The latter reports observed physical flank clipping, separately from
later sampled candidate-span failures. Omitting member= preserves the old route.
HELP
 exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [[ -h $SCRIPT ]];do
 cd "$(dirname "$SCRIPT")"
 nextScript="$(readlink "$(basename "$SCRIPT")")"
 SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=4g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "${DIR}current/" prok.CobalaminConstructionInputs "$@"
