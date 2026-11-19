#!/bin/bash
# Join accepted member/native/assembly evidence without reading sequence payloads.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_construction_members.sh in=groups.inputs.tsv.gz genomes=corpus.genomes.tsv bindings=original_prefix evidence=evidence_prefix out=fresh_prefix expectedgenomes=39979 expectedmembers=26896 expectedeval=11769 expectednew=29313 t=1 bgzip=t pigz=f zl=9 -Xmx2g -eoom
Retains every old member and native FMN event, with independent admission and
assembly-exclusion decisions. One distinct clean full core may replace old bounds
only under the existing prior/own-hit reciprocal-half correspondence rule. Weak
overlaps, multiple distinct cores and multipart associations retain priors for review;
unsupported old members remain labelled priors. Old partial/multipart/ambiguous sequences
remain seed-only priors. Competing full other-family evidence requires review.
Exact duplicate native cores share one extraction identity; every witness remains.
Evaluation assemblies and known aliases never train. Unresolved assemblies remain
separate. Reports phyla without an intact training representative; never silently
creates a coverage exception or filters identical genes from different genomes.
Writes .old.tsv.gz, .events.tsv.gz, .extract.tsv.gz, .coverage.tsv and .summary.tsv.
This is a metadata ledger requiring review and sequence/flank validation, not a
frozen split, completed extraction or resource build. Existing outputs fail.

Coverage audit of an accepted saved ledger (no admission replay or sequence reads):
mode=coverage in=members_prefix groups=groups.inputs.tsv.gz genomes=corpus.genomes.tsv
out=fresh_prefix expectedgenomes=39979 expectedmembers=26896 expectedeval=11769
expectednew=29313 expectedcores=27631
Writes .coverage.tsv, .unlabelled.tsv.gz and .summary.tsv. Counts all scanned
inputs, training inputs, inputs with intact/native-clean cores, old members and
old-positive original genomes by the original domain/phylum labels. Distinct
intact assembly groups are separate columns. Includes zero-positive categories.
Missing literal-roster labels remain UNKNOWN with every affected input listed;
an absent literal-path binding is not proof of an unseen assembly. Does not
inherit taxonomy through aliases, modify the split or create coverage exceptions.

Stage audit: mode=stages uses the coverage inputs plus evidence=joined_FMN_prefix,
imported=raw_import_prefix, raw=collection.tsv, models=scan_models.tsv,
expectedbatches=641 expectedraw=2214420. Reads saved raw tblout/maps/completions,
imported hit rows, complete joined FMN events and accepted admission outputs.
Writes .stages.tsv, .examples.tsv.gz and .summary.tsv. Counts raw discovery events,
imported events, native GA/full/CONFIDENT events, admitted distinct cores and
old-only intact priors separately. Evaluation raw/import counts are NA because
only its saved classified evidence is inspected. Validates every raw model's
completion count, including zero-hit FMN; includes RF00010 and RF00174 raw controls.
Retains up to five raw dropping examples per taxonomy category and reason.

Crossed-search selection: mode=cross-select uses the coverage inputs plus
taxonomy=original_GTDB_accession_taxonomy.tsv.gz crosscount=20.
Selects distinct Bacteroidota assembly groups per cohort in deterministic
class/order rounds. Evaluation candidates require an admitted native FMN core;
discovery candidates include zero-hit inputs. Full taxonomy binds the original
roster's assembly accession and must agree with its domain/phylum; missing
accessions or ranks remain explicit. No sequences or searches are run.
Writes .selected.tsv, .composition.tsv, .expected.tsv, .summary.tsv and
.evaluation.paths.tsv/.discovery.paths.tsv. Path files feed PrepareCmSearchBatches
include= directly; expected loci retain complete admitted metadata. This is a
small targeted procedure diagnostic, not a training split or accuracy sample.

Family/genus comparison: mode=lineages uses the coverage inputs plus taxonomy=.
Writes .lineages.tsv, .lineage_inputs.tsv.gz and .summary.tsv. Compares all
Bacteroidota discovery inputs with evaluation native-FMN-positive assembly groups,
also retaining every labelled evaluation input, including old-only and zero-core
inputs. Family/genus counts deduplicate resolved assembly groups; missing ranks
remain UNKNOWN, not absence. Original accession labels must agree across aliases.
This metadata-only census neither changes the split nor reads genome sequences.

Paired procedures on the SAME selected discovery inputs:
mode=procedure-report in=genus_audit.selection.tsv raw=procedure_collection.tsv
models=original_queries.tsv out=fresh_prefix crosscount=49.
Collection columns: procedure,input_id,prepared,native. Procedure values are
DISCOVERY_BUNDLE and EVALUATION_FMN. Every selected discovery input must have
both complete native results; cohort identity is not relabelled as procedure.
Reuses the crossed native target/model/count validator and emits .hits.tsv,
.inputs.tsv and .summary.tsv, retaining sub-GA reporting separately from GA.

Crossed-search accounting: mode=cross-report in=selection_prefix
raw=cross_collection.tsv models=original_queries.tsv out=fresh_prefix crosscount=20.
The collection columns are cohort,input_id,prepared,native (tab-delimited).
Each point must be one complete selected genome in the native batch-preparation
format. Evaluation inputs use the full discovery model roster; discovery inputs
use the single FMN evaluation procedure. Checks successful exits, complete target
maps, all model identities/dimensions and reported/raw hit conservation. Writes
.inputs.tsv, .hits.tsv, .expected.tsv and .summary.tsv. Native GA and full-model
counts remain separate from lower-threshold reporting. Expected-locus flags mean
same-contig/strand reciprocal-half presence, not one-to-one TP/FP/FN grading.
No new search, sequence extraction, resource generation or split change.
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
parseJavaArgs --xmx=2g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.NcrnaConstructionMembers "$@"
