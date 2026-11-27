#!/bin/bash
# Build source-bound cobalamin profile resources before caller experiments.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/buildncrnaddlprofiles.sh in=members.fa ledger=routing.tsv
  assignments=assignments.tsv census=membership.models.tsv
  summary=membership.summary.tsv ref=consensus.fa out=fresh_profiles.tsv
  expectedmodels=508 expectedmembers=85726 dual=f partitionbins=0 t=1 -Xmx2g -eoom

Joins the original frozen construction members, validates retained voter counts,
and votes one forward middle-masked k10/2048 profile per current model.
dual=t instead builds the fixed partitioned plain9+masked10/4096 resource,
scored as long + 0.5*short in sketch-only caller order.
Without selected=t, partitionbins=512 or 1024 builds the historical plain-7 member quarter,
masked-10 member half and shared consensus quarter in one sketch-only index.
It cannot be combined with dual=t. Other geometries and weights are not enabled
by this builder. Undefined sequence bases break words.
Seed12345/exponent4 and raw-query-key scoring are fixed by the versioned format.
The resource binds ordered model IDs and uppercase consensus bytes; changing
the consensus library requires rebuilding it. Voting occurs once, not per
worker or genome. No genome, verification, alignment or caller invocation.
This experimental resource does not enable a caller or change its defaults.

Append existing ordering profiles without voting:
  appendprofiles=t in=old.profile.tsv in2=new.profile.tsv
  ref=old.consensus.fa ref2=new.consensus.fa expectedmodels=11,2 out=fresh.tsv
Both complete model axes must be disjoint and from the same member, with equal
geometry, seed, padding set and vote policy. All old metadata and profile rows
are preserved; only new row ordinals are shifted. Output is decoded and checked
against every source consensus, voter count and key. No construction/holdout/
voting flags are accepted in this explicit mode. It never rebuilds old profiles.

New all-gene mode:
  selected=t in=construction_members.fa assignments=builder.assignments.tsv
  ref=new_consensus.fa expectedinput=N expectedmembers=KEPT expectedmodels=M
  partitionbins=1024 paddingblacklist=none out=fresh_profiles.tsv
Omit the historical ledger/census/summary inputs. Full FASTA headers bind each
actual new builder assignment, including nonvoting small/orphan records.
One native maximum key per member/bin votes, with100extra votes for every
distinct eligible consensus bin/key; ties use the greatest key. The version4
resource fixes forward8/masked10 quarter-short routing and native4+12 encoding.
Queries count matching bin maxima, preserving a caller's global DDL exponent.
Use explicit paddingblacklist=none initially, or a complete256short/2048long
table plus paddingfit=scope later. The resource retains those words so query
and model filtering agree. Both512 and1024 total bins are supported explicitly.

Ranking-only main-member comparison:
  orderonly=t in=original_members.fa assignments=original_assignments.tsv
  ref=current_main_consensus.fa expectedinput=91072 expectedmodels=354
  manifest=formal_manifest.tsv panel=panel1000.sources.tsv expectedpanel=1000
  sources=source_ids.tsv aliases=resolved_assembly_groups.tsv.gz
  partitionbins=1024 paddingblacklist=training.padding.tsv paddingfit=TRAINING
  seeds=12345,1,7,42,123,1337,20261114,8675309 out=fresh_profiles
The current retained models receive only their original construction members;
removed clusters are not reassigned. Observe all panel/alias original sequences
before excluding exact literal DNA duplicates (case folded, U=T; IUPAC symbols
match only the same symbol). A zero-voter model keeps an empty member sketch.
No consensus keys or bonus are used. The resource carries the exact seed and
uses native DDL objects; cobalaminddl=that_file changes only native-eligible
order before the ordinary cap/patience/verifier. Both minor families stay native.

Explicit in-profile bonus experiment (orderonly=t only):
  consensusbonus=100 consensusrule=maxima partitionbins=1024
  consensusbonus=60 consensusrule=maxima partitionbins=512
Defaults remain consensusbonus=0 consensusrule=none. A positive bonus requires
an explicit rule. maxima adds that many votes to each nonempty native consensus
bin maximum; allkeys adds them once to every distinct eligible consensus bin/key.
Those are different experiments. Real members still contribute one maximum vote
per bin; plurality plus the bonus selects one key, breaking ties by greatest key.
All profiles retain their REAL kept-voter counts. A zero-real-voter bonus profile
contains only the actual consensus maxima, never an invented training member.
A distinct versioned ordering-bonus header records the amount/rule. The ordinary
member-only format still requires zero-voter models to be empty. Suppression and
separate runtime consensus scoring cannot be combined with an in-profile bonus.
Joint loading accepts512or1024bins only when every member has identical width,
seed, blacklist and vote policy. No caller settings or resource defaults change.

Optional suppressconsensus=t (orderonly=t, partitionbins=1024) removes a member's
per-bin maximum vote when it equals that model's own consensus maximum. Choose
the best-supported remaining key, breaking support ties by greatest key; no
second word from a voter replaces the suppressed maximum. The explicit new
format stores only those member keys. The shared clan loader builds a separate
consensus index with the identical seed, geometry and padding blacklist, checks
binwise disjointness, and sums both indexes with at most one hit per query bin.
All three member resources must use the same mode. This experimental format
cannot be loaded by the separate-member caller.
Multiple seeds write out.SEED.tsv; one seed writes out exactly. All outputs
are fresh, and out.members.tsv records assigned/kept/excluded model counts.

Minor-member ordering profiles use the same orderonly=t geometry and exclusions.
Bind ref= to one current RF01482 or RF01689 library and expectedinput= to that
member's complete original FASTA. Supply additions=expanded.additions.tsv when
the current library retains explicit expansion singletons. Original cluster
assignments remain intact: a source can vote in its cluster and its singleton.
An explicit protectedpriors=accepted_companion.tsv expectedpriors=N supplement
admits only the companion's exact historical IDs/headers. Supply the augmented
full original FASTA and assignment overlay, counting priors in expectedinput.
The native source ledger stays unchanged. Priors must bind the same member,
formal assembly/input and canonical original FASTA; all native and historical
genes are observed before the same alias and literal-DNA holdout decisions.
This mode preserves the member-only profile format and append compatibility.
Removed singleton models receive no vote. Each original is observed and tested
for holdout once; every excluded source loses all its vote edges. Minor resources
declare their actual member. out.members.tsv counts model-vote edges, while
out.dispositions.tsv separates unique original-gene dispositions from additional
singleton edges. out.holdout.tsv and out.holdout.tsv.genomes.tsv audit the member's
source-bound exclusions. Zero-voter models retain empty profiles.

Release refit after methodology is fixed: add refit=t to orderonly=t and omit
panel=. No positive expectedpanel is required. The same checked formal manifest,
source ledger, aliases, assignments and optional singleton/prior joins remain;
panel membership, aliases and exact duplicate sequences exclude no eligible
voters. Removed/small/orphan assignments are not reassigned. Geometry, seed and
the explicitly supplied padding blacklist stay unchanged. Default refit=f keeps
the development exclusions. Refits write out.refit.tsv plus its genome audit,
out.members.tsv and out.dispositions.tsv; none is a held-out evaluation.
refit=t cannot be combined with prepareonly=t or a supplied panel=.

Before padding exists, use orderonly=t prepareonly=t without paddingblacklist
or paddingfit. With trainsample=4000 sampleseed=20261119 it writes out.training.tsv,
out.holdout.tsv and out.holdout.tsv.genomes.tsv instead of profiles. Sampling
follows the formal-roster order shuffled by the declared Java Random seed,
after panel/alias and duplicate-voter-genome exclusions. This conservative
sampling rule does not remove other distinct member votes from those genomes.
Historical alias misses remain explicit in the per-genome audit; canonical
FASTA and exact accession identities also join the current roster.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.BuildNcrnaDdlProfiles "$@"
