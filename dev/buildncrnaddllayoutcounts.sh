#!/bin/bash
# Native DDL objects and saved-window match counts; no fitting/caller.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/buildncrnaddllayoutcounts.sh \
  recipes=layouts.recipe.tsv panels=panels.tsv selection=selection.tsv \
  ref=consensus.fa members=members.fa ledger=ledger.tsv \
  assignments=assignments.tsv census=models.tsv summary=members.summary.tsv \
  out=fresh_prefix expectedmodels=508 expectedgroups=3153 t=1 -Xmx2g -eoom

Constructs every member, consensus and query as an existing native DDL object.
Each actual sketch part has a power-of-two width and uses native insertion.
wordtag=xor (default) preserves the original high-bit length tags. Optional
wordtag=leading uses bit2k above each raw word after middle masking. The explicit
legacy control always preserves its historical untagged/XOR forms.
Member winners vote once per nonempty bin; pooled consensus adds its declared
bonus to one native winner per bin. Shared DDL comparison counts equal nonempty
bins. Saves active counts then zero padding to four slots per window/model.
Consensus sketches are generated from existing consensus DNA, never stored as
resources. Existing loaders currently support the bound RF00174 data only.
No genomes, new oracle, alignment, fitting, caller integration or adoption.
study=replicates identities=prefix sourcegroups=accepted_groups.tsv.gz
genomes=original_donor_roster.tsv foldseed=20261114 ddlseed=12345
measures one complete five-fold assembly assignment at one hash seed. Original
member groups and exact sense-oriented duplicate locus sequences are excluded;
consensus DNA and the saved accepted-model oracle remain unchanged. Complete
genes come from saved windows,with original genome fallback only where needed.
Each fold rebuilds padding blacklists from fit windows only. The three word
pairs,three single-DDL layouts,blacklist off/on,and512/1024bins give36methods.
An additional bonus arm uses only unfiltered mixed7/masked9 at512bins:fixed
0,5,100,1000000 and scaled tenths0.1 through1.0 of actual remaining voters.
Scaled values round up with minimum1;fixed0stays0. No fit members leaves only
the specified consensus bonus;fixed0then gives an empty model sketch. Original
model ordinals are retained. Size thirds use remaining voters,then model ordinal.
Consensus retention counts matching NONZERO consensus keys divided by all bins.
Output includes out.replicates.tsv,out.folds.tsv,out.cohorts.tsv and out.models.tsv.
The latter tables retain the exact low-identity cut and remaining voter/size-third
assignments. Pool disjoint folds before statistics,using summarizencrnaddlreplicates.sh.
folds=5,shortprefixes=64,64,256,longprefix=4096 are the released design.
Smaller explicit values support finite witnesses;the production report rejects
incomplete folds or mismatched prefix metadata.
fixedcheck=t instead judges the original not-seen loci in their frozen two
blacklist-fit halves,with all original member voters retained. indexreference=t
also emits the unseeded original-consensus-index reference;use it in only one
hash-seed worker per assignment. These are saved-window studies,not caller tests.
study=gcdescribe describes every loaded model without held-outs. Writes
out.gc_models.tsv (member sequences,accepted padded windows,and each model's
own accepted spans) and out.gc_family.tsv (distribution of model means and
median within-model SD). GC uses Tools.calcGC: A+C+G+T denominator, ambiguity
excluded. Count is finite-GC observations; missing_gc is separate. Sample SD;
singletons have SD0; models without observations retain NA statistics.
Each accepted model/window pair counts once regardless of parent links.
No new sequence, alignment, ranking or caller run is performed.
study=gcalllosses fits GC once from ALL members and ALL accepted padded windows,
without ranking holdouts. The same member/family fallback and paired offset use
all data. Outputs out.gc_all_fit.tsv,out.gc_all_limits.tsv,out.gc_all_losses.tsv.
Limits retain maximum finite standardized distance,explicit infinite zero-SD
mismatches,the lossless SD width and lossless min/max margin. Loss curves use
z1..6 in quarter steps then6.5/7/7.5/8; margins0/.01/.02/.05; exact lossless
rows retain every accepted pair and original comparable locus. An infinite
lossless SD width is unavailable,not silently replaced by a finite width.
No DDL construction,hash-seed replication or fitting of held-out folds occurs.
study=replicategcall gcallfit=counts.gc_all_fit.tsv gcalllimits=counts.gc_all_limits.tsv
loads the same full-precision all-data statistics once for full or fixed ranking.
Neither ranking population nor any fold refits GC. DDL member and padding-word
exclusions stay unchanged. Each source/width uses weights0/.5/1/2/4/8;
z2/3/4/6/8 and its measured finite lossless z;margins0/.01/.02/.05.
These64methods reuse two model builds per fold. The shared all-data GC table
is referenced rather than copied into every worker's output. Both sources must
have a measured finite lossless width for this complete ranking roster.
gcadditions=t with study=replicategcall scores only six NEW window-statistics
combinations per bin count: observed ranges expanded by0/.01 combined with
penalty weights.5/1/2. The five previously proposed fine-SD rankings are cancelled. Gate models first,then rank survivors
by equal-bin count minus weight*absolute standardized GC distance. Missing GC,
ordinal ties,zero-SD bottom tier and failure rank retain existing semantics.
No old controls or SD4/6settings are rescored;the reporter imports saved rows.
study=gcallpower reads the same accepted gcallfit/gcalllimits without refitting
GC or building DDLs. Run once on the full configuration and once on the ORIGINAL
fixed selection with fixedcheck=t,not NOT_SEEN labels from the expanded panel.
Writes out.gc_power.tsv,out.gc_power_histogram.tsv,out.gc_survival.tsv.
Power reports both all saved windows and windows linked to comparable loci,
with explicit window denominators,mean/midpoint-median models removed,out of508,
and pair/window/locus losses. Locus survival combines acceptable models across
all its saved windows. SD1..8 by.25 plus literal8.58 and range margins0/.01/.02/.05.
The full invocation additionally writes out.gc_model_limits.tsv,
out.gc_model_quantiles.tsv and out.gc_model_histogram.tsv. Each model retains its
own accepted/finite/missing window counts,raw maximum standardized distance and
inclusive lossless SD width. Count strata<10,10..29,>=30 use accepted windows;
no-window and no-finite-GC models remain explicit. Distribution quantiles use
nearest rank;histograms have half-SD bins including zero frequencies,with
missing/infinite values separate. Literal8.58 is not the exact global lossless
8.583800604559325. These are saved-window diagnostics,not caller accuracy/time.
study=replicategc reuses the saved ranking oracle and fixed fold rosters. Uses
8/masked10 quarter-short, shortprefixes' third value and longprefix, with bonus60
at512bins and100at1024. Production lists are256/2048. Twelve settings per width:
equal-bin count minus w*abs(windowGC-fitMeanGC)/fitSD,w=0,0.5,1,2,4,8;
mean plus/minus z*fitSD filters,z=2,3,4;min/max plus/minus m,m=0,0.02,0.05.
Compare two GC sources: accepted fit-side padded windows, and surviving fit-side
member sequences. Missing window statistics use member mean/range shifted by
the family mean of paired window-minus-member model means, retaining member SD.
Models with no surviving members use fit-family mean and median within-model SD
and its observed range; no held-out member statistic is used. Missing query GC
leaves models unpenalized/unfiltered. Zero SD gives exact-mean
zero penalty; a mismatch sorts below every finite score for positive weights.
w0 bypasses GC. All ties keep original ascending model ordinal. A filter losing
every acceptable model gets rank models+1 (509 in production), an explicit failure
retained in the original judging population and percentile distribution.
Writes out.gc_fit.tsv and the usual fold/model metadata. Replicate rows add nine
window/model-pair loss and fallback counts for ALL only, with histograms including the
failure bin. Each judged window/model pair is counted once before parent links.
All24alternatives per width reuse the same native scores, not separate model builds.
study=gclosses gcfits=accepted_workers.tsv reuses the recorded seed12345 folds
once per assignment. It builds no DDLs and scores no sketches. Outputs
out.gc_curve.tsv and out.gc_choices.tsv: widths1..6 by0.25,then6.5,7,7.5,8;
both statistics sources; minimum observations0or30 from that source alone.
Fits reuse the exact sequence accumulation,never rounded GC summary text.
Fit-fold lossless widths are selected before inspecting judging losses;
unavailable through8 is NA. Pooled judging widths losing at most1,5,10loci
are descriptive tuning only. Exact pair/window/locus loss parity with the
accepted z2/3/4 outputs is required for both sources and both bin counts.
study=replicatefollowup preserves the same data,folds and exclusions but measures
only20settings:8/masked10 quarter/padding at512/1024 with fixed60/80/100/120,and
8/masked11 at512/1024 across mixed/half/quarter and plain/padding with bonus100.
The8/masked11 padding list uses the third shortprefix and the same longprefix.
Only the two expressly repeated100controls overlap the prior study. Each row
adds exact rank_histogram and low_rank_histogram fields,one comma-separated count
per original model rank,starting at rank1. Pool histograms before percentiles.
study=replicatesizes holds8/masked10,quarter1024,bonus100 and tests padding lists:
short0/16/32/64/128/256/512/1024 by long0/256/512/1024/2048/4096/8192/16384/32768.
Its71new rows omit256/4096,which is imported once from replicatefollowup.
The0/0row repeats the old unfiltered winner for percentile collection.
study=replicatedefinitions runs false-seeded-window,padding+false-pooled and
positive-member lists at shortprefixes' third value and longprefix. The latter
counts one presence per model using only members surviving the fold exclusions.
Use the previously reported one-paired-SE size choice,not a new manual cutoff.
study=partitionseeds measures one512-bin DDL with256short/256long bins or
128short/384long bins on the same eight seeds as study=seeds. One native
comparison scores the entire array. Post-hash selector bits separate lengths;
XOR pre-hash tags remain identical to mixed. Long words use three quarters
in the unequal case,uniformly routed by unused hash bits9..35.
study=partitionblacklist uses
the one explicit ddlseed and both opposite judging halves with the fixed
padding-derived64short/4096long prefix. Merge its four raw rows with mergestudy=.
These modes write new partitioned rows only;reuse completed mixed controls.
study=bestblacklist selects the paired-mean short8/masked10,mixed512,pooled100
recipe and measures only padding-derived prefixes:64/256short by1024/4096long,
all four combined settings,both opposite judging halves at one explicit ddlseed.
Lists are learned from that word form and fit half;the accepted spans and
frozen group halves are unchanged. Masked10GC excludes its masked base.

Recipe TSV has this exact header and one row per layout:
layout_id short_k long_k long_mask member_mode member_bins member_short_bins consensus_mode consensus_bins consensus_bonus exponent_bits mantissa_bits hash_namespace legacy_control query_form
Use TABS, not spaces. Lengths are explicit, with 1<=short_k<long_k<=15.
long_mask: MIDDLE_RIGHT (one base) or NONE. Member mode: MIXED (short bins zero)
or PARTITIONED (explicit short-bin boundary). Key fields are 4 and 12; exponents
WRAP identically for members and queries. No origin/support arrays are stored.
Consensus BOTH/LONG_ONLY uses separate bins and zero bonus. POOL_BOTH/POOL_LONG
uses zero separate bins and a positive bonus. NONE uses zero bins and zero bonus.
Query form: QUERY_DDL only. No raw-query scoring or custom bin mapping.

Declare exactly one legacy_control=true row: short7/long10/MIDDLE_RIGHT,
PARTITIONED/member384/short128/BOTH/consensus128/bonus0/4/12,
LEGACY_PARTITIONS/true/QUERY_DDL. It preserves the original word hashes and
smaller-key vote ties, using three separate native parts. New rows use
NATIVE_DDL/false and greater-key vote ties. Use at most100 unique recipe IDs.

Writes out.layouts.tsv (version-two,31fields) and out.counts.tsv.gz
(layout/window/panel/input/event then four slots per model).
Active terms are member_total for MIXED or member_short/member_long for
PARTITIONED, then consensus_total only when separate consensus bins exist.
Stored member keys and runtime consensus are separate byte counts; no index is
built, and all index fields are zero. Payload bytes are not heap estimates.
Use optimizencrnapartitionweights.sh layouts=... counts=... control=... to fit
weights later against the original oracle/cohort bindings. Output must be fresh.

Optional standalone index validation: checkindex=exact_single_DDL_recipe_id
selects that one recipe from the roster and uses actual DynamicDemiLog objects
for both index references and queries. Checks every saved window/model count
against direct comparison, including zeros. Writes out.layouts.tsv and
out.index.tsv instead of a count cache. Prebuilds native objects before timing.
Layout metadata describes the retained DDL model keys; transient diagnostic
index payload bytes and pointer slots are reported separately in out.index.tsv.
indexseconds=10 (default) times whole query sweeps independently for the small
index and direct comparison until each consumes at least10main-thread CPU
seconds. Reuses query scratch and consumes all returned counts in a checksum;
construction, hashing and I/O are excluded. Reports actual calls and CPU/window.
indexseconds=0 is finite correctness only: timing is NA. Values1..9 are rejected.
No new index representation, fitting, identities, alignment or caller run.

Fixed winning-recipe study: study=baseline or study=blacklist identities=prefix.
Requires native_s7_l9_mixed512_pooled_both_w100 in the roster. Uses only that
forward7/masked9 mixed512 recipe with pooled-both consensus bonus100 and unit
scores. Reads the frozen identity/cohort files; never recomputes identities.
Baseline writes out.study.tsv with the selected tag's result. Blacklist mode writes
168 rows: one library-presence list and three background lists in each of two
opposite fit/judge directions. Each has four short-only,four long-only and
sixteen joint prefixes. Short sizes16,64,256,1024; long64,256,1024,4096.
MODEL_PRESENCE counts one vote per model across ALL its construction members.
FALSE_SEEDED_WINDOWS counts occurrences with zero fully accepted models in the
complete cached oracle. TRUE_WINDOW_PADDING uses only flanks outside the envelope
of all accepted spans; POOLED_BACKGROUND combines those two separate corpora.
Background lists use only fit-half groups; ranks/nonzero means use the opposite
half. Frozen target halves are preserved; additional background groups receive
a separate deterministic split in out.background_groups.tsv. Top-prefix GC
excludes the masked base and is reported beside construction-member base GC.
Both word forms are rejected before hashing, identically in members,
consensuses and queries. Lists use frequency descending/raw-word ascending ties;
out.blacklist.tsv retains the maximum prefix for each definition/form. This is
exploratory saved-window ranking, not caller accuracy or independent validation.
Reports original not-seen/low populations, rank sums/top10/top50, and nonzero
model-score pairs per window; no count cache, weight search or new index.
Absent blacklist GC is NA. Report-only mergestudy=report.study.tsv out=prefix
sums the two disjoint judging halves per setting,weights nonzero means by
window counts,and preserves both fitted lists' GC. Requires each half exactly
once; library rows appear once. Writes out.merged.tsv without loading biology.
study=seeds runs the fixed winner for both historical high-bit and leading-bit
tags at seeds12345,1,7,42,123,1337,20261114,8675309. Writes sixteen rows with
an additional ddl_seed column and out.seeds.summary.tsv (mean/min/max of six
rank/top-k metrics per tag). ddlseed=12345 (default) sets BOTH explicit factory
arguments and process default for other modes;the paired study uses its listed
seeds. It measures hash-choice variability,not an additional recipe search.
Report-only mergeseeds=roster.tsv seedbaseline=paired.study.tsv out=prefix
For8/10 filter-only comparisons,seedbaseline may instead be the completed
recipe-comparison roster (ddl_seed,comparison_table). It selects the exact
8/masked10 mixed512 pooled100 row,checks original seed/cohort bindings,and
reports the six ranking metrics. It does not invent a nonzero-window baseline.
reads a ddl_seed/merged_table roster (absolute paths;producer ddl_seed.txt
sidecar required). Reports mean/min/max and wins/ties/losses against each
matched unfiltered XOR seed,for rank/top-k/nonzero metrics. Saved tables only.
Seed studies also save every not-seen locus rank,rank ranges ordered largest
first,and an exact tail decomposition between the lowest/highest sum seeds.
study=variance1024 is the sole requested width control:the current fixed
7/masked9/pooled100 recipe,XOR only,on seeds12345,1,7,42 at1024bins.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.BuildNcrnaDdlLayoutCounts "$@"
