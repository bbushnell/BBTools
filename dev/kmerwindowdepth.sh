#!/bin/bash
# Exact unfiltered read counts along forward-reference windows; one table per run.
# in=reads ref=fasta scaffold=record_name regions=first-last[,first-last] out=tsv k=N
# Region positions are ONE-BASED INCLUSIVE KMER ENDS, not windows of kmer starts.
# The named reference is query-only. It never contributes to the read counts.
# Outputs the central depth and A/C/G/T extensions on each side, in reference orientation.
# hollow=f by default; hollow=t uses two32bp anchors with paired hashes at packed K>64.
# queryspan=N (default k) keeps observed endpoint anchors fixed while testing another table k.
# Unequal spans require hollow=t and both spans>=66. Extra query_span column then distinguishes
# the observed word/coordinates from table k. First/last33 query bases preserve neighbor counts.
# Synthetic middle bases remain ignored. No extension decisions are changed.
# queries=seeds.fa instead of ref/scaffold/regions streams one terminal query per seed.
# Requires unpaired defined-base seeds at least queryspan long; emits ID/length/k/span/anchors/depth.
# querystart=0: with queries=, a positive value switches to batch window profiles.
# Report each context ending at querystart through length-1, immediately before its
# saved next base. Coordinates are 1-based inclusive; querystart>=queryspan is required.
# One table load serves all records. A record exactly querystart long emits no rows.
# Columns: exact ID, k, query_span, start1, end1, observed word, depth, four left
# counts, four right counts, query_previous and query_next. Counts use query orientation.
# Requires unpaired records of length>=querystart, tab-free IDs and defined context/
# suffix bases. All path policies and outseq must be off. No truth or base changes.
# Experimental path veto: queries=saved_extensions.fa pathseedlen=N outseq=guarded.fa out=decisions.tsv
# All query records must share the same original seed length N. Only right-added bases are tested.
# Requires hollow=t, k=queryspan+/-1 and both spans>=66. Default: pathseedlen=0 (off).
# pathsupport=3 pathratio=4: stop before a saved next base when the alternate-span predecessor
# has support>=3 and another next base has support>=3 and>=4*the saved base's support.
# Never corrects seeds or switches paths. Query reads never add counts; in/extra remain the evidence.
# One report row per read; evidence fields are -1/dot when no veto fires. Outputs must be new.
# pathbothspans=f: set t to test observed spans k+1 and k-1 with one loaded table.
# Requires pathseedlen>=k+1, k>=67, and queryspan omitted or equal to k.
# Stops at the earliest conflict; k+1 wins ties. Report query_span identifies the vetoing
# spacing, or -1 if neither vetoes. Single-span reports retain their previous format.
# pathrollback=f: set t to remember the earliest ambiguous added-base choice per spacing.
# Both choices and their context must meet pathsupport, and neither choice may dominate
# by pathratio. Ambiguity alone retains the path. A later hard veto in the SAME spacing
# rolls back to that checkpoint, never into the original seed. Adds conflict_end to the
# report (prefix length at the hard contradiction, -1 if absent); counts describe that
# contradiction, while retained_length describes the possibly earlier rollback point.
# pathbranchstop=f: separate experiment that stops immediately at supported ambiguity,
# using the same context/choice support floors and strict ratio comparisons above.
# Hard-conflict stops remain enabled. Does not require a later contradiction and can
# therefore discard correct extension. Requires positive pathseedlen and pathrollback=f.
# The report counts describe the stopping branch; no conflict_end column is added.
# pathselfstop=f: additive diagnostic veto when a nominal next key has positive
# global count entirely accounted for by the ORIGINAL seed's own multiplicity.
# Requires selfseeds=original.fa counted exactly once in in/extra, no in2, and
# fixed pathseedlen>=k. Both files must be FASTA. Originals must be unpaired, in query order, and match
# every saved prefix. Forces unfiltered/untrimmed/unmerged counting and all reads.
# The existing alternate-span veto is unchanged. No shared count is subtracted;
# the seed ledger only determines whether to stop earlier, never another path.
# Adds self_count/external_count columns (-1 when this veto did not fire).
# A self-only stop reports nominal query_span=k, raw context/chosen counts,
# other_count=-1/other_base=., and external_count=0. Not a live FASTQ/paired policy.
# pathbackstop=f: optional nominal-spacing veto when another backward neighbor has
# support>=pathsupport and>=pathbackratio (default2) times the ACTUAL predecessor.
# The current nominal word AND actual predecessor must also meet pathsupport.
# A weak actual predecessor may reflect an existing seed error rather than an
# unsafe extension. Balanced backward branches are not rejected by this check.
# pathbackactualsupport=0 inherits pathsupport after all arguments are parsed.
# A positive value independently sets the actual-predecessor floor; context and
# competitor still require pathsupport. Used only with pathbackstop=t. Setting1
# tests weak observed predecessors, but may also discard valid seed-error suffixes.
# Zero/default preserves the existing behavior; negative values are invalid.
# chosen_count/other_count then describe BACKWARD
# neighbors; other_base uses the saved sequence orientation, not canonical orientation.
# pathdropstop=f: optional nominal-spacing veto when context>=pathdropcontext (32),
# the saved next count<=pathdropmax (2), and context>=pathdropratio (16)*next_count.
# This detects a transition, not every low-depth path; a supported-to-zero drop counts.
# Both flags require pathseedlen>=k+1 and ordinary path-veto mode. They only shorten
# saved extensions, never edit seeds or affect Tadpole/BBMerge defaults. Thresholds
# are experimental. pathbackratio/pathdropratio>=2, pathdropcontext>=1, pathdropmax>=0.
# Enabling either flag adds stop_reason: none, alternate_span, self_only,
# backward_dominance, depth_drop. Nominal stops report query_span=k; depth_drop has
# other_count=-1/other_base=.. Existing stops win equal-position ties; backward
# dominance wins a nominal tie. With both flags off, the existing TSV is unchanged.
# pathrisestop=f: optional nominal-spacing veto before a sharp increase into a
# deeper endpoint-key class: current>=pathsupport (3), next>=pathrisemin (32),
# next>=pathriseratio (16)*current. Requires positive pathseedlen>=k+1 and ordinary
# path-veto mode. Ratios must be>=2 and pathrisemin>=1. Weak current keys abstain
# because an existing seed error leaving an anchor can also cause a depth rise.
# Stops before the transition; never edits the seed or selects an alternative.
# Enabling this flag also adds stop_reason, with depth_rise for this veto.
# Existing stops win ties, followed by backward dominance, depth drop, then rise.
# This experimental flag does not change live Tadpole/BBMerge behavior or defaults.
# pathlowbranchstop=f: separate nominal forward check for minor counts2 or3 only.
# Context and saved next count must both meet pathsupport; stop when the saved
# next count is STRICTLY LESS than pathlowbranchratio (32) times the strongest
# other count. Ratio equality retains the path. Minor counts0/1 or>=4 abstain;
# this is not a general junction classifier and does not replace the extender.
# Requires ordinary path-veto mode, pathseedlen>=k+1 and pathlowbranchratio>=2.
# Adds stop_reason=low_forward_branch and reports nominal k plus forward counts.
# Existing stops win ties, then backward dominance, drop, rise, and this check.
# Never edits a seed, adds evidence or selects another branch. Default off;
# backward support and live Tadpole/BBMerge behavior remain unchanged.
# pathsolidstop=f: separate nominal SOLID-table contradiction test on saved paths.
# Requires hollow=f, queryspan=k, pathseedlen>=k, queries= and outseq=; all other
# path policies must be off. Stop only when the current solid count and another
# next-base count are both>=pathsolidsupport (2), and the saved next count is zero.
# pathsolidsupport must be>=2. Missing context/alternative evidence abstains, as
# does ANY positive saved-next count; this does not require solid coverage.
# Never changes the original seed or chooses another base. Report stop_reason
# is solid_contradiction, with nominal k and forward evidence; none otherwise.
# Inputs remain the original evidence, not saved extensions or a truth reference.
# This default-off diagnostic does not change live Tadpole/BBMerge behavior.
# pathsolidbranchstop=f: modifier requiring pathsolidstop=t and pathratio>=2.
# Also stop at supported ambiguity: context and BOTH next choices must each
# meet pathsolidsupport, and neither choice may reach pathratio times the other.
# Use pathsolidsupport=3 pathratio=4 for the fixed ambiguity experiment.
# Missing/weak context or choice support abstains; exact ratio equality retains.
# Retains zero-choice contradiction stops and reports new stops as solid_ambiguity.
# Reuses the same Solid table; no new evidence pass, replacement base or seed edit.
# Scans the saved suffix only. With this modifier off, old reports are unchanged.
# Diagnostic only; no live Tadpole/BBMerge or default changes.
# pathnominalstop=f: separate nominal HOLLOW-table contradiction test.
# Requires hollow=t, k>64, queryspan=k, pathseedlen>=k, queries= and outseq=;
# every other path policy must be off. Stop before the first saved next base
# whose count is zero while context and strongest alternative both meet
# pathsupport (3, must be>=2). Any positive chosen count or missing support
# abstains. Continue scanning after gaps; do not require blanket coverage.
# Reports nominal_hollow_contradiction with nominal k and forward counts.
# Counts come only from original in/extra evidence; saved queries never add
# counts. Use a separately loaded span to diagnose a saved extension from
# another span. Never selects another base, edits seeds, or changes live
# Tadpole/BBMerge behavior. Cross-span benefit and false cuts require testing.
# pathspacingabstain=f: modifier for zero-choice nominal Hollow contradictions.
# Before cutting, query the saved context and next base at observed spans k-1
# and k+1 against the SAME loaded table. If either context AND saved-next count
# meet pathsupport, retain this choice and continue scanning later positions.
# Requires pathnominalstop=t, k>=67, pathseedlen>=k+1, pathnominalweak=f and
# pathnominalrise=f. Source qualification, if requested, remains unchanged.
# No extra count table, seed correction or alternative-base selection. The
# report schema and nominal evidence at an eventual cut are unchanged. With
# the modifier off, no alternate-spacing scratch or queries are added.
# Experimental protection against seed-indel spacing failures; whole-library
# validation is required. No live Tadpole/BBMerge behavior or default change.
# pathshortreads=<original.fq.gz>: optional read-only Solid protection for
# zero-choice nominal Hollow cuts. Disabled when omitted. Requires pathnominalstop=t,
# pathnominalweak=f, pathnominalrise=f, pathsourcek=0. The FASTQ must already occur
# exactly once in in/extra and must not be the saved query or seed FASTA; no in2.
# Counts ALL of that FASTQ untrimmed/unmerged, with no probability filtering;
# seeds and query sequences never contribute to this second table. One FASTQ only.
# pathshortk=64: exact support word length, either64 (default) or96; must fit in k.
# Nondefault96 requires pathshortreads and pathshortprefilter=f. Other short-table
# diagnostics remain Solid64. At96, appends short_k to every report row; default64
# preserves the existing schema. Longer support costs more heap and needs validation.
# A nominal cut is allowed only if a contained short word in its preceding K-base
# context has depth0/1 OR the proposed Solid next key has depth0/1. Otherwise
# retain the choice and continue. Optional k+/-1 spacing protection still applies.
# Uses a rolling K-pathshortk+1-window summary, reusable words and one additional count
# lookup per accepted base; next-key lookup happens only for candidate cuts.
# Adds short_low_windows/short_chosen report fields (-1 unless a cut occurs).
# Two tables remain resident: budget additional heap and measure the cost.
# Diagnostic only; no seed correction, branch replacement or live/default change.
# pathshortprefilter=f: optional two-pass singleton filtering ONLY for pathshortreads at64.
# Uses prefilter=1, prepasses=1, onepass=f; the nominal Hollow table stays unfiltered.
# Reloads all original reads to retain exact counts>=2. Removed singleton keys
# become zero, preserving the weak depth0/1 predicate; collisions may retain some
# singleton keys without changing that predicate. Sketch sampling is not allowed.
# Renames short_chosen to short_chosen_filtered: zero now means original0 OR1,
# one means a retained singleton, and counts>=2 remain exact. The -1 sentinel
# still means no cut. short_low_windows counts the same weak context windows.
# Adds one input pass and may reduce table memory; measure net cost on the input.
# Default off; no change to live extension, seed correction or branch choices.
# pathnominalweak=f: modifier requiring pathnominalstop=t and pathratio>=2.
# Also stop at a count1 saved choice when context and strongest alternative meet
# pathsupport, and alternative>=pathratio*chosen (default ratio4, support3).
# Count2+ choices abstain; zero-choice contradictions retain their old behavior.
# Reports nominal_hollow_weak_choice; source qualification still applies when
# pathsourcek is positive. Declined candidates continue scanning later positions.
# No extra table, seed edit or replacement choice. Default off; requires separate
# whole-library safety/yield testing and does not affect live extension defaults.
# pathnominalrise=f: modifier requiring pathnominalstop=t. Also stop before a
# nominal Hollow depth rise: context>=pathsupport (3), chosen>=pathrisemin (32),
# chosen>=pathriseratio (16)*context. Reuses that immutable table; no third table.
# Scans the entire saved suffix, including after gaps or declined contradictions.
# Reports nominal_hollow_rise, context/chosen counts, and other_count=-1/base=.
# By default a rise is independent of pathsourcek: that comparison qualifies only
# zero-choice contradictions. source_chosen is -1 on a rise because no source
# comparison was performed. Preserves seed/quality prefixes and never chooses
# an alternative. This default-off experiment requires separate whole-library
# validation at each auxiliary span; no live Tadpole/BBMerge/default change.
# pathrisesource=f: modifier requiring pathnominalrise=t and positive pathsourcek.
# Qualify a nominal rise only when its saved NEXT key has source-span count>0
# and<pathsupport. Source context depth is not used: a weak context can recover
# to a strongly supported correct choice. Absent or supported choices abstain
# and scanning continues, so later rises or contradictions can still stop.
# Reuses the second table; no third table or new threshold. On an accepted rise,
# source_chosen reports the positive weak count. Other report fields/reason and
# seed/quality prefixes are unchanged. With this modifier off, old behavior and
# report bytes are preserved. This experiment needs broad safety/yield testing.
# pathsourcek=0: optional shorter Hollow source span for pathnominalstop=t.
# A positive value requires 64<pathsourcek<k and loads a second immutable table
# from the SAME original in/extra inputs. Only after a nominal contradiction,
# count the saved next key at this source span. Veto only when source_chosen>0
# and auxiliary other_count>source_chosen. Equal counts abstain: different spans
# may support different repeat copies. Otherwise continue scanning later
# positions; no stop at all preserves the original saved output. This does not
# select an alternative or edit the seed. Appends source_k/source_chosen (-1
# chosen when no source comparison supports a veto, including an unqualified rise).
# Zero retains previous output and one-table cost.
# This raw-count comparison is experimental; differing spans have different
# coverage. Budget additional heap for the second table. No live/default changes.
# pathsolidbackstop=f: separate SOLID-table confirmation of an existing Hollow
# weak-backward seed-start candidate. Requires pathcandidates=decisions.tsv,
# hollow=f, queryspan=k, k>=67, pathseedlen>k and ALL other path policies off.
# Candidate input is the native13-column self/backward report (including
# self_count/external_count/stop_reason) from pathbackactualsupport=1. Rows must
# match queries by exact ID/order, seed length and k; raw original_length may
# exceed an already-guarded query length. Candidate decisions NEVER add counts.
# Only backward_dominance at exactly the seed boundary is eligible, with fixed
# context/competitor>=3, actual predecessor1 or2, and competitor>=2*actual.
# Truncate to the seed only if its FULL Solid terminal depth is<2. Depth>=2
# retains the saved suffix. Later candidates and all other stop reasons abstain.
# Report preserves Hollow backward counts and adds solid_seed_depth (-1 when
# not evaluated); stop_reason=solid_weak_backward_seed only when this veto fires.
# pathsolidsupport remains2 in this mode; thresholds are fixed for the experiment.
# Never edits seed bases or extends sequence. Not a blanket coverage requirement,
# live-extension change, default change, or proof that weak Solid seeds are unsafe.
# pathbackfork=f: modifier requiring pathsolidbackstop=t. Instead of immediate
# seed truncation, inspect at most k-64 saved choices after an eligible weak seed.
# Load a second explicit Solid64 table from the same original evidence. While
# full SolidK context depth<2, stop before the first short fork with context>=3,
# chosen>=3, other>=3 and 2*other>=chosen. Full-key depth>=2 ends this extra
# check permanently. No fork before recovery/horizon/end retains the whole path.
# The horizon is exclusive: end<seed_length+k-64. Seeds are never edited.
# Adds short_context/short_chosen/short_other/short_other_base/fork_end columns;
# -1/dot unless a fork fires. fork_end is the retained prefix length. Existing
# count columns remain the Hollow seed candidate; solid_seed_depth stays initial.
# Stop reason is solid_seed_short_fork. Requires additional heap for two tables.
# Fixed experimental thresholds; no live extension or default changes.
# pathshortstop=f: separate two-table diagnostic for weak Hollow choices opposed
# by stronger Solid64 support. Requires hollow=t, k>64, queryspan=k,
# pathseedlen>=k, queries= and outseq=; every other path policy must be off.
# Loads a second explicit Solid64 table from the SAME original in/extra evidence.
# Stops before a saved next base only if Hollow chosen<3, Solid64 chosen>=3,
# and strongest Solid64 alternative>=2*chosen. Missing/weak short support abstains.
# Never chooses the alternative or edits the seed. These fixed thresholds are
# experimental, not proof of an unsafe join. Two tables remain resident: budget
# additional heap. Report: ID, seed/original/retained lengths, k, short_k,
# hollow_chosen, solid_chosen, solid_other, other_base, stop_reason (none or
# weak_hollow_short_dominance). Evidence is -1/dot when no stop occurs.
# No change to Tadpole/BBMerge live extension or defaults.
# pathshortdrop=f: modifier requiring pathshortstop=t. Replaces its weak-choice
# predicate with a conjunction: Hollow context>=32, context>=4*saved-next count,
# Solid64 context>=2, saved-next count=0, and strongest alternative>=2.
# Neither Hollow collapse nor short contradiction suffices alone. All saved
# suffix positions are tested until the first conjunction; a failed early test
# does not end inspection. No seed edits or alternative-base selection.
# Appends hollow_context and solid_context to the short-path report; -1 when
# no stop fires. Stop reason is hollow_drop_solid_contradiction. Thresholds are
# fixed experimental choices, not a validated default. Original counting inputs
# and the two-table memory cost are unchanged; flag off preserves the old mode.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=8g --xms=256m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.KmerWindowDepth "$@"
