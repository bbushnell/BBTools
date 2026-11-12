#!/bin/bash
# Integrates completed CM searches and historical annotations for reference review.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden
Last modified November 12, 2026

Usage: dev/buildncrnareference.sh in=collector.tsv census=census.tsv memberships=members.tsv \
       scope=scope.tsv assets=models.tsv queries=queries.tsv out=new_prefix t=8

Builds evidence and parameterized tiers for an ncRNA evaluation reference.
Reuses completed whole-genome cmsearch results; performs no searches or calling.
Most modes are descriptive. Tiers applies the explicitly selected policy below;
successful serialization is not independent biological validation.
Annotation rows, raw model hits, overlapping intervals and genes are different units.

Required:
in=            Evidence mode: eight-field completed-search collector table.
               Profile mode: prefix of saved .summary.tsv, .coverage.tsv and
               .sources.tsv.gz files from one completed evidence run.
               Geometry mode: the same prefix, plus .hits.tsv.gz and .overlaps.tsv.gz.
               Attributes mode: the same prefix, plus .historical.tsv.gz and
               .relations.tsv.gz; also opens the original GFF paths for retained rows.
census=        Five-field path-based genome census.
memberships=   Ten-field source-membership table; original claims remain authoritative.
scope=         Six-field reviewed model scope, including reference_family/evidence_role.
assets=        Six-field captured-model metadata, including release, CLEN, GA and sha80.
queries=       Exact searched accession/name/CLEN roster.
out=           Fresh output prefix, or a fresh directory for tiers mode.

Options:
mode=evidence  Join completed searches and historical annotations.
mode=batchevidence Import completed per-model searches of whole-genome batches.
               Uses the complete batch collector and captured preparation maps;
               does not reopen original or prepared FASTA payloads. See schemas below.
searchspacemb= Required only for batchevidence: the original full search-space -Z
               in Mb, including both strands. Must match every native command;
               it may exceed the bases in an explicitly imported timing subset.
nativepaths=   Optional batchevidence-only TSV for recovered native directories:
               batch_id original_native_directory physical_native_directory
               It must cover every imported batch exactly once, with no extras.
               Original paths must match the collector; physical paths must be
               distinct existing absolute directories. Read recovered bytes from
               the physical location while validating command footers and emitting
               file/line witnesses against the ORIGINAL paths. No footer rewriting
               or line renumbering. Default: open the original directories.
               This map declares access locations; archive identity and verified
               recovery remain the producing workflow's responsibility.
mode=profile   Summarize saved tables only; does not open FASTA, GFF or CM payloads.
mode=scores    Describe bound batch scores, exact model spans and 80-90% GA floors.
               Requires completed batch evidence; legacy GA_pass is insufficient.
               Reads saved hits/overlaps and their population metadata only.
mode=retention Compare in=current_batch_prefix with compare=previous_legacy_prefix.
               Both require the same complete catalog and model metadata, saved
               hits/overlaps and profile inputs. Reopens old native reports to
               recover actual inclusion flags; does not rerun searches or read FASTA.
mode=tiers     Construct scored/masked per-input GFFs and complete event provenance.
               Requires complete bound batch evidence, attributes=, bound priors=,
               marginalpercent=, saved contigs/history/relations/batch_sources and
               the normal census/model metadata. Reads captured completion metadata
               to verify the original 75%-GA reporting floor; no alignment reports.
mode=tierreport Report an existing tiers directory without changing decisions.
               Uses in=reference_directory out=new_prefix; no census/model inputs.
               Separates raw-event causes from serialized physical mask reasons,
               reconciles both to saved totals, and counts measured/hitless priors.
samplemodel=   Optional accession for bounded original event and other-family
               overlap examples (for example RF03166); no accession by default.
examplespergroup=3 Maximum examples per tier/primary-reason/cross-family group,
               and qualifying other-family witnesses per example; range0..10.
marginalpercent= Explicit selected floor in [75,100], required only for tiers.
collapsepercent=80 Same-family grouping requires at least this percentage of BOTH
               intervals to overlap for EVERY pair of qualifying group members.
               This is a chosen geometry convention, not a measured optimum.
               Supports four decimal places within (0,100]; saved in rebuild config.
               Same-family support does not cause cross-family discordance.
               Broad tandem bridges and weak evidence remain witnesses only.
mode=geometry  Join saved model-hit events and overlap edges into review groups.
               Uses one thread; retains events, primitive line indexes and edges.
mode=attributes Join current original-GFF attributes onto all saved retained rows.
               Reads each distinct GFF path once, preserving every input/file/line key.
mode=strata    Join saved evidence and captured attributes into nonexclusive review strata.
attributes=    Completed capture prefix; required for strata or tiers.
attributesevidence= Optional original saved-evidence prefix for an earlier capture.
               Requires the same complete catalog and identical saved source claims,
               availability states and counts, plus every exact current source key.
               Reuses recorded claims; does not re-read live GFF payloads.
mode=prior     Bind unchanged prior decisions to saved census/contig observations.
priors=        Prior mode: six-column source manifest. Tiers mode: completed bound
               .prior_decisions.tsv.gz table (including an explicitly empty table).
mode=anchors   Join declared prior evidence to current saved hits; no admission.
mode=sourceanchors  Describe annotation-concordant hits from accepted source links.
sources=       Completed strata prefix; required only with mode=sourceanchors.
               Uses its .source_links.tsv.gz and .strata_summary.tsv alongside
               the original saved evidence and normal input/model metadata.
anchors=       Required six-field manifest: source_id format decisions bindings
               authority expected_rows. Formats: bound_prior_v1 (bindings=NA),
               tpp_construction_v1 (bindings=original eleven-field genome roster).
               Uses saved profile/hit/overlap/contig tables and normal model metadata.
               Accepted cores, accepted whole intervals and construction-eligible
               cores remain separate evidence bases. Literal contigs/lengths must
               match; outside input paths are retained without alias guessing.
mode=cost      Report a terminal native timing pilot; no searches.
               Uses only in=pilot_directory census=full_population_summary out=prefix.
               Legacy seven-column costs.tsv uses candidates.tsv, five phase times,
               completion records and final native file sizes. Expanded eleven-column
               costs.tsv uses model_assets.tsv, report.exit and saved native sizes;
               failed/short/unstarted rows remain unpriced. Both use population.tsv
               and start/end records. Expanded successful pilots use prepare/startup/
               searches phase times; total CPU is unknown after a failed pilot.
               Other reference metadata arguments are not used in this mode.
compare=       Retention: required previous legacy saved-evidence prefix.
               Cost: optional prior per-model cost table; require matching
               identities and numeric agreement at six-decimal reporting precision.
               Cost comparison is supported only for the legacy timing format.
mode=selectpanel Select original input IDs from a completed reference inventory.
               Requires census= memberships= reference= inventorymembers=
               selectioninventory=completed_inventory_prefix selectionpolicy=
               selectiontaxonomy= selectionfamily= out=fresh_prefix.
               Protects explicit reference-member support and retained controls,
               then fills by domain/phylum and least represented claimed TaxID.
               Taxonomy joins use exact declared source paths; unknowns remain
               explicit and eligible conflicting claims fail. No FASTA/GFF payload
               reads, sequence-equivalence collapse, calls or adjudication.
               Writes .selection.tsv (all original inputs and reasons),
               .inputs.config (input-only fragment in original-ID order),
               .reference.tsv (unchanged original-ID pointers), .summary.tsv,
               .support.tsv, .strata.tsv and final .provenance.tsv.
               Use per-input native evaluation for duplicate contig names across
               genomes; never manufacture historical source memberships.
               Reference-member coverage is not caller-subtype precision.
               The project recipe authenticates all inputs; this tool checks
               metadata consistency, not recorded digests against file bytes.
               Policy TSV header (one row per inventory family):
               family positive_target positive_ceiling protected_members
               retained_sources keep_protected_priors keep_multi_supported
               Lists are comma-separated or '.'; booleans are true/false.
               Mandatory positive IDs may exceed the target but never the ceiling.
               Retained masked-only/zero-target controls are additional inputs.
               Exact-total policy appends two columns to that header:
               total_inputs member_min_inputs
               Use positive_target=0, total_inputs>0 and a ceiling at least that
               total; protected_members='.' and both protection booleans=false.
               member_min_inputs is a declared ordered list, for example
               RF00174:300,RF01482:100,RF01689:100. Include every selected member.
               Count unique fixed controls C, including K positive controls and
               C0 nonpositive controls. Add total_inputs-C new positives; the
               native positive target is total_inputs-C0. Supported controls and
               overlapping member support receive credit without duplicate slots.
               Meet member minima in declared order with phylum/TaxID balancing,
               then fill the remaining positive budget. Insufficient support or
               an over-budget quota fails; the total is never silently expanded.
               Extended summary/support tables append control/budget accounting,
               quota order and nonadditive member input/scored-locus counts.
               A zero total with member_min_inputs='.' keeps ordinary policy
               behavior. Full original reference/mask pointers remain unchanged.
               Taxonomy manifest header: schema path source_root expected_sha80.
               gtdb_domain_phylum_basename_v1 reads headerless3-column maps and
               binds source_root/domain/basename.fna.gz exactly.
               x4_domain_phylum_basename_source_v1 reads headerless4-column maps
               and binds source_root/source/basename.fna.gz exactly; source is
               bacteria4 or archaea4. Supplied catalog aliases may bind these
               exact paths. Unmatched inputs remain UNKNOWN; no name guessing.
mode=inventory Inventory existing canonical reference GFFs by original input and
               explicitly mapped family/member. Requires census= memberships=
               reference= inventorymembers= witnesses= events= out=fresh_prefix.
               Member map header: family, member, reference_model (tab-separated).
               Keeps zero-target and masked-only inputs, original parent IDs and
               coordinates, qualified support, protected priors and source strata.
               Cobalamin recovery/exclusive counts share the native grader's rules.
               Member counts are nonadditive; positive inputs are not locus counts
               or proven independent organisms. No FASTA reads, panel selection,
               new adjudication or hash authentication. Authenticate inputs first.
               Writes .inputs.tsv, .members.tsv, .loci.tsv.gz, .summary.tsv,
               .strata.tsv and final .provenance.tsv; fresh outputs only.
mode=panels    Count source membership rows and distinct exact input IDs separately.
               Requires only census= memberships= out=prefix; no sequence reads.
               Writes .panel_population.tsv. Sources are declared panel labels;
               distinct input files are not automatically independent organisms.
mode=panel     Prepare one ordered evaluation panel from saved metadata only.
               Requires in=flat_caller_config census= memberships= sourcepanel=
               reference=full_reference_manifest observed=saved_contigs out=prefix.
               Every in=/in1= line selects one exact catalog member in sourcepanel.
               Nested configs, paired/comma-list inputs, repeated inputs and
               ambiguous selected contig names fail. No FASTA/GFF payload reads.
               Writes .genomes.tsv and .contigs.tsv for BindNcrnaCalls, a
               .selection.tsv crosswalk and finally .reference.tsv pointing to
               unchanged canonical GFFs. Unknown roster fields remain NA.
               No five-column truth or biological decisions are generated.
               This checks metadata consistency, NOT registry authorization:
               the project preparation recipe must verify the registered identity.
t=auto         Independent genome workers in evidence mode, source-file workers in
               attributes mode. Saved-table modes, including tiers, use one thread.
maxinputs=     Positive prefix length for a diagnostic sample; default all inputs.
               The ENTIRE collector/census join is validated even for a sample.
               Saved-table modes reject maxinputs; they summarize the saved population.
ow=f           Outputs must be fresh. Overwrite and append are unsupported.
pigz=f         Disables pigz. Other BBTools compression defaults remain active;
               use bgzip=f as well when explicitly requesting Java-only gzip.

Evidence outputs:
.contigs.tsv.gz       Current literal target names/lengths; not proof of past scan bytes.
.hits.tsv.gz          Every native model-hit event with model/range/truncation/edge data.
.sources.tsv.gz       Unique input/GFF claims, reused-membership counts and exclusions.
.historical.tsv.gz    Retained rows and explicit overlap-assessability states.
.relations.tsv.gz     Original Parent IDs/RFAM claims, dangling links and cycle witnesses.
.overlaps.tsv.gz      All historical/model and model/model overlaps. Opposite-strand
                     geometry is diagnostic, not support. Ancestor-only context has
                     no gene-overlap test. Components remain rows, not separate genes.
.coverage.tsv        Exact per-accession span/truncation event counts and denominators.
                     Zero-hit models have NA span/truncation and zero event counts.
.summary.tsv         Written only after successful workers and output closure.

Batch evidence also writes .batch_sources.tsv, declaring captured map/FASTA/member
paths and fingerprints once per batch. Retain that table and its referenced metadata.
Batch hit rows add native_target, native_included and GA_decision_basis; GA_pass uses
the native flag only after verifying --incT equals the original model GA. Displayed
bits are rounded, so they may disagree at the boundary. Neither native inclusion
nor a printed score establishes full length, a clean family or biological truth.
Batch history adds best_native_table: a physical line number alone is not a key.
Geometry pair/bridge reports retain each witness's own native table and physical line.
Summary status explicitly distinguishes the captured batch evidence format.
Profile, geometry, attributes and strata support that format. Prior/anchor/source-
anchor modes and the old GA-only draft constructor currently reject it; they must
not silently reinterpret native inclusion as the old printed-score comparison.

Cost outputs: .model_costs.tsv and .cost_summary.tsv. CPU scales by each sample's
nucleotide bases. Native-report size scales the same way; sample header bytes are
included. The legacy format's separate planning allowance gives each genome the entire multi-genome
sample report for each model; it is deliberately generous, NOT a statistical bound.
Bytes are logical file lengths, not allocated disk blocks. Neither size estimate
includes extra per-input metadata, copied models or temporary decoded genomes.
Queue wait, startup/I/O overhead and parallel efficiency are unmeasured.
Expanded-format projections include only completed points lasting at least ten
CPU seconds. If any model is unpriced, the summary labels the aggregate as a
measured-subset projection, never the whole-roster cost. Missing bytes are NA,
empty files are zero, and partial output status survives. Reused point CPU is not
charged to this pilot. Earlier sample points remain in the producer's points.tsv;
this report uses selected final points only. Expanded reporting consumes saved
telemetry without reopening native payloads and makes no historical-proof claim.

Anchor outputs: .anchor_rows.tsv.gz, .anchor_links.tsv.gz, .anchor_models.tsv,
.anchor_summary.tsv. Every prior row remains accounted for; exact repeated anchors
do not multiply observations. Missing, unsearched, multi-hit and weak-overlap cases
are explicit. Only a single own-model hit with reciprocal overlap >=0.5 contributes
observed score/coverage quantiles; full untruncated model margins are also separate.
Other-model overlaps remain evidence, never automatic subtype decisions. Quantiles
are min/p05/median/p95/max, with nearest-rank percentiles. They do NOT select cutoffs.
Thresholded searches cannot reveal the censored low-score tail. A missing hit never
revokes a prior decision, and construction eligibility is not scored reference truth.

Source-anchor outputs: .source_anchor_loci.tsv.gz, .source_anchor_models.tsv,
.source_anchor_summary.tsv. Each physical input/contig/strand/interval/model is
counted once, regardless of annotation-link count. Identical repeated native
observations give one score; conflicting observations give no chosen score.
Exact accession and family-only support remain distinct. CLEAR means no recorded
warnings on matching source-link rows, not independently reviewed truth. Warnings
and source kinds come from those same rows, never event-wide combined flags.
Source-kind counts are nonexclusive; locus and score-observation counts differ.
First exact/family witness lines refer to the saved source-link file, header included.
Literal tRNA/tmRNA types with NO_OVERLAPPING_SOURCE_RFAM supply a separate TYPE_FAMILY_ONLY
basis and type-witness line. Generic ncRNA, unknown RFAM and conflicting claims
cannot use this fallback. A type-only annotation never identifies a subtype.
No missing low-score observations are recovered; inherited search censoring remains.
No precision/FPR, subtype decision, admission cutoff or GFF is produced.

Profile outputs:
.model_profile.tsv      Per-accession event totals, native truncation counts,
                        exact full-model-span counts, p05/median/p95 and modal span.
.input_availability.tsv Each observed input's available/unavailable GFF claims.
.availability.tsv       Input counts for no available, available-only and mixed claims.
.profile_summary.tsv    Written last after all input checks and output closure.
                        Retains observed/catalog counts and full-versus-sample status.

Tier outputs (mode=tiers, out=new_directory):
reference.tsv               Completed manifest, published last; one GFF per input.
input_N.gff.gz               Scored/confident and dont_care/marginal physical loci.
event_provenance.tsv.gz      EVERY captured hit, including below-floor/uncertain hits.
annotation_witnesses.tsv.gz  Every primary annotation's model and prior witnesses.
prior_provenance.tsv.gz      All original prior decisions, without reclassification.
source_provenance.tsv.gz     Original source-row dispositions; source claims, not genes.
source_event_links.tsv.gz    Direct source/event links, including wrong-strand evidence.
event_counts.tsv            Per-model raw-event partition; never a gene count.
family_counts.tsv           Separate raw-event and physical-annotation counts.
source_family_counts.tsv    Source-row counts by declared family; unknown stays unassigned.
summary.tsv                 Per-input scored/masked annotations and protected priors.
parameters.tsv              Selected floor, score precision and explicit tier rules.
capture_thresholds.tsv      Original per-batch/model reporting and inclusion metadata.
sources.tsv                 Durable input dependencies for later re-thresholding.
rethreshold.config          Standard BBTools config retaining those exact dependencies.

Confident score eligibility uses native inclusion bound to the original model GA.
Full model span, no native truncation and clean evidence remain separate requirements.
Marginal masking requires the printed floor; a +/-0.05-bit interval touching that
floor stays explicitly uncertain and UNMASKED. Partial or discordant evidence below
the floor does not create a mask. Named accepted prior intervals/assignments remain
protected; new overlaps cannot widen them or excuse a known wrong-family call.
Exact input/target/family/strand/bounds duplicates collapse to one primary annotation
with every model measurement retained. Direct qualifying competitors use their own
intervals, never a transitive union. Original unsupported source rows stay provenance.

Each primary annotation has cm_models, cm_bits, cm_ga, cm_percent_ga, cm_score_basis
and cm_percent_precision. cm_models identifies measurement order independently of
reference_model (which preserves a prior's resolved assignment). Bits and GA remain
lossless saved decimals. Percent is rounded to six places, with the exact operand
pair retained; it never claims unrounded native precision. Hitless priors use '.'.
Complete CM attributes are required when present; no-score legacy GFF rows remain
supported unchanged. Events keep native inclusion and printed-GA disagreement flags.

To change the floor without another search:
  dev/buildncrnareference.sh config=old_reference/rethreshold.config \
      marginalpercent=85 out=new_reference
This recomputes tiers and qualifying competition, not just old labels. Keep the
saved evidence, its batch completion metadata, captured attributes/original binding,
bound priors and model/catalog metadata listed in sources.tsv. Original FASTA, GFF
and CM payloads and native alignment reports are not opened by this operation.

Tier-report outputs (mode=tierreport; completion requires .summary.tsv):
.events.tsv   Per-model/family tier counts with NON-ADDITIVE overlapping incidences.
.reasons.tsv  Exclusive primary reasons for RAW marginal events; bridges separate.
.physical.tsv Physical scored/masked counts, exclusive stored mask reasons, and
              protected-prior counts split into missing versus measured scores.
.examples.tsv Bounded original provenance rows for samplemodel, unchanged fields.
.example_overlaps.tsv Sample event ID plus qualifying same-strand other-family
              overlapping event; original intervals and measurements are retained.
.summary.tsv  Conserved populations and bounded-example counts, written last.
Zero-hit models remain in the event table; their family is '.' if no event binds it.
Counts are descriptive, not independent biological validation or new admission.

Retention outputs (mode=retention):
.retention_summary.tsv        Per-model and whole-population raw/native-included
                              retained, lost and added counts; written last.
.retention_differences.tsv.gz  Changed coordinate groups with duplicate counts and
                              representative old/new native-table witnesses.
Coordinates include input, literal target, model, strand, sequence bounds, model
bounds and native truncation. Paths/line numbers are witnesses, never join keys.
Duplicate multiplicity is preserved. Raw retention, native-included retention and
equal printed-score multiplicity are independent comparisons. Same-flag matching
takes priority when classifying additions; witnesses do not identify uniquely
paired duplicates or recover unrounded scores. Old saved GA_pass is not a native
flag. Old positive native tables must match saved rows and original --cut_ga
footers; zero-hit inputs remain covered by completed saved population metadata.
Historical model-file provenance remains the accepted generation's responsibility.
No biological admission or explanation of changed search output is inferred.

Score outputs (mode=scores):
.score_histogram.tsv    One-percent printed-score/GA bins, exact model span, native
                       truncation/inclusion and contig edge. Last bin is >=200% GA.
                       Direct other-family overlap counts are separate from totals.
.floor_curve.tsv        Per-model/inclusion counts at 80 through 90 percent GA.
                       Definite-pass/fail and rounding-sensitive counts partition
                       ALL events. Other columns describe printed-score passes.
.score_summary.tsv      Per-model denominators, native inclusion and disagreements
                       between that inclusion and printed score>=GA; written last.
Native scores have one decimal: +/-0.05 bits defines conservative rounding bounds;
equality at either bound remains uncertain. Native inclusion proves a pass at these
sub-GA floors, but a printed score alone does not recover an unrounded score.
Competitors are DIRECT valid same-strand other-family overlaps. Any reported
competitor, a printed-floor pass, a definite pass and uncertain qualification are
reported separately. A transitive review group is never a competing interval.
Exact model span and native truncation are not biological completeness decisions.
No floor is selected, no mask is applied, and event counts are not gene counts.

Span quantiles use weighted nearest ranks over raw hit events; tied modes choose
the smallest span. Zero-hit models retain explicit zero totals and NA quantiles.
Model span is measured in model positions, not genomic bases or gene completeness.
An available GFF was parsed, possibly empty; it is not proof of complete annotation.
Supply the SAME census, memberships, scope, assets and queries as the evidence run.
For example, reuse evidence.* with mode=profile in=evidence out=profiles plus those
five metadata arguments. No biological cutoff or gene admission follows a profile.

Geometry outputs:
.groups.tsv.gz         Connected same-input/target/strand review groups and counts.
.group_members.tsv.gz  Every original event with source key and model/range evidence.
                       Invalid intervals have NA group IDs, never inferred negatives.
.model_pairs.tsv.gz    Every saved model pair: equality, containment or partial overlap,
                       exact overlap/lengths and rounded reciprocal fractions.
                       Opposite-strand and unassessable edges remain diagnostic.
.bridges.tsv.gz        One deterministic disjoint full-model/untruncated witness pair
                       per bridging event and witness family. Retains both overlaps,
                       disjoint gap and whether the bridge covers the outer span.
.geometry_summary.tsv  Population/count conservation, original full/sample status.

Connected groups are NOT genes; chaining does not imply every pair overlaps.
Full model span does not prove biological completeness or independent genes.
Bridge witnesses do not prove fusion and never suppress their component hits.
All original edges/members survive; witness rows are not exhaustive triple lists.
No score winner, subtype assignment, overlap cutoff or biological span floor is used.
Native source/table/line pointers retain access to E-values and raw evidence.
The saved historical/model rows are counted but not biologically interpreted here.

Strata outputs (no genome/GFF/CM payload reads):
.event_strata.tsv.gz   Every native event, its group and separate model-span,
                      truncation, edge and source-origin counts/flags.
.source_strata.tsv.gz Every historical source key, row kind, selection and claim flags.
.source_links.tsv.gz  Every historical/model overlap with source kind and association.
.strata_counts.tsv    Nonexclusive event/group counts per queried accession and union.
.joint_strata.tsv     Joint coverage/truncation/GA-margin/source/alternative/edge bins.
.strata_summary.tsv   Written last; retains observed/catalog and original scope.

RFAM claims use the complete captured model roster, including unqueried models.
Agreement, same-family alternatives, other families, related classes and unknown
accessions remain distinct. These are source claims, not independent confirmations.
Qualified RFAM comparisons include regulatory and tagged-gene rows. Source links
retain their literal original feature type and spatial/strand association.
Pseudo/biotype/partial/exception event flags use directly overlapping RNA-type rows
only. Component, gene/transcript, other and opposite-strand claims remain in source
rows/qualified links; they are not inherited as RNA claims. Parent relationships
remain in the original evidence; this mode does not infer inherited attributes.
Repeated contradictory booleans survive. Free-text notes do not become claims.
Groups may occur under several accessions: accession group counts are NON-ADDITIVE.
The joint table counts each event once per accession and once in ALL_QUERIED.
Each group is counted once within a displayed bin but can occupy several bins.
Coverage and score-margin bins are descriptive intervals, NOT acceptance cutoffs.
Exact bits, model GA and margin remain in the event table. Availability, RNA-row
overlap and RFAM agreement are separate; an untagged RNA row is not an absent GFF.
Same-family CM alternatives use the supplied reference_family, not formal clan IDs.
Only queried-model alternatives are observed; unsearched models are not excluded.
No span cutoff, biological admission, junk exclusion or prior-decision import occurs.

Attribute outputs:
.attributes.tsv.gz     Original input/GFF/line key plus the current raw ninth field.
.attribute_sources.tsv Counts per source path, including shared input memberships.
.attribute_summary.tsv Population and matched-row counts; written after worker closure.

Attribute mode matches saved target/type/coordinates/strand/ID and complete decoded
Parent/RFAM claim lists. Missing files/lines or mismatches fail without a completion
summary. Original evidence is unchanged. Comments and blanks count toward source
line numbers; embedded FASTA is not consumed as annotation. No FASTA file is opened.
Raw attribute bytes receive the normal TSV percent/tab/CR/LF escaping: decode that
outer framing before interpreting the original GFF escapes. Other attributes are
captured NOW; matching saved identities does not prove uncaptured historical bytes
were unchanged. Pseudogene/biotype/partial claims are preserved, not adjudicated.
The inherited saved ID token NA conflates absent ID with a literal ID named NA.
Matching that rendering cannot establish historical ID presence; the current raw
attribute field retains current presence. Input/GFF/line is the identifying key.

Prior-decision mode reads the saved profile inputs plus .contigs.tsv.gz, and the
tables named by priors=. It never reads FASTA/GFF/CM payloads or requires a new hit.
Manifest header: source_id format decisions bindings authority expected_rows
The decisions/bindings paths are absolute. Supported versioned formats:
  cobalamin_panel_v1: decisions is the ten-field CobalaminPanelPrep locus table;
    bindings is its eleven-field input roster, joining original genome to FASTA.
    These are the preserved outputs of dev/cobalamin_panel_prep.sh, not new
    decisions made by this importer. Exact headers (spaces below show columns):
    decisions: locus_id genome seqid source_family member_family source_start1 source_stop1 core_start1 core_stop1 strand
    bindings: index genome assembly taxid domain phylum gff fasta loci genome_bases contigs
  fmn_adjudication_v1: decisions is the original sixteen-field FMN decision table;
    bindings has header: decision_line fasta seqid start1 stop1 strand classification
    Every physical decision line (header is line 1) must match all binding fields.
    Decisions are reviewed input, with this exact header:
    seqid start1 stop1 strand classification winning_model verifier_identity model_length original_annotated_members best_CM_bits CM_GA_bits CM_status decision_source independent_reviewer trace_packet_sha80 CM_packet_sha80
    Supply the reviewed table and its explicit binding table; this tool does
    not invent a reviewer, classification or generator for those decisions.
    Accepted classification tokens are ACCEPTED_NOVEL_FMN and UNRESOLVED.
Only exact census or explicitly recorded supplied-FASTA paths may bind an input.
Unknown/ambiguous paths, missing literal contigs, invalid intervals and changed
row counts fail. Original classifications and complete source rows survive.
Outputs: .prior_decisions.tsv.gz, .prior_counts.tsv, .prior_summary.tsv (last).
FMN intervals remain prior-reviewed whole intervals, not newly adjudicated aptamer
bounds. All cobalamin states survive; zero core bounds retain absence semantics.
New-event association and new biological admission are separate from this import.

For example, import preserved decisions using the SAME legacy evidence metadata:
  dev/buildncrnareference.sh mode=prior in=evidence census=census.tsv \
      memberships=members.tsv scope=scope.tsv assets=models.tsv queries=queries.tsv \
      priors=prior_sources.tsv out=bound_priors
prior_sources.tsv has the six-column manifest header above, one row per source,
absolute decisions/bindings paths, a declared authority and exact data-row count.
The evidence prefix must already have completed profile outputs and saved contigs.

For anchors, bound_prior_v1 consumes the resulting bound_priors.prior_decisions.tsv.gz
as decisions, with bindings=NA. It does not consume the original decision table.
tpp_construction_v1 instead consumes the preserved construction-review table and
the same eleven-field input roster shown above. Its exact decision header is:
locus_id genome assembly taxid domain phylum seqid strand intervals1 length source_line source_ids partial pseudo status contig_length edge undefined_bases lflank rflank source_lines cm_state raw_hits overlapping_sense_hits cm_model_from cm_model_to cm_seq_from cm_seq_to cm_trunc cm_bits overlap_bp left_excess right_excess geometry selected_phase selected_state selected_context_start selected_context_stop selected_bits selected_genomic_start selected_genomic_stop construction_candidate disposition coverage_critical_note initial_construction_candidate initial_disposition
This is a preserved construction input; mode=anchors does not create it.
The anchors manifest uses source_id format decisions bindings authority expected_rows,
with absolute file paths except the literal NA binding for bound_prior_v1:
  dev/buildncrnareference.sh mode=anchors in=evidence census=census.tsv \
      memberships=members.tsv scope=scope.tsv assets=models.tsv queries=queries.tsv \
      anchors=anchor_sources.tsv out=anchors
Keep the same saved evidence/profile/contig inputs and declared source authority.

All source/native line numbers are one-based. Coordinates are one-based inclusive.
The best marker orders same-strand, in-range evidence by overlap, model coverage,
bits, E-value and lexical ties; it neither resolves biology nor discards alternatives.
The original GFF and native table are referenced by path/line and must be retained.
Selection bits are RNA=1, RFAM=2, unknown=4, descendant=8, ancestor=16.
Native truncation is 0=none, 1=5-prime, 2=3-prime, 3=both, not a gene-fragment label.
Text fields escape percent/tab/CR/LF as %HH and otherwise preserve UTF-8.
Record order may differ with worker count. Compare keys/rows, not compressed bytes.
No compartment is inferred from a model name or a source panel's label.

Input headers (literal tab-separated fields; the following spaces show columns):
in: input_id fasta result_directory models raw_model_hits search_CPU_seconds decode_CPU_seconds max_native_RSS_KiB
census: input_id fasta file_bytes nucleotide_bases fasta_records
memberships: source member_index input_id supplied_fasta source_genome source_assembly taxid declared_compartment reported_gff reported_sha80
scope: CM_accession biological_label search_role scoring_scope reference_family evidence_role
assets: accession release CLEN GA_bits sha80 path
queries: accession name CLEN

scope evidence_role must be exactly TARGET_CLASS or RELATED_CLASS. These are
reviewed evidence roles, not biological admission decisions. Each accession has
one scope row; biological_label, search_role and scoring_scope are nonempty,
and reference_family uses only letters, digits, underscore or period. All assets
must have one matching scope row and the same declared release. A minimal scope
row for a target model (tabs between fields) is:
RF00174	cobalamin	target	genome	cobalamin	TARGET_CLASS

Each result_directory contains hits.tbl (native Infernal 1.1.5 cmsearch table),
input.tsv and completion.tsv. input.tsv must bind the original search to the
same input_id, exact FASTA path, nucleotide_bases and FASTA_records as the census;
matching hit counts or dimensions alone do not identify a search. Its header is:
input_id fasta nucleotide_bases FASTA_records
The completion.tsv header is:
input_id model FASTA_records input_bases searched_residues reported_hits status
There must be one COMPLETE row per queried model, including zero-hit models.
searched_residues is twice input_bases for the required whole-genome two-strand scan.
Source indices in memberships are one-based per source. Unknown claims use NA.
reported_sha80 remains a claim, not a new actual-file identity check.

Batch input schemas (literal TSV fields; spaces below separate column names):
collector: batch_id prepared_directory result_directory models raw_model_hits
Each prepared directory contains the batch writer's members.tsv, targets.tsv and
summary.tsv, plus a completed one-row provenance binding.tsv:
batch_id fasta fasta_sha80 mapping mapping_sha80 members members_sha80 genomes nucleotide_bases fasta_records status
status is PREPARATION_BOUND. Paths name that directory's batch.fa, targets.tsv and
members.tsv. Its prepared summary must conserve every original input, contig and
base. The search directory's input.tsv retains exactly that accepted binding row.
Fingerprints identify the recorded preparation; this reader does not recompute
payload hashes or certify arbitrary altered metadata from matching totals alone.

The batch search directory requires exit.txt=0 and completion.tsv:
batch_id model FASTA_records input_bases searched_residues reported_hits original_GA_bits report_bits status
Every queried model has exactly one COMPLETE row, including zero-hit models.
Each model's subdirectory (named by Rfam accession) contains exit.txt=0, hits.tbl
and search.out. Actual query/name/CLEN, target dimensions, hit counts, native
program/version/completion and original-GA command binding are checked. All searches
must use the declared full -Z, --cpu 1, --nohmmonly, --noali and default filters.
Unknown synthetic targets, foreign/duplicate/incomplete mappings and invalid bounds
fail. Real native paths, physical lines, model, strand, truncation and raw records
are retained; no per-genome native report or CPU receipt is fabricated.
Per-input/model hit counts are derived only after every batch/model output succeeds.
Missing or failed output never becomes a zero-hit pair. Repeated original names in
different inputs remain distinct; repeated names within one input are rejected.
HELP
}
if [ $# -eq 0 ] || [ "$1" = -h ] || [ "$1" = --help ]; then
	usage
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ]; do
	cd "$(dirname "$SCRIPT")"
	RESOLVED_SCRIPT="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$RESOLVED_SCRIPT"
done
DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)/"
popd > /dev/null
if [ -f "${DIR}bbtools.jar" ]; then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=4g --xms=256m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.BuildNcrnaReference "$@"
