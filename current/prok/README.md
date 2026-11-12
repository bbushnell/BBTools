# Prok Package - Prokaryotic Genome Annotation and Gene Prediction

BBTools prokaryotic package providing comprehensive tools for bacterial and archaeal genome annotation, gene prediction, ribosomal sequence processing, and statistical modeling of genomic features.

## ncRNA resource development

`dev/summarizencrnabenchmark.sh in=times.tsv gradeinputs=grade_inputs.tsv grades=grade_summary.tsv raw=grade_raw.tsv policies=shipped,selected workers=1,16 repeats=3 out=report`
summarizes completed ordinary-caller benchmarks with a declared complete roster.
Native grade tables carry three additional leading fields: policy, threads and
repeat. It validates input/reference bindings, scored and masked parent universes,
raw query/reference counts and endpoints against the native summaries. Exact
input/family/reference keys identify lost/gained detections, including repeat and
worker-count differences; equal TP totals do not imply equal detected loci.
Optional `taxonomy=` binds literal original input IDs/FASTA paths to declared
domain, phylum, category and taxonomy-status labels. Additional class reports
retain zero-query inputs, exact first-detection/last-loss input IDs and unknown
labels; their counts reconcile with the existing family totals. This metadata
does not change scoring or infer taxonomy, rarity or biological independence.
Paired optional `subtyperaw=` and `subtypegrades=` inputs add exact changes in
the existing cobalamin subtype views. Both use the same policy/threads/repeat
prefix as the ordinary reports. Native subtype totals are reconciled with their
reference and query rows; exact original query lines, models and coordinates
bind each witness to the ordinary union. Separate reports distinguish a lost
subtype match that remains union-detected through another model, a same-subtype
rematch, and an ordinary union miss. Gains are also reported, rather than inferred
from net TP differences. Both sides' witnesses remain visible; views stay
nonadditive and do not establish biological subtype identity. Existing ordinary
outputs are unchanged; no caller, regrading or new reference admission occurs.
CPU medians use user+system per replicate. Optional `counts=` preserves original
CDS/tRNA/rRNA totals from the binder; projected generic-RNA counts are separate.
Optional `profile=` plus `profiletime=` reports disjoint measured worker CPU over
that profile invocation's process CPU, with shared CDS/RNA path excluded from
RNA-only totals. Unmeasured loading and other-thread costs remain unallocated.
The tool reads no sequences and performs no calling, alignment, regrading or
reference admission. It does not authenticate producer declarations or imply
independent biological genes. See the launcher for exact schemas and outputs.

Explicit `timingonly=t` instead accepts the complete timing roster and optional
profile pair, rejects grade/count inputs, and emits no grade or comparison tables.
Its timing table explicitly states that caller outputs were not compared.
It reports that accuracy and detection stability were not evaluated. This permits
timing a larger workload without repeating a separately completed accuracy study.

`dev/selectncrnabenchmarkinputs.sh in=genomes.tsv memberships=memberships.tsv taxonomy=taxonomy.tsv sizes=1000,11769 out=panel`
selects whole-input benchmark metadata without consulting reference/call counts or
opening sequences. Domain/phylum groups each receive one input, with remaining
slots apportioned by population using largest remainders. Within each group,
systematic sampling over joint log2 genome-size/fragmentation bins uses seeded
ranks within bins. Explicit unknown taxonomy is retained. The minimum-per-group
rule overrepresents small groups; inspect the joint composition report against
the full supplied census. The full-population size preserves every original path;
all outputs retain census order and literal UTF-8 paths. This represents the
declared corpus, not an unbiased census of prokaryotes or independent assemblies.
The owning recipe records corpus choice, source identities, seed and sample rung.
The `.totals.tsv` table gives corpus and selected input, nucleotide-base and
FASTA-record totals from the saved metadata, without opening sequence payloads.

`dev/partitionncrnabenchmarkinputs.sh in=panel.selection.tsv census=genomes.tsv memberships=memberships.tsv contigs=evidence.contigs.tsv.gz out=new_directory`
keeps each selected original input intact and assigns it to the first batch with
no repeated literal contig name. It reads saved metadata only, preserving exact
paths, decoded names, lengths and census order. Duplicate names within selected
inputs fail; path/base/record dimensions are checked for the complete catalog.
Each batch has flat caller inputs, binder-compatible
genome/contig rosters and an explicit binder file budget; zero-call inputs remain
represented. The summary reports achieved greedy batches and the maximum per-name
input multiplicity lower bound without claiming an optimal partition. Reuse the
same recorded partition across caller builds. Calling, binding, grading and
authentication remain the owning workflow's responsibility. See its help for
schemas and generated-label semantics; no FASTA headers are changed.

`dev/ncrnachangeevidence.sh in=sweep.changes.tsv inputs=sweep.inputs.tsv policies=policies.tsv calls=sweep family=cobalamin_ out=evidence.tsv`
extracts original per-input/policy stage and GFF rows around every changed reference
locus. It binds sidecar ordinals through the native input report and explicit
joint-policy roster, reads each selected file once, converts oriented inclusive
coordinates to genomic coordinates, and retains opposite-strand features. Raw
rows are percent-framed with their source paths and line numbers. `padding=0`
adds no inspection context; `models=t` includes individual model observations.
This is a read-only evidence view, not a caller, grader or causal classifier.
Verifier acceptance precedes snapshot/path selection; missing observations of
reused work are not evidence of rejection. No FASTA payloads are opened.

`dev/selectncrnascenarios.sh in=global_ga.tsv.gz in2=global_margin.tsv.gz calls=calls.tsv bindings=bindings.tsv models=models.tsv scope=scope.tsv callerconfig=original.config out=new_directory`
selects the same exact inputs and enabled/supported families from paired global
scenario tables. The original flat configuration must explicitly declare each
required canonical gate and euk=f; unsupported aliases/only-modes are rejected.
The complete input list must match the saved-call adapter bindings, including
original-to-canonical FASTA aliases and zero-call inputs. Every generic query is
attributed through the declared model-prefix/display-label map; unknown models,
conflicting labels and disabled-family calls fail. Stable input-local decision
keys must be unique and identically ordered across scenarios. Included row bytes
and original witness pointers survive unchanged; exclusions are counted by family,
reason and before/after state. The final summary is written after all other outputs
close. Original global tables and their source companions remain required. This
tool selects diagnostic scope; it does not select biological thresholds or truth.

`dev/draftncrnascenarios.sh in=evidence_prefix census=genomes.tsv memberships=members.tsv scope=scope.tsv assets=assets.tsv queries=models.tsv policy=draft.tsv priors=bound_prior.tsv.gz links=source_links.tsv.gz out=new_directory`
constructs two hypothetical disposition tables from saved GA-filtered evidence and an
explicit draft policy: structural checks at each model's GA, and the same route
with its declared additional margin. Saved GA flags must match the captured model
thresholds, and any below-GA native event is rejected. Prior accepted positives are preserved;
other historical classifications gain no inherited masking exemption. Every
queried model participates in competition checks. Exact native duplicates collapse,
conflicting duplicate observations stay unresolved, and same-family masks merge
only at identical physical bounds. All original events and prior rows retain
witnesses; the original source-link table is a required companion. Counts separate
raw events, unique model intervals and emitted mask rows. The tool does not choose
an operating point, publish a reference GFF or certify the inferred family as
biological truth. Global tables still require exact per-baseline input/family
selection before use with the saved-call comparator. Run without arguments for
the supported draft actions and schemas.

`dev/bindncrnacalls.sh in=combined.gff.gz genomes=panel.tsv panelcontigs=panel_contigs.tsv observed=contigs.tsv.gz census=genomes.tsv memberships=memberships.tsv callerconfig=original.config baseline=name inputsha80=20hex out=new_directory`
binds a saved combined CallGenes GFF to explicit original panel metadata. It keeps
every generic RNA row, including duplicates and unknown/excluded models, in
per-input query GFFs. Literal contig names must uniquely identify panel inputs;
the original caller configuration and saved contig observations must agree with
the complete panel. Empty RNA inputs remain explicit. A physical-line map links
projected rows back to the original GFF. Other feature types are validated and
counted but remain only in the original file; projected zero counts are not
evidence of their absence. Retained arrays, logical output bytes and file counts
have configurable limits. The final query manifest is published only after all
outputs close. The supplied fingerprint is recorded as a declaration, not checked
by this tool. No FASTA payloads, calling or biological reference decisions occur.

`dev/comparencrnascenarios.sh in=before.tsv in2=after.tsv calls=calls.tsv contigs=contigs.tsv familymap=families.tsv out=new_directory`
compares two explicitly supplied hypothetical dispositions against the same saved
calls. An input-local decision identifier declares correspondence across coordinate
or family changes; missing keys remain additions/removals. It retains every row,
reports model/family counts, transitions and hypothetical call outcomes, and uses
the existing scored-first matcher. Run without arguments for schemas and outputs.
It evaluates no biological predicates, requires no invented approval claims, and
publishes no reference GFFs or reference manifest. Metric changes describe the
changed assumed references, not measured improvements in biological accuracy.

`dev/writencrnareference.sh in=dispositions.tsv inputs=inputs.tsv contigs=contigs.tsv rules=rules.tsv out=new_directory`
serializes explicit dispositions into per-input reference GFFs, including empty
views. It preserves all original decision rows and rule-authority declarations in
companions. Only resolved scored/masked rows with valid observed contig bounds
are emitted; duplicate or contradictory physical dispositions fail. The completed
grader-compatible manifest is published last. Run without arguments for schemas.
The tool requires declared rule approval but cannot independently verify it. It
does not select biological thresholds, infer families or decide which loci to admit.

`dev/gradencrnareference.sh in=queries.tsv reference=references.tsv familymap=families.tsv out=prefix`
grades per-input CallGenes output against explicit scored/don't-care reference
GFFs. Exact input populations and original FASTA paths must agree, including empty
query files. It retains duplicate/excluded calls, reports family TP/FP/FN and
endpoint offsets, and gives scored matches precedence over same-family masks.
Run without arguments for manifest schemas and report columns. Manifest agreement
does not independently verify producer provenance; combined saved GFFs need their
original input mapping. This tool neither admits biological loci nor estimates
genomic false-positive rates; masks encode uncertainty and can hide false calls.

Optional `subtypes=map.tsv witnesses=annotation_witnesses.tsv.gz events=event_provenance.tsv.gz`
adds cobalamin subtype reports without altering the ordinary union reports.
The map names each subtype, caller prefix and reference model. Protected prior
assignments remain authoritative; only qualifying confident group members supply
new scored subtype support. Marginal support remains unscored and weak/bridge
observations confer none. One invocation computes model-support recovery and
exclusive subtype accuracy independently from the same loaded loci. The exclusive
view masks multiply supported parents in their supported subtypes before matching,
allowing reassignment to a single-supported reference; unsupported third subtypes
gain no ambiguity exemption. Duplicate accounting remains independent of masking.
View-tagged raw/summary reports, retained parent support and distinct ambiguity
counts support later presentation choices without new grading. Recovery misses
are not missed physical genes or proof of wrong subtype. These are reference-
conditional, nonadditive results. The project authenticates both consumed sidecars
against its registered manifest; the native reader checks their internal bindings.

`dev/buildncrnareference.sh mode=sourceanchors sources=strata_prefix` measures
annotation-concordant score/coverage distributions from saved qualified links.
It counts physical intervals once, separates exact-accession/family-only support
and recorded source warnings, and retains source-row witness pointers. Repeated
links do not inflate observations; conflicting native scores are not arbitrated.
These are source-supported observations, not independently reviewed truth or FPR.
Literal tRNA/tmRNA types with no RFAM claim form a separately labeled family-only
basis; unknown/conflicting RFAM claims are never bypassed through that fallback.

`dev/buildncrnareference.sh mode=anchors` joins declared prior decision tables to saved
model evidence, using the normal input/model metadata plus `anchors=manifest.tsv`.
Accepted cores, accepted whole intervals and construction-eligible cores remain
distinct. Missing, multi-hit, weak-overlap and outside-path cases are reported;
single-hit score/coverage quantiles are observations, never admission thresholds.
The launcher documents versioned inputs, literal-contig checks and censoring limits.

`dev/buildncrnareference.sh mode=panels census=genomes.tsv memberships=members.tsv out=prefix`
counts membership rows separately from distinct catalog input IDs and bases for
each source panel. Repeated memberships do not inflate input counts. This uses
saved metadata only; source labels are not inferred biological families and
distinct input files do not establish organism independence.

`dev/buildncrnareference.sh mode=cost in=pilot_directory census=population_summary.tsv out=prefix`
reports completed native search timings and logical native-report sizes. CPU and
byte projections scale each model by its own sample bases. The legacy format's generous
per-genome report allowance is a planning scenario, not a statistical bound.
Expanded timing tables also accept stopped pilots with completed telemetry:
failed, short and unstarted points remain unpriced, partial byte counts remain
partial, and incomplete rosters receive only a measured-subset projection.
This format reads saved telemetry rather than reopening native payloads; legacy
seven-column reports retain their original behavior and comparison support.
Temporary genome decoding, metadata files, queue time and parallel efficiency are
not included. This mode makes no searches or biological reference decisions.

`dev/buildncrnareference.sh` joins completed whole-genome CM searches with a path-based
genome census and original annotation memberships. It writes descriptive hit,
source, overlap, contig and exact model-span tables; this implementation does not
admit biological loci or produce an accepted truth GFF. Missing annotations,
unknown strands and unresolved coordinates remain explicit. Original GFF/native
files and the supplied model-scope tables must be retained. Run without arguments
for the schemas, required paths and output semantics. The finite native regression
suite is `dev/ncrna_reference_builder_test.sh out=new_fixture_directory`.
`mode=profile in=evidence_prefix` reads the saved summary, coverage histogram and
source-claim table with the original five metadata arguments. It reports weighted
model-span quantiles and per-input annotation availability without opening genome,
annotation or model payloads. Zero-hit models remain explicit. A parsed empty GFF
counts as available, never as a biological negative; full model span is not proof
of a complete gene. The profile mode chooses no cutoffs and admits no loci.
`mode=geometry` uses that saved evidence prefix plus its hit and overlap tables.
It preserves every event while reporting same-strand connected review groups,
exact/contained/partial overlaps, opposite-strand diagnostics and possible bridges
across disjoint full-model hits. Groups are not gene counts, and no hit is removed
or chosen as a winner. Original source keys, range, truncation and model roles
remain available for subsequent biological review.
`mode=attributes` joins every retained historical row to its original GFF path and
line, matching the saved identity and Parent/RFAM claims before capturing the raw
attribute field. Shared paths are read once and dispatched to each original input
key. It reads source annotations without rereading genomes. Newly captured flags
are current source claims; agreement with saved identity fields does not certify
uncaptured historical attributes or admit/exclude genes.
The saved ID rendering `NA` cannot distinguish an absent ID from a literal ID
named NA. Matching it does not verify historical ID presence; the current raw
attribute field preserves current presence and the input/file/line key stays unique.
`mode=strata in=evidence_prefix attributes=capture_prefix` describes nonexclusive
model/span/truncation/edge/source-claim strata using saved evidence only. Every
native event, source row and historical/model link remains keyed. It separates
RNA-row claims from component, gene/transcript context and other row kinds;
context pseudo flags and opposite-strand claims never become RNA-event flags.
Parent-associated claims remain available through the original relationships,
without inferred inheritance. RFAM interpretation uses the full captured roster,
including unqueried models, and retains unknown/conflicting claims. Qualified RFAM
comparisons remain visible for regulatory or tagged-gene source rows; source links
preserve the literal original feature type and association. Group counts are
non-additive across accessions, not gene totals. Prior decisions and biological
admission rules are separate future inputs, not implemented by this mode.
The joint-strata table cross-tabulates model coverage, native truncation, margin
above the model's gathering threshold, source association/RFAM claims, alternative
hits and physical edges. Its display bins are not acceptance thresholds. Events
are counted once per accession/union; a group is counted once per bin and may
occupy multiple bins. Exact scores and margins remain in the event table. Named
same-family alternatives use the supplied family map, not formal clan membership;
models that were not searched cannot be ruled out by this report.

`mode=prior in=evidence_prefix priors=manifest.tsv` binds existing cobalamin-panel
and FMN decisions to saved input and literal-contig observations. It preserves
every original row, assignment, authority and source/core interval, including
accepted loci absent from the new thresholded scan. Original FASTA paths must
match the census or an explicitly recorded source membership; missing or conflicting
identities fail. FMN bounds remain prior-reviewed whole intervals, not new aptamer
bounds. This mode reads metadata only and neither requires hits nor admits new loci.
The launcher documents the versioned source schemas and three output tables.

`dev/samplegenomicbackground.sh` uses `CobalaminShredSampler` to select whole source
groups and sample annotation-masked genomic windows. Optional lineage/GC quotas
and nonoverlapping cores preserve a fixed empirical-length schedule; default
cobalamin behavior remains available. These backgrounds are provisional until
independently screened. See `docs/guides/SampleGenomicBackgroundGuide.txt` for
input schemas, population semantics, outputs and failure handling.

`reducencrnamodels.sh` selects a compact consensus/HBM library using exhaustive
gene-to-model support, optionally protects whole-genome detections, and can
regenerate retained models with paired rollback. `reduce=f regenerate=t` rebuilds
the complete supplied library without reduction, with caller verification on by
default. Explicit `callerverify=f` produces a labeled, unverified fixed-library
proposal using corpus retention alone. `graph=` reuses saved ORIGINAL or REDUCED
support instead of repeating the initial alignments. `mode=restore callerverify=f`
restores paired models from supplied GFF evidence without running a caller and
labels the proposal as not caller-verified. The full user guide is
`docs/guides/ReduceNcrnaModelsGuide.txt` in the BBTools release.

`mode=scan` performs a pure exhaustive Quantum window assay. It scores both
orientations by default in that mode (`rcomp=f` is an explicit diagnostic), retains
the ORIGINAL scored graph with a window/strand mapping, and summarizes requested
`reportcuts=` without repeating alignments. Acceptance is counted once per physical
window. No model reduction, regeneration, caller or library export runs; an optional
HBM validates pairing only. A background FPR claim still requires independently
screened negative inputs. Ordinary reducer modes retain forward-only behavior.

`buildncrnaboundarynets.sh` learns endpoint features, compares native classifier
and regression networks through ordinary CallGenes, and exports the accepted
family resources. Its guide is `docs/guides/BuildNcrnaBoundaryNetsGuide.txt`.

`dev/ncrnaendpointsummary.sh` summarizes existing nine-column matched-locus reports
without calling or matching again. It reports signed/absolute endpoint errors,
nearest-rank quantiles and joint both-end reachability for chosen radii. Supply
the report's exact family label and overlap rule. These are coordinate diagnostics,
not predicted neural-network recovery or caller-window containment.

`dev/ncrnagccomparison.sh` measures numerical differences through the actual native
boundary vector writer, classifier/regression readers and Scavenger GC input.
It selects at most 60 real contigs from an existing admitted training population,
retaining empty bins and missing-model counts. `requireexact=t` checks every
written/read feature against the raw values and scorer with matching leave-one-out
counts; it does not equate those counts with ordinary full-count inference.
No network is trained. The boundary writer uses eighteen fixed decimal places
to preserve these features through native decimal readers; this is not an
arbitrary-float encoding API or a measured accuracy improvement. Run the launcher
without arguments for input schemas, sampling bins and output descriptions.

`dev/ncrnaendpointscreen.sh` screens endpoint word lengths and inside/outside
compositions on the complete original population. Tables learn only the native
internal training partition; training uses key-level leave-one-out and validation
uses unchanged training counts. It reports smoothed count-decisiveness curves
across valid nearby endpoints, with explicit ties, abstentions and missing/ambiguous
denominators. It performs no HBM alignment, source-genome reread or training.
These curves screen candidate geometries; final choices still require native
classifier/regression endpoint results and whole-genome caller evaluation.

Both guides contain complete short commands using
`testdata/ncrna_model_tools/`, a small synthetic example derived from the shipped
6S library. Run them from the BBTools installation directory containing the shell
launchers and resources. The examples demonstrate program behavior; they are not
biological accuracy estimates or recommended production training recipes.

Experimental endpoint implementations can install an
`NcrnaBoundaryRefiner.Factory` on a loaded `NcrnaFamily` before workers start.
The factory creates independent worker state. The caller checks for a refiner
after alignment and model selection, before path selection; no check is added
inside the alignment loops. A null factory preserves ordinary endpoint handling.
Use `CallGenes.run(timer, setup)` to configure the family inside the caller's
resource-cleanup scope. Production code depends only on the generic interface.
The preserved mitochondrial experiment and its driver live in `test/prok` and
are invoked explicitly through `dev/mito_endpoint_callgenes.sh`; they are not
enabled by ordinary CallGenes. Their class headers describe why they were kept.

## Optional TPP and glycine riboswitch calling

`callgenes.sh tpp=t` enables the packaged 425-model TPP riboswitch family.
`tppcdsstop=t` separately limits additional three-prime endpoint extension across
a downstream same-strand CDS start, preserving existing overlap. Both options
are off by default. With both enabled on 3,210 genomes, the final family detected
5,763 confirmed references, missed 101 and produced 500 unmatched scored calls;
355 calls overlapping uncertain reference regions remained unscored. Both ends
were exact for 2,362 matched calls. Compared with the initial 1,442-model library,
two reference detections were lost and five gained: one lost locus stays shortened
and bypasses refinement, and one is removed by gene-path selection. These are
in-corpus measurements, and unmatched calls are not all established false positives.
Adaptive shortlisting preserves all 5,763 detections while reducing reported
model-verification alignments from 3,359,821 to 1,725,570 (49%). Both-exact remains
2,362; one retained locus changes from offsets (0, -6) to (-1, -1) bases at its
five-prime and three-prime ends. The selected mapping policy uses seven-base
words, a minimum of five shared words, at least 0.1875 of the best model's count,
and at most 361 models. The endpoint network was fitted against the 425-model
library before this shortlist selection and evaluated with the selected policy.

`callgenes.sh glycine=t` enables the packaged 329-model glycine riboswitch
family. Its models represent individual aptamers, including separately supported
adjacent aptamers; there is no composite tandem model. Calling is off by default.
Within the enabled family, the packaged nine-base endpoint network refines the
reported coordinates. Original alignment bounds govern acceptance, neighboring
candidate discovery and final-path selection. Applying refinement only after
selection prevents endpoint movement from suppressing neighboring detections.
Reported intervals can consequently overlap. Glycine WINDOW/PATH diagnostics
describe raw selection coordinates; final GFF and sequence output use reporting
coordinates.

On the 198-genome evaluation panel, the network refitted on all 43,013 eligible
occurrences retained all 223 detected references out of 227 (98.24% recall), with
four characterized misses, 33 unmatched calls and two unresolved calls unscored.
With the reduced library, refinement improves both-exact endpoints from 74 to 122
and leaves other families' calls unchanged. Mean absolute five-prime error increases
from 2.06 to 2.36 bases; three-prime error decreases from 2.17 to 1.28. Maximum
errors change from 46 to 36 bases and from 13 to 19 bases, respectively.

The larger 1,075-model library and its network reached 131 both-exact calls:
14 previously exact loci lose exactness and five become exact with the smaller
library. Among those 14 losses, median/maximum absolute errors are 5/36 bases at
the five-prime end and 0/18 bases at the three-prime end. All four separately added
references remain detected, with one both-exact rather than the larger library's
two. Adaptive mapping reduces reported model-verification alignments from 71,837
to 13,436 (81%) relative to the reduced-library baseline. These counts exclude
additional work inside boundary-feature calculation and are not wall-time speedups.
These are all-corpus coverage measurements against covariance-model reference
bounds, not independent generalization or physical transcript-end accuracy.
Unmatched calls have not all been established to be false positives.

Resources live in `resources/callgenes_TPP/` and `resources/callgenes_Glycine/`.
Selecting a family with missing or incompatible resources fails loudly. Omitted
options and explicit `tpp=f glycine=f` preserve ordinary production-family defaults.

## Core Gene Prediction

### CallGenes (CallGenes.java)
**Purpose**: Primary executable class for prokaryotic gene-calling and annotation pipeline.
- **Core Function**: Predicts and processes open reading frames (ORFs) in prokaryotic genomic sequences
- **Key Features**:
  - Multi-pass gene model refinement with iterative training
  - Supports multiple output formats: GFF, amino acid FASTA, rRNA FASTA
  - Configurable gene prediction parameters (length, scoring thresholds)
  - Parallel processing of input genome files
  - Extensive statistical tracking and reporting
- **Usage**: Comprehensive gene annotation for prokaryotic genomic data, including identifying coding sequences (CDS), detecting rRNA genes (16S, 18S, 23S, 5S), detecting transfer RNA (tRNA), and translating gene sequences

### GeneCaller (GeneCaller.java)
**Purpose**: A thread-safe gene caller that identifies and annotates genes within prokaryotic DNA sequences.
- **Core Function**: Detects and predicts genes (protein-coding and RNA) across six reading frames with advanced scoring and filtering mechanisms
- **Key Features**:
  - Handles both protein-coding (CDS) and RNA gene types (tRNA, 16S, 23S, 5S, 18S)
  - Supports multi-strand gene detection with configurable thresholds
  - Applies probabilistic gene model for sophisticated scoring
  - Breaks long ORFs at internal high-scoring start sites
  - Provides thread-local aligners for efficient parallel processing
- **Usage**: Used in prokaryotic genome annotation to systematically identify and characterize genes within DNA sequences

### Orf (Orf.java)
**Purpose**: Represents an Open Reading Frame (ORF), a DNA sequence segment with start and stop codons divisible by 3.
- **Core Function**: Models and scores potential genes within prokaryotic DNA sequences across six reading frames
- **Key Features**:
  - Handles multiple gene types (CDS, tRNA, rRNA: 16S, 18S, 5S, 23S)
  - Calculates sophisticated gene scoring using start/stop codon context, k-mer frequencies
  - Supports coordinate flipping for minus-strand genes
  - Provides detailed codon-level binary encoding
  - Tracks gene prediction path and path score for complex gene models
- **Usage**: Used in prokaryotic genome annotation to predict and score potential genes, enables advanced gene model training and refinement

## Gene Model Training and Statistics

### GeneModel (GeneModel.java)
**Purpose**: Store and analyze k-mer frequencies for prokaryotic gene starts, stops, and interiors
- **Core Function**: Performs statistical modeling of prokaryotic genomic sequences, tracking gene annotations and nucleotide composition across multiple gene types
- **Key Features**:
  - Processes multiple gene types: Protein-coding sequences (CDS), tRNA, rRNA (16S, 23S, 5S, 18S)
  - Supports multi-frame k-mer analysis for gene start and stop codon contexts
  - Tracks nucleotide frequency distribution and GC content
  - Handles both positive and negative genomic strands
  - Supports loading gene annotations from GFF files
- **Usage**: Train gene models for different prokaryotic clades or genomic conditions, analyze gene start/stop patterns across different gene types

### GeneModelParser (GeneModelParser.java)
**Purpose**: Static class for loading and parsing Probabilistic Gene Model (PGM) files containing k-mer frequency statistics.
- **Core Function**: Deserializes and reconstructs complete gene prediction models from specialized file format
- **Key Features**:
  - Handles parsing of complex PGM file headers with metadata about training sequences
  - Supports parsing of different feature types (CDS, tRNA, 16S, 23S, 5S, 18S rRNA)
  - Static factory method `loadModel()` for convenient model reconstruction
  - Robust line-by-line parsing with error checking and assertions
- **Usage**: Used for loading pre-trained gene prediction models in BBTools prokaryotic genome analysis pipeline

### PGMTools (PGMTools.java)
**Purpose**: Static utility class for manipulating Probabilistic Gene Model (PGM) files and merging gene models.
- **Core Function**: Provides methods for loading, merging, and processing probabilistic gene prediction models with advanced scaling and normalization capabilities
- **Key Features**:
  - Supports loading multiple PGM files with optional per-file multipliers
  - Implements model normalization based on training dataset size
  - Provides methods for merging and mixing gene models with configurable weighting
  - Handles various gene types: CDS, tRNA, rRNA (16S, 18S, 23S, 5S)
  - Supports file input/output with flexible overwrite and duplicate handling
- **Usage**: Used in prokaryotic genome analysis pipeline for combining gene prediction models from multiple sources, normalizing and scaling models based on training data

### AnalyzeGenes (AnalyzeGenes.java)
**Purpose**: Analyzes paired prokaryotic fna and gff files to calculate patterns in coding and noncoding frames, start and stop sites.
- **Core Function**: Process FASTA and GFF genomic files to extract k-mer frequency statistics for training prokaryotic gene prediction models.
- **Key Features**:
  - Parallel processing of multiple genome files using thread pools
  - Multi-threaded file reading and processing
  - Configurable processing parameters (verbose mode, ribosomal alignment)
  - Supports compressed (.gz) input files
  - Automatic FASTA/GFF file pairing
  - Performance statistics generation
- **Usage**: Command-line tool for bioinformatics researchers developing prokaryotic gene prediction models, particularly useful for analyzing multiple genome datasets simultaneously

## Statistical Framework

### StatsContainer (StatsContainer.java)
**Purpose**: Manages statistical tracking and aggregation for different genomic feature types with frame-specific k-mer analysis.
- **Core Function**: Collects and processes k-mer frequency statistics for various genomic feature regions (inner, start, stop) across multiple reading frames
- **Key Features**:
  - Supports multiple feature types: CDS, tRNA, 16S, 18S, 5S, 23S rRNA
  - Configurable k-mer size and frame analysis for different genomic regions
  - Tracks length statistics and distribution for features
  - Enables merging and scaling of statistical data from multiple sources
  - Provides serialization and deserialization of statistical models
- **Usage**: Used in prokaryotic genome annotation for training gene prediction models, aggregates statistical information about genomic feature characteristics

### FrameStats (FrameStats.java)
**Purpose**: Stores frame-relative k-mer counts for genomic feature analysis.
- **Core Function**: Tracks k-mer frequencies across multiple reading frames for specific genomic features
- **Key Features**:
  - Supports configurable context windows for k-mer pattern analysis
  - Handles multiple reading frames simultaneously
  - Validates k-mer frequencies across true and false positive feature instances
  - Provides advanced scoring mechanisms for genomic sequence feature detection
- **Usage**: Used in computational genomics to analyze sequence context around specific genomic features like coding start sites, stop sites, and ribosomal regions

### ScoreTracker (ScoreTracker.java)
**Purpose**: Creates a score tracker for collecting and aggregating statistical metrics for specific genomic feature types.
- **Core Function**: Accumulates scoring statistics for different genomic features like CDS, tRNA, ribosomal RNAs across multiple processing threads or analysis batches
- **Key Features**:
  - Tracks statistical metrics for gene start, stop, and inner scores
  - Supports multiple feature types (CDS, tRNA, 16S, 18S, 5S, 23S rRNA)
  - Enables merging statistics from multiple processing threads
  - Calculates average scores for gene start/stop contexts and k-mer frequencies
  - Computes approximate genic fraction of genome coverage
- **Usage**: Used in prokaryotic genome analysis pipelines to aggregate and report gene prediction performance metrics

## Ribosomal Sequence Processing

### MergeRibo (MergeRibo.java)
**Purpose**: Picks one ribosomal (16S or 18S) sequence per taxID for further processing.
- **Core Function**: Processes and selects the highest quality ribosomal sequence for each taxonomic ID
- **Key Features**:
  - Supports processing of both 16S and 18S ribosomal sequences
  - Uses advanced sequence alignment and scoring mechanisms
  - Supports multi-threaded processing of input files
  - Implements consensus sequence generation
  - Supports taxonomic tree-based sequence filtering
- **Usage**: Used in bioinformatics pipelines for selecting representative ribosomal sequences, helps consolidate multiple sequences from the same taxonomic group

### MergeRibo_Fast (MergeRibo_Fast.java)
**Purpose**: Picks one ribosomal (16S/18S) sequence per taxonomic ID using best alignment quality.
- **Core Function**: Deduplicates ribosomal sequences by selecting the highest-quality representative sequence for each taxonomic group.
- **Key Features**:
  - Multi-threaded ribosomal sequence processing
  - Consensus-based sequence selection (16S and 18S rRNA)
  - Alignment quality scoring via length × identity metric
  - Configurable input/output file handling
- **Usage**: Used in bioinformatics pipelines to reduce redundancy in ribosomal sequence datasets by selecting the most representative sequence for each taxonomic group

### SplitRibo (SplitRibo.java)
**Purpose**: Splits a mix of ribosomal sequences (such as Silva) into different files per type (16S, 18S, etc).
- **Core Function**: Classifies and routes ribosomal sequences to type-specific output files using alignment-based identification
- **Key Features**:
  - Multi-threaded ribosomal sequence processing
  - Two-phase alignment classification (broad consensus and clade-specific refinement)
  - Supports multiple ribosomal types: 16S, 18S, 23S, 5S, with mitochondrial and prokaryotic variants
  - Configurable minimum identity thresholds for sequence classification
  - Generates type-specific output files using replaceable filename pattern
- **Usage**: Processing large ribosomal sequence databases like Silva, separating mixed ribosomal sequence collections

### RiboMaker (RiboMaker.java)
**Purpose**: Makes a consensus ribosomal sequence using raw reads as input.
- **Core Function**: Generates consensus ribosomal sequences through advanced multi-threaded alignment and filtering processes
- **Key Features**:
  - Multi-threaded ribosomal sequence processing
  - K-mer based read filtering against reference sequences
  - Alignment quality evaluation using identity metrics
  - Supports paired-end read processing
  - Configurable alignment and filtering parameters
- **Usage**: Used in bioinformatics pipelines for generating high-quality ribosomal consensus sequences, processes raw sequencing reads to create representative ribosomal sequence

### FilterSilva (FilterSilva.java)
**Purpose**: Removes unwanted sequences from Silva database, particularly bacteria flagged as eukaryotes due to name misidentification.
- **Core Function**: Performs intelligent filtering of biological sequence data by taxonomic classification and header analysis
- **Key Features**:
  - Parses Silva database headers for taxonomic classification
  - Filters out misclassified bacterial and archaeal sequences from eukaryotic sets
  - Removes organellar sequences (Chloroplast, Mitochondria)
  - Supports concurrent read processing of FASTQ/FASTA files
  - Provides detailed processing statistics and performance tracking
- **Usage**: Used for cleaning and curating Silva taxonomic databases, ensuring accurate phylogenetic sequence collections

## Data Management

### ScafData (ScafData.java)
**Purpose**: Tracks and manages information about genomic scaffolds during gene annotation processing.
- **Core Function**: Provides a data structure for storing and manipulating scaffold-level genomic sequence information, supporting gene calling and annotation workflows
- **Key Features**:
  - Stores genomic sequence data with strand-specific annotations
  - Tracks frame annotations and open reading frame (ORF) processing state
  - Supports reverse complementing scaffold sequences
  - Manages separate containers for coding (CDS) and RNA features
  - Provides methods for subsequence extraction and strand orientation tracking
- **Usage**: Used in prokaryotic genome analysis to represent and process individual genomic scaffolds

### PFeature (PFeature.java)
**Purpose**: Represents a base implementation for genomic features with coordinate and strand management in prokaryotic sequence processing.
- **Core Function**: Provides an abstract base class for storing and manipulating genomic feature coordinates, with specialized strand processing capabilities
- **Key Features**:
  - Supports 0-based inclusive coordinate system for genomic features
  - Implements coordinate flipping between forward and reverse strands
  - Provides methods for calculating feature length and current strand
  - Implements Comparable interface for consistent feature ordering
  - Maintains scaffold name, start/stop positions, and strand orientation
- **Usage**: Used as a base class for different types of genomic features in prokaryotic genome analysis

### ProkObject (ProkObject.java)
**Purpose**: Static utility class for centralized gene type processing and validation in prokaryotic feature detection.
- **Core Function**: Provides comprehensive static methods for configuring, processing, and validating different prokaryotic gene types during genomic analysis
- **Key Features**:
  - Manages gene type constants for CDS, tRNA, rRNA (16S, 23S, 5S, 18S)
  - Handles k-mer loading for RNA gene validation
  - Manages consensus sequence loading and filtering
  - Provides boundary tolerance configuration for different RNA types
  - Supports configurable gene type processing flags
- **Usage**: Used as a central configuration and utility class in BBTools' prokaryotic genome annotation pipeline

## ncRNA shortlist measurement

`callgenes.sh ncrnapolicy=one_policy.tsv` applies a measured complete shortlist/
verifier tuple to each named shared-boundary family member. It reads the same
nine-column policy table as the sweep, requires exactly one policy and does not
enable callers. Rows override scalar shortlist controls for those members; other
families retain their settings. Missing, repeated or incompletely configured
loaded members fail before workers or output. Registry, loader and ordinary
worker construction preserve each packaged stopping policy when the option is
omitted. SAM-I uses the public fixed 238-model shortlist with patience16;
cobalamin RF01482/RF01689 keep the public fixed 100-model shortlists and legacy
stopping. RF00174 uses adaptive cap384/floor2 with patience8. These mapping
policies do not change family enablement or endpoint refinement.
Mapping K, seed and resource selection, identity cutoffs, endpoint models and
path scoring remain independent; the tRNA `patience=` flag keeps its meaning.

`dev/optimizencrnashortlist.sh` performs a measured local search over coupled model
shortlist controls. `cost=alignments` retains the existing logical-alignment
objective and resume contract; `cost=modelbases` uses summed aligned model length
and binds that unit into the saved contract. Both use the same TP target, optional
loss budgets and tie-breakers. The tool does not adopt production defaults.
Model bases omit query-window length and do not measure alignment cells or
whole-caller runtime. See the launcher's help for the exact selection order.

`dev/ncrnasweepfrontier.sh` reads completed native search metadata and a complete
original-input taxonomy selection ledger. It reports every saved policy and both
the TP/model-base and coverage/model-base frontiers, with exact search/round
provenance, apparent FP, CDS, first-gain and last-loss counts. It keeps different
panels, runtimes, references and caller settings separate. Repeated tuples merge
only when all reported summary and per-input counters agree; conflicts remain
visible. The taxonomy view describes representation within the supplied corpus,
with unknown/unclassified labels separate. It never opens genome/GFF payloads,
runs a caller, estimates missing counters or selects a production operating point.

`dev/summarizencrnasweep.sh` joins saved paired changes to reference support and
original-input coverage. Optional `taxonomy=` reports first/last detection counts
using saved full-corpus categories. `compare=left,right` adds exact lost/gained
reference keys, category totals, exclusive detections and subset flags for two
policies. Net TP differences do not establish set inclusion. This metadata view
does not compare complete emitted-call sets or infer biological causes of losses.

## Data Acquisition

### FetchProks (FetchProks.java)
**Purpose**: Crawls NCBI's FTP site to download bacterial genomes and annotations with phylogenetic sampling.
- **Core Function**: Parallel crawler that retrieves genomic and annotation files from NCBI RefSeq bacterial genome repositories
- **Key Features**:
  - Genus-based parallel processing with thread distribution
  - Configurable species-per-genus limit for balanced sampling
  - Quality-based assembly selection algorithm
  - Supports Candidatus organism classification
  - Configurable file renaming and sequence standardization
- **Usage**: Command-line tool for bulk bacterial genome retrieval with fine-tuned filtering and processing options
