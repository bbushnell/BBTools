#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified September 17, 2026

Description:  Finds orfs and calls genes in unspliced prokaryotes.
This includes bacteria, archaea, viruses, and mitochondria.
Can also predict 16S, 18S, 23S, 5S, 5.8S, LSU, ITS, and tRNAs.

Usage:  callgenes.sh in=contigs.fa out=calls.gff outa=aminos.faa out16S=16S.fa

File parameters:
in=<file>       A fasta file; the only required parameter.
out=<file>      Output gff file.
outa=<file>     Amino acid output.
out16s=<file>   16S output.
out18s=<file>   18S output.
outits=<file>   ITS1, ITS2, and combined ITS sequence output.  Implies its=t.
model=<file>    A pgm file or comma-delimited list.
                If unspecified a default model will be used.
stats=stderr    Stats output (may be stderr, stdin, a file, or null).
hist=null       Gene length histogram.
compareto=      Optional reference gff file to compare with the gene calls.
                'auto' will name it based on the input file name.

Formatting parameters:
json=false      Print stats in JSON.
binlen=21       Histogram bin length.
bins=1000       Maximum histogram bins.
pz=f            (printzero) Print histogram lines with zero count.

Genetic-code parameters:
transl_table=   Explicit NCBI table 4, 11, or 25 for CDS calling and translation.
                Omitted preserves historical ATG/GTG/TTG starts and translation.
                Explicit 11 uses its full NCBI initiation set. The selected code
                also governs start/stop candidates during multipass training.
codefile=      Alternative complete TSV: codon, amino_acid, start header;
                64 unique DNA codons, one canonical residue or *, start 0/1.
                Mutually exclusive with transl_table. Real terminal stops are
                omitted from proteins; truncated ends retain complete codons.
                Initial M requires a known complete allowed start. Synthetic
                contig-edge starts retain their elongation residue.
                Existing PGMs are not retrained by selecting a code; accuracy
                with a reassigned code requires separate evaluation.
                Explicit codes currently reject recode/detranslate modes.

rRNA alignment parameters:
rrnafallback=f  For 16S, 18S, 23S, and 5S: after the universal rRNA consensus
                fails alignment, try later consensuses in shipped-file order.
                Default false preserves historical universal-only behavior.



Taxonomy parameters:
taxonomy=t      Use QuickClade to classify the input and select a phylum-
                specific gene model for improved accuracy.  Intended for
                ISOLATES (one organism per file): a metagenome would be
                misclassified to a single phylum, so set taxonomy=f for
                mixed samples.  Requires a QuickClade server; if it is
                unreachable, callgenes falls back to the general model with
                a warning (it never fails on this account).
percontig=f     Classify each contig separately (for metagenomes).
                Default is per-file (classify once for all contigs).
taxaddress=     QuickClade server address.  Default: refseq.
local=f         Use QuickClade's LOCAL reference database (large multi-GB
                files under resources/) instead of the network server. This
                is the REVERSE of QuickClade's own default (which favors
                local, falling back to the server only if local files are
                missing) -- callgenes stays lightweight by default, and
                local=t opts INTO the heavy local files if present (falls
                back to the server, with a warning, if they're not).
server=t        Alias for local=f (server=f is the same as local=t).
                With taxonomy=t, the detected domain, phylum, and which path
                (local/server) served the classification are written as
                "##Domain"/"##Phylum"/"##ClassificationSource" GFF header
                lines, and are available to an in-process Java caller via
                CallGenes.lastDetectedDomain()/lastDetectedPhylum()/
                lastClassificationSource().

Experimental CM verification (default off; tRNA and 5S only):
cm=f           Opt into exact local-CYK best-interval search for small families
                around the candidate extent; retain caller endpoints and
                path scores. Each candidate gets its own best interval with
                overlap at least half the shorter of the candidate and CM hit.
                Existing candidate and borderline filters remain unchanged.
                Small-family defaults: +/-12 search, conditional tRNA-Sec,
                retry the full existing window, then add +/-64 on rejection.
                With euk18s=t, 18S still uses the same fast alignment as cm=f;
                no 18S covariance model is loaded or attached.
cmlocal=       Explicit t selects local-model CYK for the enabled small families;
                explicit f selects global scoring and disables omitted retries.
                Omitted: local for small families. Requires cm=t
                when true. Changes model begin/end scoring, not candidate windows
                or caller endpoints. Small families use exact interval search.
                Cutoffs are unchanged provisional model GA values, not calibrated
                operating points for local scoring. Diagnostic routes gain a _local marker.
cmlocalwide=   Defaults to true for enabled small-family local scoring; f disables.
                For small families, retry a failed local decision over the full
                existing padded window. Requires cm=t cmlocal=t. Spliced tRNA
                candidates also search the unspliced genomic ORF plus its existing
                padding; any mature-product pass is retained. No added padding,
                endpoint changes or18S changes. Uses half-of-the-shorter overlap.
                Wide diagnostic routes contain _wide_local; all attempts retain
                their costs, with one selectedModel row across both windows.
cmlocalextend= Defaults to 64 when local/wide retries are enabled, otherwise 0.
                Explicit 0 disables extension. After all local-window attempts reject a small-family
                candidate, retry once with this many additional bases per side
                of the full genomic window, clipped to the sequence boundaries.
                Requires cm=t cmlocal=t cmlocalwide=t. Existing passes, candidate
                endpoints, association and model cutoffs remain unchanged.
                Extended diagnostics use _extended_local and retain all costs.
cmtrnaextendedgreedy=t  Nonoverlapping local-hit association only
                for the extended tRNA retry. When omitted, enabled only with
                cm=t, local/wide scoring and positive cmlocalextend.
                Explicit f retains the earlier constrained-span control.
                With cmdiag, CM_GREEDY_HIT stderr rows retain each extracted hit.
cmsmallgreedy=   Defaults to true with cm=t and local small-family scoring.
                Explicit f disables this rule; cmtrnaextendedgreedy still
                independently controls the extended tRNA retry.
                Greedy association covers all tRNA/5S local windows,
                including spliced mature products and every retry. Association
                uses the sequence supplied to the CM; diagnostics retain genomic
                coordinates. Includes extended retries even if the narrower
                cmtrnaextendedgreedy flag is false. Requires cm=t cmlocal=t;
                Global CM mode is unchanged; disabled with cm=f.
euk18sjoinedproposals=  Experimental replay of frozen joined-window geometry.
                Requires euk18s=t, rolling PacBio, raw alignment ends,
                and minlen=916. Ordinary candidate rules remain unchanged;
                joined spans use the same identity table and final caller path.
euk18slivepositions=  Experimental live joined discovery from the shared 17-mer
                hit stream and an exact-consensus-bound position resource.
                Mutually exclusive with saved proposal replay; same fast-only
                controls. A joined call must realize a deletion of at least 256
                bases and retain at least 916 aligned core bases.
euk18sjoinedsupport=1  Select side-support 1 or 2 for either joined source.
euk18sjoined=t  Explicit control for joined discovery; t requires rolling PacBio.
                Omitted Quantum/native-MSA calls use ordinary windows only.
                With t and no explicit position path, load euk18S_seed_positions.tsv.
                Use f in both backends for a matched aligner-only comparison.
cmmodeldir=    Directory containing RF00005.cm and RF01852.cm (tRNA and tRNA-Sec)
                and RF00001.cm (5S). RF01960 is not used by this caller.
                Any model passing its own cutoff admits a
                candidate; diagnostics select the highest-scoring passing model
                among those evaluated.
                Without this, look in resources/. Models are loaded only for
                enabled families and cm=t. Model GA cutoffs are provisional,
                not calibrated for this scoring mode. IUPAC residues use equal-
                weight averages of canonical emission bit scores; gaps are invalid.
cmwindowslack=12 Limit the search to the candidate extent plus this many bases
                per side, intersected with the original padded window. Applies
                to small-family CM verification; -1 restores the full window.
cmsecfallback=t Evaluate RF01852 only if RF00005 fails its own cutoff.
                Set false to evaluate both tRNA models for every candidate.
cmmaxcells=100m Maximum DP cells per model/window (one float score plus one int
                traceback choice per cell) for small-family CM verification.
                This is not a whole-JVM memory bound and does not apply to 18S.
cmdiag=        Optional TSV of every scored candidate, including rejected ones,
                source route, caller/window/CM extents, score and compute cost.
                All model attempts are retained; selectedModel marks one row
                per candidate. Unrestricted best hits diagnose neighboring loci.
                These candidates may still lose final gene-path competition.
                Exact 18S rows have an exact_extent route suffix, NA unrestricted
                search fields and zero interval ties; a finite parse consumes
                its supplied extent and does not refine the caller endpoints.

tRNA detection parameters:
scavengeonly=t  Call tRNAs with the kmer-guided scavenger only (finds tRNAs at
                conserved long-kmer positions).  This is the shipped default;
                set scavengeonly=f to use the PGM candidate generator instead.
scavenge=f      Run the scavenger IN ADDITION to the PGM candidate generator
                (augment mode) rather than replacing it.
trnaintron=t    Intron-aware pass: splice candidate spans before verification
                to recover intron-containing (mainly archaeal) tRNAs.
                (Not used while scavengeonly is on.)
mintrnakhits=1  Min conserved tRNA long-kmers a candidate must carry before it
                is aligned; a cheap pre-filter that rejects non-tRNA windows.
                0 disables the pre-filter.

tRNA alignment parameters:
trnaalign=t     Align predicted tRNAs to a consensus library to verify
                and annotate them.  Reduces false positives and adds
                anticodon/amino-acid annotation.  Uses built-in library
                by default; override with trnalib= and trnamodel=.
trnalib=<file>  Custom tRNA consensus library (fasta).
trnamodel=<file> Custom tRNA HBM model file.
indexk=7        Kmer length for the library shortlist index.  Longer is more
                selective (fewer models aligned per candidate); shorter is more
                permissive.  7 is the shipped default.
indextopn=60    Max library models aligned per candidate (search breadth).
indexminhits=12 Min shared index-kmers for a model to enter the shortlist.
                Fixed fallback; used only when adaptiveminhits=f.
adaptiveminhits=t  Adapt the shortlist cutoff per candidate instead of the
                fixed indexminhits (the shipped default).  Cutoff = ceil(max(
                adaptfloor, adapttopfrac*maxSharedKmers, adaptqfrac*queryKmers)).
adaptfloor=11   Absolute floor for the adaptive cutoff (the constant term).
adapttopfrac=0.48  Adaptive cutoff as a fraction of the best model's shared-
                kmer count.
adaptqfrac=0.072   Adaptive cutoff as a fraction of the candidate's kmer count.
patience=20     Stop aligning after this many models without improvement
                once a passing hit has been found (with earlyexit).
earlyexit=t     Enable the patience-based early exit.
idpass=0.75     Min alignment identity to accept a tRNA.
hbmpass=0.75    Min HBM model score to accept a borderline tRNA.
acextract=t     Extract the anticodon directly from each verified tRNA's
                structure (anticodon loop position projected through the
                alignment).  Adds an anticodon: attribute to the GFF.
acvalidate=12   Min structural score (anticodon stem + U33 + purine-37,
                max 15) to trust an extracted anticodon; failures fall
                back to model-name annotation.
acmargin=3      Min score margin between the best and runner-up anticodon
                register; ambiguous positions fall back to the model name.
trimtrna=t      Trim verified tRNA boundaries to the consensus alignment
                extent, removing scanner slop for accurate coordinates.
maxtrna=120     (Experimental) Raise the tRNA candidate length cap, enabling
                relaxed length scoring for over-length candidates such as
                intron-containing archaeal tRNAs.  Measured neutral: unspliced
                candidates still fail alignment verification.

Conserved ncRNA engineering parameters:
euk18s=f       Enable the resource-bound development 18S recipe: seven cleaned
                models, 385 seeds, distinct-2 support, voted windows and rolling MSA.
euk18sconsensus=   Consensus FASTA; requires matching euk18smodelcutoffs= table.
euk18sseeds=       Forward K17 FASTA; omitted uses euk18S_17mers.fa.gz.
euk18srawends=alignment  Or vote; omitted offsets use euk18S_seed_offsets.tsv.
euk18saligner=pacbio  Alternative quantum retains the earlier primary aligner.
                PacBio MSA primary with sqrt-discounted deletion-run identity.
                Keeps model cutoffs and optional local rescue unchanged.
                Requires alignment ends and endpoint NN off. Small-family cm=t
                may run concurrently. Native
                out-of-window alignments are rejected without clamping.
                With ncrnadiag, also runs and counts a diagnostic Quantum
                alignment; supplemental NCRNA_MODEL_PACBIO stderr rows retain
                both identities, native score/endpoints/padding, separate costs
                and matrix payload. Diagnostics are off unless requested.
euk18smsadel4=-2  Development-only PacBio penalty at deletion residues21-80.
euk18smsaengine=rolling  Uses six score rows plus
                compact traceback, retaining PacBio scores and paths. Requires
                euk18saligner=pacbio. Explicit native selects the full matrices.
                Rolling matrix-payload diagnostics count rows plus ancestry;
                neither engine's payload is a whole-JVM memory requirement.
euk18smsadel5=-1  Long-deletion charge every fourth residue beyond80.
euk18smsains4=-8  Long query-insertion extension penalty (genomic inserts use DEL).
euk18smsains=-205 Insertion-open penalty; euk18smsadel=-292 deletion-open.
euk18smsasubscale=0.66  Scale SUB/SUBR/SUB2/SUB3 together, rounding magnitudes
                to nearest integer with half ties upward. Range (0,1].
                All cost flags require euk18saligner=pacbio. Initial penalty
                domain is [-292,0]. Default substitution scale is 66%; explicit
                1 restores the former baseline. Other penalties retain native constants.
euk18svotewindows=t euk18svotebeforepadding=t  Use repaired-core votes.
euk18swindowpad=3109  Fallback window radius. Mapping uses K9.
euk18sflank=100  Retained explicit engineering control.
euk18sminlen=916 euk18smaxlen=4000 euk18sseedminhits=2
euk18sseeddistinct=t euk18scollapsefrac=0.9
euk18smodelclip=rescue euk18srescueid=0.66  Retry primary rejection only.
euk18sindexminhits=96 euk18stopn=<uncapped>
euk18sindexmargin=121  Nonnegative values replace topN with an inclusive shared-K9
                score gap below the best model in each window. Keeps all positive
                scores >= best-minus-margin and >= indexminhits; no model-count cap.
                Query words are deduplicated; model occurrences retain multiplicity.
                Do not also specify euk18stopn. -1 keeps the existing topN policy.
                First accepted alignment exits.
euk18sindexfrac=0.8  B in [0,1]; also require shared K9 >= ceil(B * best).
                Combined cutoff is max(indexminhits, best-indexmargin, ceil(B*best));
                omit the margin term when indexmargin=-1. B=0 preserves prior behavior.
                Nonzero B uses all models, so do not also set euk18stopn.
euk18svoteslack=60 euk18svotemaxsd=5 euk18svoteseedsha80=<20hex>
euk18sendpoint=f  Optional resource-sized refinement; requires explicit
                euk18sendpointtables= and euk18sendpointnets= directories.
euk5s=f        Independently enable eukaryotic 5S with two seed hits per window.
                Needs neither ncrna=t nor rrna17=; leaves legacy flags unchanged.
                All calls compete in the shared DP by score and overlap penalties.
                Shortlist all models (currently 7); load model-specific cutoffs
                (currently .68/.68) and refine endpoints with packaged networks.
euk5sconsensus= Explicit euk5S consensus FASTA for a frozen resource comparison.
                Requires euk5s=t or an euk5S rrna17 profile. Seed overrides use
                ncrnafamily=euk5S ncrnakmers=; other families are unchanged.
euk5sendpoint=t Refine enabled euk5S calls; f retains raw alignment endpoints.
euk5sendpointtables=  Override the directory containing TABLES.tsv and tables.
euk5sendpointnets=    Override the directory containing MANIFEST.tsv and networks/.
                Defaults resolve resources/euk5S_endpoint/. Missing or incompatible
                resources fail before calling; start and stop may use different k.
euk5smodelcutoffs=  Override resources/euk5S_model_cutoffs.tsv; f uses scalars.
                Each model needs family/model/idpass/idborderline TSV columns.
                Explicit model tables take precedence over scalar identity flags.
                Targeted scalar overrides bypass the default table, retaining .68
                for either unset scalar. euk5s=f loads none of these resources.
ncrna=f         Enable generic conserved-ncRNA families (off by default).
                Alias: generalncrna=.
rrna17=f        Enable an integrated 17-mer rRNA profile: prok loads
				16S/23S/prokaryotic-5S; euk loads 18S/eukaryotic-5S.
				euk18s and euk5s explicitly load either eukaryotic family alone.
				Families in one
				profile share one seed scan per contig and strand.
                Requires 16s=f 23s=f 5s=f 18s=f so legacy PGM and 17-mer rRNA
				calls cannot mix.
                Alias: rrna17mer=.
rrna17prok5shbm=t  Enable prok5S HBM rescue for alignments with identity
                from 0.60 through 0.80.  Set false for a matched detection
                ablation; requires rrna17=prok and leaves the full profile on.
rrna23sendpoint=off  23S endpoint refinement inside rrna17=prok is disabled.
                The discarded g28/g16 candidate-scoring networks are not shipped;
                explicitly requesting either fails if its resources are missing.
v2refresh=f     Recheck claimed loci before each ncRNA window when enabled.
tmrna=f         Add the experimental tmRNA family; requires ncrna=t (or
                generalncrna=t).  Endpoint refinement follows the global
                ncrnaboundarynet flag below.
sixs=f          Add the experimental paired 6S/SsrS families (RF00013 and RF01685);
                requires ncrna=t (or generalncrna=t).  Aliases: ssrs=, 6s=.
                Both families must load together, with complete consensus, HBM,
                and kmer resources and index-aligned consensus/HBM entries.
                If ncrnaboundarynet=t, complete index-aligned boundary resources
                are also required for every loaded family; missing or misaligned
                resources fail loudly.
r58lsu=f        Add the bundled eukaryotic 5.8S and LSU families.  This implies
                ncrna=t; the family remains off unless r58lsu=t is specified.
                Explicit resource overrides below remain available.
                The global ncrnaboundarynet gate remains incompatible with this
                pair; use the R58-specific endpoint gate below.
r58ncrnaboundarynet=f  Refine only R58 endpoints with its boundary network and
                9-mer tables.  Requires r58lsu=t; LSU remains boundary-NN-off.
r58ncrnaboundary3primeonly=f  Hold the aligned R58 5' coordinate fixed and
                refine only the 3' endpoint.  Requires r58ncrnaboundarynet=t.
r58ncrnaboundarymarginstop=0  Require this confidence improvement before moving
                the R58 3' endpoint.  The validated candidate value is 0.01.
r58ncrnaboundarymarginstart=0  Corresponding 5' confidence margin.  The current
                validated 3'-only configuration leaves this unused.
its=f           Derive ITS1, ITS2, and combined ITS annotations from compatible
                18S, 5.8S, and LSU calls.  outits=<file> also enables this.
s18=f           Add the experimental generic 18S development family (pilot);
                requires ncrna=t and explicit s18kmers/s18consensus/s18models.
                Distinct from 18s=t, which enables the legacy PGM-path 18S caller.
                This is a sweep harness, off by default, not a release caller.
ncrnaboundarynet=f  Refine generic ncRNA endpoints with the family-specific
                    boundary networks.  Requires ncrna=t (or generalncrna=t).
                    When enabled, every loaded family must have complete,
                    index-aligned boundary resources; missing or misaligned
                    resources fail loudly.
ncrnafamily=    Select rnasep, srp_small, srp_large, tmrna, sixs_rf00013,
                sixs_rf01685, r58, lsu, s18, or euk5S as the generic sweep
                overrides below; requires ncrna=t (or generalncrna=t).  A tmrna
                target also requires tmrna=t; a sixs_* target also requires
                sixs=t; an r58 or lsu target also requires r58lsu=t; an s18
                target requires s18=t.  euk5S instead requires euk5s=t,
                rrna17=euk or rrna17=euk5s; it does not require ncrna=t.
                Accepted aliases include rnase-p, rnase_p, srp-small, srpsmall,
                srp-large, srplarge, tm-rna, tm_rna, and ssra (case-insensitive).
                Other resource-available families remain loaded and compete
                normally; this is not an isolation switch.
ncrnakmers=     Explicit conserved-kmer fasta for the selected ncrnafamily=;
                requires that family's enabled profile and exactly one target.
ncrnaseedminhits=2  euk5S conserved 17-mer occurrences required per window.
                Standalone euk5s=t requires at least 2; rrna17 profiles retain
                explicit one-hit engineering controls for baseline comparisons.
ncrnaindexk=    Targeted family mapping k (1-15); seed k remains unchanged. Larger values
                allocate exponentially larger indices; 6/7 are tuning choices.
ncrnaindextopn=  Maximum models shortlisted per candidate (>=1).
ncrnaindexminhits=  Shared mapping-kmer cutoff (>=0). A multi-model index keeps
                its best model below this cutoff unless strict filtering is on.
                Zero bypasses mapping only for a one-model library.
ncrnastrictindexcutoff=f  Reject euk5S windows below the mapping cutoff.
ncrnarankedfallback=f  Stop after the first identity-passing euk5S alignment
                whose span meets the requested maximum length.
ncrnamaxlen=2147483647  Maximum inclusive euk5S alignment/call span (>=90).
                Does not cap seed-window size.  Overlong models are skipped
                before acceptance, allowing later shortlisted models to pass.
                Seed minhits, strict filtering, ranked fallback and maxlen
                require ncrnafamily=euk5S. Index controls support every listed
                ncrnafamily target. Omitted values
                preserve packaged behavior; no calibrated cutoff is implied.
                euk5S has no HBM: ncrnahbmpass is rejected for this target.
ncrnaadaptiveminhits=  Enable adaptive mapping cutoff for the selected family.
ncrnaadaptfloor=       Minimum adaptive cutoff (finite, >=0).
ncrnaadapttopfrac=     Best-model shared-count multiplier (finite, >=0).
ncrnaadaptqfrac=       Unique query-kmer multiplier (finite, >=0).
                Aliases: ncrnamapk, ncrnafixedminhits, ncrnaadaptiveshortlist.
                Prefix any index/adaptive flag above with r58 or lsu for an
                independent family override (e.g. r58ncrnaindexk=7). These
                require r58lsu=t and take precedence over targeted overrides.
r58ncrnaminlen=140  R58 minimum candidate length (1-1000); alias r58minlen.
                Only R58 changes; omitted retains 140. Requires r58lsu=t.
ncrnadiag=     Optional stage-diagnostic TSV path; alias r58diag. Reports all
                enabled ncRNA families, before final path selection. Omitted
                produces no diagnostics and adds no alignments.
rnasepkmers=    Explicit RNase-P conserved-kmer fasta; requires ncrna=t.
srpsmallkmers=  Explicit SRP-small conserved-kmer fasta; requires ncrna=t.
srplargekmers=  Explicit SRP-large conserved-kmer fasta; requires ncrna=t.
tmrnakmers=     Explicit tmRNA conserved-kmer fasta; requires ncrna=t and tmrna=t.
sixsrf00013kmers=  Explicit RF00013 16-mer fasta; requires ncrna=t and sixs=t.
sixsrf01685kmers=  Explicit RF01685 17-mer fasta; requires ncrna=t and sixs=t.
r58kmers=       Optional R58 17-mer fasta override; requires r58lsu=t.
lsukmers=       Optional LSU 17-mer fasta override; requires r58lsu=t.
r58consensus=/r58models=  Optional R58 consensus/HBM overrides; require r58lsu=t.
lsuconsensus=/lsumodels=  Optional LSU consensus/HBM overrides; require r58lsu=t.
r58boundarynet=           Optional R58 boundary-network override; requires
                          r58ncrnaboundarynet=t.
r58boundarystarttable=    Optional R58 start 9-mer table override; requires the
                          R58-specific boundary gate.
r58boundarystoptable=     Optional R58 stop 9-mer table override; requires the
                          R58-specific boundary gate.
s18kmers=       Explicit 18S development 17-mer fasta; requires ncrna=t and s18=t.
s18consensus=/s18models=  Explicit 18S development consensus/HBM; require s18=t.
tmrnaconsensus= Explicit tmRNA consensus fasta; requires ncrna=t and tmrna=t.
tmrnamodels=    Explicit aligned tmRNA HBM models; requires ncrna=t and tmrna=t.
                Generic sweep overrides below require exactly one ncrnafamily=;
                family-specific overrides do not select or isolate a family.
ncrnaidpass=    Alignment identity that accepts without HBM rescue.
ncrnaidborderline=  Lowest identity eligible for HBM rescue.
ncrnahbmpass=   HBM score needed to rescue a borderline alignment.
ncrnascorea=    Constant A in A+B*length*identity^2.
ncrnascoreb=    Multiplier B in A+B*length*identity^2.
rnasepscorea/b= Per-family A/B overrides; require ncrna=t.
srpsmallscorea/b= Per-family A/B overrides; require ncrna=t.
srplargescorea/b= Per-family A/B overrides; require ncrna=t.
tmrnascorea/b=  Explicit tmRNA A/B values (provisional shipped default 5/3);
                require ncrna=t and tmrna=t.  Family score overrides do not
                require ncrnafamily=, may be combined in one run, and take
                precedence over generic ones.
ncrnacollapsefrac=  Overlap fraction used to collapse candidate windows.
ncrnawindowpad= Candidate padding around conserved-kmer hits.
rnaseppad=     RNase-P-only candidate padding override; requires ncrna=t.
srpsmallpad=   SRP-small-only candidate padding override; requires ncrna=t.
srplargepad=   SRP-large-only candidate padding override; requires ncrna=t.
tmrnapad=       tmRNA-only padding override (provisional shipped default 370);
                requires ncrna=t and tmrna=t.
sixsrf00013pad= RF00013-only padding override (provisional shipped default 250);
                requires ncrna=t and sixs=t.
sixsrf01685pad= RF01685-only padding override (provisional shipped default 125);
                requires ncrna=t and sixs=t.
                When omitted, shipped defaults are rnasep=320, srp_small=75,
                srp_large=250, tmrna=370, RF00013=250, and RF01685=125.
                These flags override only candidate-window padding; run one
                family and one varied dimension at a time.

Advanced tRNA candidate-generation thresholds (rarely changed):
trnaregion=20   Region-open score cutoff.
trnacand=36     Composite candidate score cutoff (tRNA is very sensitive here).
trnastart=2.4   Start point-model score cutoff.
trnastop=1.5    Stop point-model score cutoff.
trnainner=2.2   Average inner-kmer score cutoff.

Other parameters:
minlen=60       Don't call genes shorter than this.
trd=f           (trimreaddescription) Set to true to trim read headers after
                the first whitespace.  Necessary for IGV.
merge=f         For paired reads, merge before calling.
detranslate=f   Output canonical nucleotide sequences instead of amino acids.
recode=f        Re-encode nucleotide sequences over called genes, leaving
                non-coding regions unchanged.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
"
}

if [ -z "$1" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	usage
	exit
fi

resolveSymlinks(){
	SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		RESOLVED_SCRIPT="$(readlink "$SCRIPT")"
		SCRIPT="$RESOLVED_SCRIPT"
		if [ "${SCRIPT#/}" = "$SCRIPT" ]; then
			RESOLVED_SCRIPT="$DIR/$SCRIPT"
			SCRIPT="$RESOLVED_SCRIPT"
		fi
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=6g" "--xms=6g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prok.CallGenes $@"
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" prok.CallGenes "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"
