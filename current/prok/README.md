# Prok Package - Prokaryotic Genome Annotation and Gene Prediction

BBTools prokaryotic package providing comprehensive tools for bacterial and archaeal genome annotation, gene prediction, ribosomal sequence processing, and statistical modeling of genomic features.

## ncRNA resource development

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
