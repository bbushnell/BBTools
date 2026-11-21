# JGI Package - BBTools Comprehensive Utilities Collection

The JGI package contains utilities for sequencing data analysis, read processing, quality control, and bioinformatics workflows.

## Package Overview

This package represents one of the largest collections in BBTools, offering utilities that span:
- Sequencing read quality control and preprocessing
- Assembly analysis and statistics
- Coverage analysis and normalization 
- Deduplication and clustering algorithms
- Primer and adapter sequence handling
- Barcode processing and demultiplexing
- Contamination detection and removal
- Statistical analysis and reporting

## Utility Classes

#### RenameReads (RenameReads.java)
**Purpose**: Utility for renaming and processing sequencing read files with multiple renaming strategies
**Core Function**: Transforms read identifiers using various renaming techniques like insert size, coordinates, or custom prefixes
**Key Features**:
  - Support for paired-end read renaming
  - Multiple renaming modes (by insert, trim, coordinates)
  - Quality score quantization
  - SRA header fixing
  - Flexible prefix and suffix addition
**Usage**: Used in bioinformatics pipelines to standardize read identifiers across different sequencing platforms and experiments

#### BBMergeOverlapper.java (BBMergeOverlapper.java)
**Purpose**: Perform advanced read merging and overlap detection for paired-end sequencing reads
**Core Function**: Implements multiple sophisticated algorithms to determine optimal read overlap, merge reads with quality-aware error correction
**Key Features**:
  - Quality-weighted overlap detection algorithm
  - Multiple overlap scoring methods (ratio-based, quality-based)
  - Supports variable read lengths and multiple insertion sizes
  - Handles ambiguous overlaps with detailed statistical tracking
  - Configurable mismatch and minimum overlap thresholds
  - Probabilistic error correction during read merging
**Usage**: Critical component in paired-end read processing for sequencing data cleanup, error correction, and read reconstruction in bioinformatics pipelines

#### AssemblyStats3 (AssemblyStats3.java)

`stats3.sh` writes a compact nine-column row per input: filename, total length,
nonempty sequence-record count, GC fraction, maximum length, and total bases in
records at least 5k, 10k, 25k and 50k bases long. It uses `Assembly` and an `IntList`;
it does not calculate N50/N90 or split scaffolds at gaps. `stats.sh` continues to use
`AssemblyStats2` for its more detailed report.

`reads=` (alias `maxreads=`) limits nonempty sequence records in input order,
separately for each file, before sorting lengths. Negative is unlimited (the default);
zero produces an empty statistics row. Empty header records do not count, and the
loader's existing acceptance of sequence before the first header is retained.

#### Assembly (Assembly.java)
**Purpose**: Parses and analyzes FASTA genome assembly files to extract sequence statistics.
**Core Function**: Loads FASTA file, computes nucleotide counts, tracks contig lengths, and calculates genome composition metrics.
**Key Features**:
  - Reads FASTA files with support for various base types (DNA/RNA)
  - Tracks total assembly length and contig sizes
  - Computes nucleotide base composition (GC content)
  - Supports case-insensitive base detection
  - Handles extended IUPAC base codes
**Usage**: Used in bioinformatics pipelines to analyze genome assembly characteristics and perform preliminary sequence statistics.

#### AddAdapters.java (AddAdapters.java)
**Purpose**: Legacy generator of insert-unaware adapter-contaminated read data, with an optional length-grading mode
**Core Function**: Replaces selected read regions with adapter sequence and random filler, preserving read lengths and quality arrays and recording original/expected lengths in generated headers
**Key Features**:
  - Selects a 5' prefix or 3' suffix without using insert size
  - Mates may share the sampled replacement position while independently choosing adapter entries
  - Retains undefined input bases; optional adapter-region substitutions use existing qualities or Q30 when quality is null
  - Optionally adds reverse-complement adapter entries
  - `minlen` classifies expected-remainder statistics and never filters generation output
  - Legacy grade mode compares surviving annotated lengths, without checking sequence identity or accounting for discarded reads
**Usage**: The deprecated `addadapters.sh` generates data for subsequent external trimmer evaluation; it does not run trimming tools. Use `randomreads.sh` when paired adapter placement must respect insert size.

#### AdjustHomopolymers.java (AdjustHomopolymers.java)
**Purpose**: Applies controlled expansions or contractions to homopolymer runs in sequence reads
**Core Function**: Adjusts each completed run, including the final run, using exact equality of the bytes supplied by input/validation
**Key Features**:
  - Advertised sequence inputs are single-end or paired-end FASTA/FASTQ; editing uses one caller with background stream I/O
  - Adjusts A/C/G/T/U in either case, with case included in byte equality; undefined runs remain unchanged
  - Requires a finite, nonzero rate at least -1: `-1` removes defined runs exactly, other rates use float multiplication and int truncation toward zero, and values above 1 are allowed
  - Expansion copies the run's final original quality; contraction drops qualities with their bases; null qualities remain null
  - Checks run ownership and projected builder capacity before applying an adjustment; allocation success for huge reads is not guaranteed
  - Requires sequence bases and does not recalculate alignment metadata; mapped SAM/BAM consistency is outside the reviewed correctness scope
**Usage**: Use `adjusthomopolymers.sh` with the desired rate for controlled run-length changes, consulting its help for input/output options.

#### AssemblyStatsWrapper (AssemblyStatsWrapper.java)
**Purpose**: Wrapper for processing multiple assembly statistics files with customizable input handling
**Core Function**: Dynamically processes input files through AssemblyStats2, managing file parsing, garbage collection, and output generation
**Key Features**:
  - Supports single and multiple input file processing
  - Dynamically handles comma-separated file lists
  - Manages memory by forcing garbage collection between large files
  - Configures output headers and appending for multi-file processing
  - Flexible parameter handling for assembly statistics generation
**Usage**: Command-line tool for generating comprehensive assembly statistics across multiple genome or sequence assembly files

#### BBMask (BBMask.java)
**Purpose**: Masks genomic sequences by identifying and replacing low-complexity or repetitive regions with 'N' characters.
**Core Function**: Advanced sequence masking algorithm that identifies and masks low-complexity, repetitive, or mapped genomic regions.
**Key Features**:
  - Supports multiple masking modes: repeat masking, low-complexity masking, and SAM-based masking
  - Configurable k-mer size and window parameters for detecting repeats and low-complexity regions
  - Multithreaded processing for improved performance
  - Optional entropy-based masking for sequence complexity detection
**Usage**: Preprocessing genomic sequences to remove low-quality or repetitive regions before further analysis, such as sequence alignment or assembly

#### AssemblyStats2 (AssemblyStats2.java)
**Purpose**: Comprehensive analysis of genomic assembly sequence statistics
**Core Function**: Processes FASTA and FASTQ files to extract detailed assembly metrics
**Key Features**:
  - Calculates scaffold and contig length distributions
  - Generates GC content histograms and statistics
  - Supports configurable length cutoffs and N-break detection
  - Outputs results in multiple format options (including JSON)
  - Handles both FASTA and FASTQ input formats
  - Tracks base composition (A, C, G, T, N) across assemblies
**Usage**: Used for comprehensive genomic assembly quality assessment and reporting, typically in bioinformatics pipelines for genome sequencing projects

#### BBDuk (BBDuk.java)
**Purpose**: Sequence processing tool for separating, trimming, or masking sequences based on k-mer matching against reference sequences
**Core Function**: Performs advanced sequence filtering and manipulation using k-mer-based reference comparison
**Key Features**:
  - Supports Hamming and edit distance sequence matching
  - Handles k-mer sizes from 1-31 with emulation for K>31
  - Configurable trimming modes (left, right, tips)
  - Supports sequence masking and filtering
  - Flexible reference-based sequence processing
  - Handles paired-end and single-end sequencing data
**Usage**: Used in bioinformatics preprocessing to clean, filter, and prepare sequencing reads for downstream analysis

#### CalcUniqueness (CalcUniqueness.java)
**Purpose**: Analyzes k-mer uniqueness patterns in sequencing data using parallel hash table techniques.
**Core Function**: Calculates unique and non-unique k-mer percentages across sequencing reads through interval-based sampling.
**Key Features**:
  - Multi-threaded k-mer uniqueness analysis for single and paired-end reads
  - Parallel 31-table hash distribution for efficient k-mer tracking
  - Quality-aware k-mer filtering with error probability thresholds
  - Interval-based and cumulative statistical reporting modes
  - Optional spike-fixing algorithm to smooth uniqueness percentage curves
**Usage**: Used in genomic sequencing analysis to assess the novelty and complexity of sequencing data by tracking unique k-mer frequencies.

#### BBQC (BBQC.java)
**Purpose**: Wrapper for BBDukF, BBMap, and BBNorm to perform comprehensive quality control and artifact removal for sequencing reads.
**Core Function**: Executes a multi-stage read processing pipeline including adapter trimming, quality filtering, artifact removal, and optional normalization
**Key Features**:
  - Supports DNA and RNA read processing
  - Removes adapters, low-quality bases, and contaminating sequences
  - Optional human read removal
  - Quality and length-based read filtering
  - Error correction and read normalization
**Usage**: Used in bioinformatics preprocessing to clean and standardize sequencing data before further analysis

#### CheckStrand2 (CheckStrand2.java)
**Purpose**: Determines the strandedness of RNA-seq reads through multiple computational methods
**Core Function**: Performs comprehensive strand orientation analysis of sequencing reads using k-mer sketching, stop codon analysis, poly-A/T tail detection, and gene calling.
**Key Features**:
  - Supports multiple RNA-seq strandedness detection techniques (depth analysis, stop codon analysis, poly-A/T counting)
  - Can process both genome and transcriptome inputs
  - Implements parallel processing with multi-threaded read analysis
  - Generates detailed strandedness reports with quantitative metrics
  - Supports SAM/BAM and FASTQ input formats
**Usage**: Analyzes RNA-seq data to determine library preparation strand orientation, critical for downstream differential expression and transcript quantification workflows

#### CallPeaks (CallPeaks.java)
**Purpose**: Peak calling algorithm for k-mer coverage histograms in genomic sequencing data
**Core Function**: Identifies significant peaks in k-mer coverage distribution to estimate genome characteristics like ploidy, genome size, and repeat content
**Key Features**:
  - Supports multiple peak detection modes (by volume or height)
  - Calculates genome size, ploidy, heterozygosity rate
  - Handles progressive smoothing of coverage data
  - Supports GC content analysis in peaks
  - Configurable peak filtering parameters
**Usage**: Analyzes sequencing data coverage histograms to provide insights into genome structure and composition, particularly useful for genome assembly and characterization

#### Consect (Consect.java)
**Purpose**: Generates a consensus from multiple error correction results for sequencing reads.
**Core Function**: Iteratively process batches of reads from multiple error correction tools, applying a consensus algorithm to correct base calls.
**Key Features**:
  - Supports multiple input error-corrected read files
  - Validates read identities and consistency
  - Tracks detailed correction and disagreement statistics
  - Handles quality score adjustments
  - Supports multithreaded read processing
**Usage**: Used in genomic sequencing workflows to refine read quality by comparing multiple error correction results

#### CountDuplicates (CountDuplicates.java)
**Purpose**: Probabilistically count and remove duplicate sequencing reads with minimal memory overhead
**Core Function**: Uses hashcode-based sampling to track read duplicates with configurable hash parameters and sampling rates
**Key Features**:
  - Probabilistic duplicate detection using 20 bytes per unique read
  - Configurable hash modes: bases, names, qualities
  - Optional separate output streams for duplicate and non-duplicate reads
  - Multi-threaded processing with concurrent read streams
  - Supports bit-mask sampling for large datasets
  - Generates detailed duplication statistics report
**Usage**: Preprocessing sequencing data to remove redundant reads, estimate duplication rates, and clean datasets for downstream bioinformatics analysis

#### CheckStrand (CheckStrand.java)
**Purpose**: Determines the strandedness of RNA-seq reads using k-mer sketching techniques
**Core Function**: Computes the strand bias of sequencing reads by comparing canonical and forward k-mer sketches to assess library strand specificity
**Key Features**:
  - Multi-threaded strandedness calculation with flexible sampling
  - Supports both normalized and depth-aware strandedness metrics
  - Adaptive k-mer sketching with configurable sketch size
  - Optional gene-based strand validation using reference transcriptome
  - Handles single-end and paired-end sequencing libraries
  - Detailed strand bias statistical reporting
**Usage**: Used in RNA-seq library preparation quality control to verify strand specificity and characterize sequencing library orientation

#### CountDuplicatesBuffered (CountDuplicatesBuffered.java)
**Purpose**: Probabilistically counts duplicate reads with minimal memory overhead
**Core Function**: Tracks read uniqueness using a configurable hashing strategy with thread-local buffering.
**Key Features**:
  - Uses 12-mer precision hash table for duplicate detection
  - Supports optional hashing of bases, names, and qualities
  - Configurable sampling rate for computational efficiency
  - Optional duplicate read output
  - Thread-safe with concurrent read/write processing
**Usage**: Identifies and optionally filters/reports duplicate reads in sequencing datasets

These are selected entries; this page is not a complete package inventory.

#### Reformat (reformat.sh)

Read-name uniqueness is disabled by default. Enable it with `uniquenames=t`:

```bash
bash ./reformat.sh in=reads.fq out=renamed.fq uniquenames=t trd=f ow=f
```

An unused name keeps its spelling. Collisions add `_2`, `_3`, and so on, avoiding
names already encountered, including names with natural numeric suffixes.
`uniquenamescope=identifier` is the default: it compares the identifier before
the first whitespace and preserves the following description. To compare the
complete retained header, add `uniquenamescope=header`. Selecting a scope does
not enable uniqueness by itself.

Set `trd=f` to retain descriptions. With `trd=t`, description trimming occurs
before uniqueness comparison in either scope. For these three input headers,
`<TAB>` represents one actual tab character:

| Input order | Identifier scope, `trd=f` | Header scope, `trd=f` | Either scope, `trd=t` |
|---|---|---|---|
| `x a` | `x a` | `x a` | `x` |
| `x<TAB>b` | `x_2<TAB>b` | `x<TAB>b` | `x_2` |
| `x a` | `x_3 a` | `x a_2` | `x_3` |

First and second mates have independent naming histories. When a terminal mate
annotation remains in the header, generated suffixes precede it: a repeated
`p/1` becomes `p_2/1`; retained `p /1` and `p 1:` become `p_2 /1` and `p_2 1:`.

#### CovStatsLine (CovStatsLine.java)

This mutable helper parses coverage-table rows after `initializeHeader` sets the
column map. The map is shared across the process and unsynchronized; callers must
coordinate schema changes with parsing. Constructors reject nonfinite scalar
values. `reads()` checks the strand-count sum, and `add()` checks base, strand and
read totals and rejects nonfinite weighted results before changing the receiver.
These summaries do not retain sufficient information to pool medians or standard
deviations reliably. Successful merges emit literal `NA` for both, retaining their
columns and positions; the parser accepts `NA` in `Median_fold` and `Std_Dev`.
Known singleton values stay numeric. Unavailable numeric fields are zero; check
`medianAvailable` or `stdDevAvailable`, or use `medianText()` or `stdDevText()`, to
distinguish unavailable values from known zeros. Missing columns retain the legacy
known-zero default.
Finite ranges, count relationships
and arbitrary public-field changes are not generally validated. See the class
Javadoc for field units, serialization order and the long-field compatibility note.

## Command Line Usage

Use the selected tool's shell launcher from the BBTools installation directory.
For example, filter an assembly using a coverage table with matching contig IDs:

```bash
bash ./filterbycoverage.sh in=assembly.fa cov=coverage.tsv out=clean.fa \
    outdirty=removed.fa mincov=5 minpercent=40 ow=f
```

This example writes passing contigs to `clean.fa` and filtered contigs to
`removed.fa`, with minimum mean coverage 5 and minimum covered percent 40.
Contigs missing from the coverage table are filtered. Other thresholds retain
their defaults; run `bash ./filterbycoverage.sh --help` for the available options.
The launcher supplies the Java environment, memory settings and assertions.

Common parameters across utilities:
- `in=file`: Input file (FASTQ, FASTA, SAM, etc.)
- `out=file`: Output file
- `overwrite=t`: Overwrite existing output files
- `threads=N`: Number of processing threads
- `verbose=t`: Enable detailed logging

## Integration with BBTools Ecosystem

The JGI package utilities are designed to work seamlessly with other BBTools components:
- Input/output compatibility across the suite
- Consistent parameter naming conventions
- Shared underlying libraries for file I/O and processing
- Optimized for high-performance genomic data processing

## Performance and Scalability

Threading, streaming and compressed-format support vary by tool. Consult that
tool's launcher help for its supported inputs, worker-thread options and heap
settings. Choose resources for the actual tool and workload, and measure that
configuration when comparing performance.
