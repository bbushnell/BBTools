# JGI Package - BBTools Comprehensive Utilities Collection

The JGI package contains 92 specialized utilities designed for genomic sequencing data analysis, read processing, quality control, and bioinformatics workflows. These tools provide comprehensive functionality for modern genomic research and data processing pipelines.

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
**Purpose**: Generates statistical metrics for assembly files using the Assembly class.
**Core Function**: Processes multiple assembly files, extracting key metrics like total length, contig count, GC content, and contig length distributions.
**Key Features**:
  - Processes multiple input assembly files
  - Calculates assembly size and contig statistics
  - Computes GC content
  - Generates summary report with metrics
  - Supports custom output file configuration
**Usage**: Generates assembly quality report for genomic assemblies, useful in genome research and bioinformatics pipeline quality control.

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
**Purpose**: Programmatically add and simulate adapter sequences to sequencing reads
**Core Function**: Randomly inject adapter sequences into reads with configurable location, probability, and error simulation
**Key Features**:
  - Supports both 5' and 3' adapter injection with configurable probabilistic placement
  - Optional reverse-complement adapter generation
  - Quality-aware error simulation during adapter insertion
  - Flexible read processing for single and paired-end sequencing data
  - Tracks detailed statistics about adapter addition and read modifications
**Usage**: Used in sequencing data preprocessing to simulate adapter contamination, evaluate adapter removal algorithms, and generate synthetic training datasets for bioinformatics tools

#### AdjustHomopolymers.java (AdjustHomopolymers.java)
**Purpose**: Command-line utility for processing sequencing reads by adjusting homopolymer runs
**Core Function**: Modifies read bases and qualities by expanding or contracting homopolymer sequences based on a configurable rate
**Key Features**:
  - Supports single and paired-end read processing
  - Configurable homopolymer adjustment rate
  - Handles read validation and quality tracking
  - Preserves read quality scores during adjustments
  - Supports various input/output file formats (FASTQ)
**Usage**: Used in bioinformatics pipelines to normalize sequence reads by standardizing homopolymer regions

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

*[Remaining utilities follow the same detailed documentation pattern.]*

## Command Line Usage

Most JGI utilities follow a standard BBTools command-line interface:

```bash
java -cp bbmap.jar jgi.UtilityName [parameters]
```

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

All utilities in this package are optimized for:
- Multi-threaded processing capabilities
- Memory-efficient algorithms
- Support for compressed file formats
- Streaming data processing where applicable
- Large-scale genomic dataset handling

This comprehensive utility collection serves as the foundation for complex bioinformatics workflows, providing researchers with robust, efficient tools for genomic data analysis and processing.