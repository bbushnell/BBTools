# PacBio Package

This package contains specialized tools for processing and analyzing Pacific Biosciences (PacBio) long-read sequencing data. The tools cover adapter removal, read partitioning, contig manipulation, genomic site analysis, coverage calculation, and synthetic data generation.

## Classes

## CalcCoverageFromSites
**Purpose**: Calculate coverage statistics from alignment site information
**Core Function**: Processes genomic alignment sites to compute detailed coverage metrics across chromosomes
**Key Features**:
- Calculates total and correct coverage for genomic bases
- Tracks coverage statistics for defined and N bases
- Supports minimum coverage threshold filtering
- Generates comprehensive coverage reports with percentage breakdowns
**Usage**: Analyze alignment coverage for genomic sequencing experiments, providing insights into alignment accuracy and depth

## GenerateMultiChrom
**Purpose**: Generates synthetic multi-chromosome genome sequences with configurable mutations
**Core Function**: Creates multiple chromosome copies with controlled variations and optional N-region insertions
**Key Features**:
- Generates synthetic chromosomes with configurable substitution and indel rates
- Adds random N-regions to simulate assembly challenges
- Supports multiple chromosome copy generation
- Allows custom genome build and input source specification
**Usage**: Used for creating synthetic genome datasets with controlled genetic variations

## MakePacBioScript
**Purpose**: Generates customized workflow scripts for processing PacBio reads
**Core Function**: Parses command-line arguments to create shell scripts for PacBio read error correction and processing
**Key Features**:
- Supports multiple input file types (raw, gzipped, bzipped)
- Configurable parameters for threads, memory, and genome size
- Template-based script generation for different workflow modes
**Usage**: Automates script creation for PacBio read preprocessing and error correction workflows

## MergeFastaContigs
**Purpose**: Main utility for merging FASTA and FASTQ contigs into artificial chromosomes
**Core Function**: Processes input sequence files, concatenates contigs with N-padding, and generates structured FASTA output
**Key Features**:
- Supports merging multiple FASTA or FASTQ input files
- Adds configurable N-padding between contigs
- Generates chromosome-level output with position tracking
- Handles contig size and count limits
**Usage**: Combining genomic contigs into longer, more manageable sequences for genomic analysis

## MergeReadsAndGenome
**Purpose**: Main entry point for merging reads and genome with extensive configuration
**Core Function**: Processes command-line arguments to generate and merge synthetic reads with genome data
**Key Features**:
- Supports synthetic read generation with configurable parameters
- Handles multiple input read files
- Applies quality corrections to reads
- Assigns sequential numeric IDs to reads
**Usage**: Used for preprocessing and merging genomic read datasets in bioinformatics pipelines

## PartitionFastaFile
**Purpose**: Utility for splitting large FASTA files into manageable partitions
**Core Function**: Reads input FASTA file and creates multiple output files based on specified base pair count
**Key Features**:
- Preserves complete sequence integrity during splitting
- Dynamically creates new output files when partition size is reached
- Tracks total bases written and provides processing statistics
- Supports configurable maximum output length
**Usage**: Used for breaking down large genomic sequence files into smaller, processable chunks

## PartitionReads
**Purpose**: Main entry point for read partitioning with configurable parameters
**Core Function**: Distributes input reads across multiple output files/partitions
**Key Features**:
- Supports multiple input file formats (FASTQ, FASTA)
- Configurable number of output partitions
- Handles paired-end and single-end read datasets
- Flexible output file naming with partition indexing
**Usage**: Used for splitting large read datasets into smaller, manageable partitions for parallel processing or analysis

## ProcessStackedSitesNormalized
**Purpose**: Processes and normalizes stacked genomic sites with advanced scoring and filtering
**Core Function**: Processes genomic site scores across intervals, applying normalization and retention strategies
**Key Features**:
- Normalizes site scores based on read length and position
- Filters sites using configurable thresholds
- Supports weighted retention of top-scoring sites
- Handles site scoring across genomic intervals
**Usage**: Used in genomic analysis to process and filter site scores from sequencing data

## RemoveAdapters2
**Purpose**: Specialized tool for detecting and removing PacBio sequencing adapters from long-read data
**Core Function**: Multi-threaded adapter detection and removal pipeline with high-sensitivity alignment
**Key Features**:
- Supports forward and reverse complement adapter searches
- Configurable alignment score thresholds for adapter detection
- Optional read splitting at adapter sites
- Locality-aware adapter detection algorithm
**Usage**: Preprocessing PacBio sequencing reads to clean adapter contamination before downstream genomic analysis

## RemoveAdapters3
**Purpose**: Increased sensitivity to nearby adapters
**Core Function**: Multi-threaded PacBio read adapter removal with reverse-complementation analysis
**Key Features**:
- Processes reads concurrently with multiple threads
- Detects adapters in both forward and reverse complement orientations
- Supports splitting reads at adapter locations
- Configurable alignment score thresholds
**Usage**: Preprocessing PacBio sequencing reads to remove adapter contamination

## SiteR
**Purpose**: Compact bit-packed representation of genomic alignment site with efficient encoding of read and chromosome information
**Core Function**: Stores alignment site details using bit-level optimization to minimize memory usage
**Key Features**:
- Bit-packs chromosome, strand, read ID, and pair number
- Supports linked list chaining of related sites
- Efficient extraction of encoded site metadata
- Provides comparison and text representation methods
**Usage**: Used in genomic alignment tracking and processing in PacBio sequencing data analysis

## SortSites
**Purpose**: Sorts genomic sites across temporary files with configurable processing modes
**Core Function**: Distributes and sorts site data using position or ID-based keys, managing temporary file storage and merging
**Key Features**:
- Supports sorting by genomic position or read ID
- Filters and manages perfect/semiperfect sites
- Dynamically creates and manages temporary file storage
- Configurable block size and interval settings
**Usage**: Used in genomic data processing to organize and sort site information efficiently, particularly in variant analysis and read mapping pipelines

## SplitOffPerfectContigs
**Purpose**: Main entry point for perfect contig separation with extensive genome configuration
**Core Function**: Processes command line arguments to separate and filter contigs based on coverage and quality metrics
**Key Features**:
- Supports genome build input and output configuration
- Performs coverage-based contig filtering
- Enables splitting of low-quality contigs
- Generates detailed processing statistics
**Usage**: Separates high-quality contigs from reference genome sequences using coverage thresholds

## StackSites
**Purpose**: Processes and aggregates alignment sites for genomic read mapping
**Core Function**: Reads input files of alignment sites, processes site quality, and outputs stacked site data with coverage arrays
**Key Features**:
- Supports paired-end and single-end read processing
- Validates alignment perfection and correctness
- Generates chromosome-specific coverage arrays
- Sorts and writes site scores with interval-based organization
**Usage**: Used in genomic read alignment and mapping pipelines to analyze and organize alignment sites across chromosomes

## StackSites2
**Purpose**: Process and filter genomic site alignments across multiple chromosomes with distributed temporary file management
**Core Function**: Analyzes read alignments, calculates coverage, and selectively retains genomic sites based on alignment quality
**Key Features**:
- Processes paired-end reads with comprehensive site scoring
- Implements distributed temporary file management for large genomic datasets
- Supports genomic coverage tracking across chromosomes
- Filters sites based on alignment perfection and coverage thresholds
**Usage**: Used for advanced genomic alignment and site filtering in bioinformatics research, particularly with PacBio sequencing data