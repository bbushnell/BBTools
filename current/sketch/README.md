# Sketch Package - MinHash Genomic Sketching and Similarity Analysis

BBTools sketch package providing comprehensive MinHash sketching implementation for rapid genomic sequence comparison, taxonomic classification, and large-scale similarity analysis.

## Core Sketching Infrastructure

### Sketch (Sketch.java)
**Purpose**: Represents a computational sketch of a genomic sequence using MinHash algorithm for efficient sequence comparison
- **Core Function**: Manages k-mer hashes, taxonomic metadata, and comparison capabilities for genomic sequence analysis
- **Key Features**:
  - Stores sorted hashcodes with specialized encoding (Long.MAX_VALUE-(raw hashcode))
  - Supports k-mer count tracking with optional frequency information
  - Handles base composition tracking (ACGT counts)
  - Includes taxonomic metadata like taxID, genome size, and sequence information
  - Supports 16S and 18S ribosomal RNA sequence attachment
  - Implements set operations like intersection and union of sketches
- **Usage**: Genomic sequence comparison and similarity analysis, taxonomic classification and identification, efficient representation of large genomic datasets

### SketchMaker (SketchMaker.java)
**Purpose**: Creates MinHashSketches rapidly for genomic sequence analysis
- **Core Function**: Generates MinHash sketches from input sequences with flexible processing modes
- **Key Features**:
  - Multiple sketching modes (per-taxa, per-sequence, single sketch)
  - Multithreaded processing with concurrent read input streams
  - Configurable taxonomic and IMG (Integrated Microbial Genomes) filtering
  - Supports paired-end read processing
  - Flexible output and metadata generation
- **Usage**: Generating compact genomic signatures for comparison, classification, and clustering of biological sequences
- **Processing Modes**:
  - **ONE_SKETCH**: Creates a single sketch from all input sequences
  - **PER_TAXA**: Generates separate sketches for each taxonomic group
  - **PER_IMG**: Creates sketches for each Integrated Microbial Genome entry
  - **PER_SEQUENCE**: Generates individual sketches for each input sequence

### SketchMakerMini (SketchMakerMini.java)
**Purpose**: Creates MinHashSketches rapidly for genomic sequence analysis.
- **Core Function**: Generates computational sketches from genomic reads using advanced k-mer processing techniques, supporting multiple input modes and sequence types
- **Key Features**:
  - Supports processing of nucleotide, amino acid, and translated sequences
  - Configurable MinHash sketching with flexible entropy and quality filtering
  - Multi-mode processing: one sketch per file, per-sample, or multi-sketch
  - Advanced k-mer tracking with both forward and reverse complement handling
  - Integrated support for 16S rRNA sequence detection
- **Usage**: Generate genomic sketches for sequence comparison and taxonomic classification, process reads from various sequencing technologies

### SketchTool (SketchTool.java)
**Purpose**: Converts k-mer tables into sketches for genomic sequence comparison
- **Core Function**: Generates compact representative summaries of genomic sequences using k-mer sampling techniques
- **Key Features**:
  - Multithreaded sketch generation from k-mer tables
  - Configurable sketch size and key occurrence thresholds
  - Supports different sketch encoding modes (A48, HEX, NUC)
  - Flexible loading of sketches from files and sequence data
  - Handles paired-end read processing
- **Usage**: Genomic sequence comparison, metagenomics, and taxonomic classification through k-mer signature generation

## Sketch Comparison and Analysis

### CompareSketch (CompareSketch.java)
**Purpose**: Compares one or more input sketches to a set of reference sketches.
- **Core Function**: A command-line tool for processing and comparing computational sketches, supporting various input and output configurations for genomic data analysis.
- **Key Features**:
  - Multi-threaded sketch comparison 
  - Flexible input/output file handling
  - Supports various comparison modes (per-file, per-taxa)
  - Customizable metadata overriding
  - JSON and text output formats
  - Parallel compression/decompression support
- **Usage**: Used for genomic sequence comparison, taxonomic classification, and computational sketch analysis across multiple input files

### Comparison (Comparison.java)
**Purpose**: Manages comparative metrics between two genomic sketches for taxonomic and sequence similarity analysis
- **Core Function**: Calculates comprehensive similarity metrics between two genetic sketches using multi-dimensional computational approaches
- **Key Features**:
  - Multi-stage k-mer hit counting and analysis
  - Taxonomy identification through k-mer comparison
  - Advanced score calculation using depth, hits, and contamination metrics
  - Flexible comparison using multiple comparator strategies (score, depth, volume, hits)
  - Integrated support for 16S and 18S rRNA sequence identity calculation
  - Supports dual k-mer ANI (Average Nucleotide Identity) estimation
- **Usage**: Provides detailed comparative analysis for genomic sketches, enabling taxonomic classification, sequence similarity assessment, quality evaluation

### SketchSearcher (SketchSearcher.java)
**Purpose**: Manages sketch-based sequence comparison and reference database searches for genomic and protein data.
- **Core Function**: Performs computational comparisons between genetic sketches using multi-threaded techniques, supporting various reference databases and taxonomic filtering.
- **Key Features**:
  - Multi-threaded sketch comparison across reference databases
  - Supports multiple reference datasets (NT, RefSeq, SILVA, IMG, Protein)
  - Flexible taxonomic filtering and level-based search
  - Concurrent sketch index management
  - Customizable comparison parameters (minimum hits, minimum ANI, etc.)
- **Usage**: Used in bioinformatics for rapid sequence similarity searches, particularly for taxonomic classification of genomic sequences, comparative genomics, protein sequence matching

### SendSketch (SendSketch.java)
**Purpose**: Compares input sketches to a set of reference sketches with customizable processing and output configurations.
- **Core Function**: Command-line tool for processing and sending computational genomic sketches to remote servers for comparison, supporting multi-mode analysis with complex input handling
- **Key Features**:
  - Multi-threaded sketch processing and comparison
  - Support for local and remote sketch sending modes
  - Flexible address resolution for multiple genomic databases (NT, RefSeq, Silva, IMG, etc.)
  - Customizable metadata overriding capabilities
  - Supports JSON and text output formats
  - Configurable blacklist and taxonomy filtering
- **Usage**: Send genomic sketches to remote servers for large-scale taxonomic comparisons, process and filter input sketches with complex parameter configurations

## Data Structures and Management

### SketchHeap (SketchHeap.java)
**Purpose**: Specialized heap data structure for managing genome sketch information and associated metadata.
- **Core Function**: Manages a collection of long-valued keys with support for set and map-based storage modes, optimized for genomic sketching operations.
- **Key Features**:
  - Supports flexible storage modes (set or map) based on key occurrence count
  - Tracks genomic metadata like genome size, sequences, and taxonomic information
  - Provides methods for adding, checking, and converting sketch data
  - Implements genome size estimation algorithms
  - Supports blacklist and whitelist filtering for keys
- **Usage**: Used in genomic sketching and comparison processes, helps estimate genome characteristics and manage hash-based representations

### SketchIndex (SketchIndex.java)
**Purpose**: High-performance indexing mechanism for sketch-based data processing
- **Core Function**: Creates and manages an efficient multi-threaded index for sketches, enabling rapid similarity searches across large datasets
- **Key Features**:
  - Concurrent indexing via IndexThread inner class
  - Supports multiple lookup strategies (list and map-based)
  - Handles taxonomy-level filtering
  - Configurable thread usage and index limits
  - Automatic allocation of k-mer hash tables
- **Usage**: Used in genomic sketch comparison to quickly find similar sketches across reference datasets

### SketchObject (SketchObject.java)
**Purpose**: Handles complex k-mer hashing and transformation for genomic sequence sketching
- **Core Function**: Provides advanced algorithmic transformations for generating genomic sequence signatures
- **Key Features**:
  - Multiple k-mer hashing strategies (hash1, hash2, hash3)
  - Support for reverse complement and amino acid sequences
  - Dynamic parameter configuration for k-mer size and hash functions
  - Genomic size estimation and sketch size calculation
  - Specialized hashing for different genomic contexts
- **Usage**: Core computational engine for generating compact genomic representations used in similarity comparisons and taxonomic classification

## Result Processing and Analysis

### SketchResults (SketchResults.java)
**Purpose**: Manages and processes results from genomic sketch comparisons, providing comprehensive result tracking and analysis capabilities.
- **Core Function**: Stores, filters, and processes comparison results between genomic sketches, supporting complex taxonomic and sequence similarity evaluations
- **Key Features**:
  - Supports concurrent result collection via `ConcurrentHashMap`
  - Dynamic result list management with flexible sorting
  - Advanced recomparison and filtering mechanisms
  - Multi-level taxonomic hit tracking
  - Configurable result display and output formatting
- **Usage**: Used in genomic sketch comparison workflows to collect and manage comparison results across multiple sketches, apply complex filtering based on taxonomic levels

### AnalyzeSketchResults (AnalyzeSketchResults.java)
**Purpose**: Command-line tool for analyzing genomic sketch comparison results across multiple taxonomic levels
- **Core Function**: Processes and validates sketch comparison results from various genomic comparison methods
- **Key Features**:
  - Supports multiple sketch comparison modes (MASH, Sourmash, BLAST, BBSketch)
  - Comprehensive result parsing and statistical analysis
  - Configurable taxonomic level reporting
  - Multi-threaded accuracy processing
  - Detailed error handling and logging
- **Usage**: Validate and analyze genomic similarity measurements, including Average Nucleotide Identity (ANI) and Alignment Accuracy Index (AAI)

### SummarizeSketchStats (SummarizeSketchStats.java)
**Purpose**: Command-line tool for summarizing and processing genomic sketch comparison results with advanced taxonomic filtering.
- **Core Function**: Processes and analyzes genomic sketch results files, extracting detailed comparative metrics across taxonomic levels
- **Key Features**:
  - Supports flexible input processing from multiple sketch result files
  - Configurable taxonomic level filtering with `taxLevel` parameter
  - Handles various result parsing modes (BBSketch, MASH)
  - Supports detailed output customization with header and formatting options
  - Can filter results based on taxonomic similarity and unique hits
- **Usage**: Summarize and validate genomic sketch comparison results, extract comprehensive statistical metrics from genomic comparisons

### ResultLineParser (ResultLineParser.java)
**Purpose**: Parses tabular result lines from genomic comparison sketches, extracting taxonomic and similarity metrics.
- **Core Function**: Processes raw result lines from different sketch comparison modes (BBSketch and MASH), extracting and storing taxonomic, sequence, and similarity information
- **Key Features**:
  - Supports multi-mode parsing for BBSketch and MASH result formats
  - Extracts detailed taxonomic identifiers, sequence sizes, and comparison metrics
  - Handles complex parsing of ANI (Average Nucleotide Identity) and SSU (Small Subunit) data
  - Supports optional text preservation for result tracking
- **Usage**: Processes genomic sketch comparison results from different computational pipelines, extracts and normalizes taxonomic and similarity information

## Filtering and Quality Control

### Blacklist (Blacklist.java)
**Purpose**: Manages centralized k-mer blacklisting for sequence filtering and identification across multiple genomic databases.
- **Core Function**: Provides a flexible system for loading, tracking, and querying k-mer blacklists from various genomic reference databases
- **Key Features**:
  - Supports multiple database-specific blacklists (NT, RefSeq, Silva, IMG, etc.)
  - Dynamic runtime blacklist loading and resolution
  - Thread-safe k-mer tracking with hash table storage
  - Canonical k-mer representation using complementary storage
  - Efficient memory allocation for large-scale k-mer sets
- **Usage**: Used in genomic sequence filtering to exclude known repetitive, low-complexity, or undesired sequence k-mers across bioinformatics preprocessing workflows

### BlacklistMaker (BlacklistMaker.java)
**Purpose**: Creates taxa-aware blacklists from genomic reads by processing k-mer distributions.
- **Core Function**: Generates a filtered list of k-mers based on taxonomic occurrence and processing parameters, useful for identifying and filtering out common or low-quality genomic elements.
- **Key Features**:
  - Multiple processing modes: per-sequence, per-taxa, and per-IMG (Integrated Microbial Genomes)
  - Configurable k-mer occurrence thresholds via `minTaxCount`
  - Multi-threaded, parallel processing for efficient large-scale genomic analysis
  - Adaptive memory allocation for prefiltering and k-mer counting
  - Supports both nucleotide and amino acid k-mer processing
- **Usage**: Generate genomic blacklists for filtering sequence sketches, identify and remove k-mers that are taxonomically uninformative

### BlacklistMaker2 (BlacklistMaker2.java)
**Purpose**: Creates blacklists from existing sketches by processing k-mer distributions.
- **Core Function**: Generates a filtered list of k-mers based on taxonomic or sequence occurrence thresholds, useful for filtering out common or low-quality genomic elements.
- **Key Features**:
  - Supports multiple processing modes: per-taxa, per-sequence, or per-IMG (Integrated Microbial Genomes)
  - Configurable k-mer occurrence thresholds with `minTaxCount` parameter
  - Multi-threaded processing for efficient large-scale sketch analysis
  - Ability to filter out unranked or low-quality taxonomic assignments
  - Flexible taxonomic level promotion during processing
- **Usage**: Generate genomic blacklists for filtering sequence sketches, identify and remove k-mers that are too common or taxonomically uninformative

### Whitelist (Whitelist.java)
**Purpose**: Manages a set of k-mer values for filtering genomic sketches using a specialized hashing mechanism
- **Core Function**: Provides a method to check if a k-mer is present in a pre-defined collection of k-mer tables, supporting sketch filtering
- **Key Features**:
  - Uses a hash-based method to distribute k-mers across 31 different k-mer tables
  - Supports initialization of k-mer tables for whitelist checking
  - Implements methods to check k-mer presence with both raw and hashed values
  - Provides a flag to check if a whitelist has been initialized
- **Usage**: Filter genomic sketches by keeping only whitelisted k-mers, enables selective retention of k-mers based on predefined criteria

## Specialized Processing Tools

### CompareBuffer (CompareBuffer.java)
**Purpose**: Tracks and calculates various hit and comparison metrics for k-mer sequence analysis.
- **Core Function**: Manages computational buffer for tracking sequence hits, depth, and similarity metrics across different k-mer lengths
- **Key Features**:
  - Tracks multiple hit types: unique, multi-hit, contamination hits
  - Supports dual k-mer length analysis with advanced ANI (Average Nucleotide Identity) calculations
  - Provides methods to compute depth, reference hit averages, and genome similarity
  - Optional BitSet allocation for k-mer tracking with memory efficiency
- **Usage**: Used in genomic sequence comparison and hit analysis, calculates genome similarity metrics across different k-mer sizes

### MergeSketch (MergeSketch.java)
**Purpose**: Combines multiple genomic sketches into a single, representative sketch using multithreaded processing.
- **Core Function**: Loads multiple sketch files, merges their k-mer information using a heap-based algorithm, and generates a unified sketch output
- **Key Features**:
  - Multithreaded sketch loading for performance
  - Configurable sketch size with automatic genome-wide representation
  - Flexible metadata preservation and override
  - Supports various command-line metadata customizations
  - Preserves total k-mer information across input sketches
- **Usage**: Consolidate multiple genomic sketches into a single representative sketch, merge sketches from different samples or genomic sources

### SubSketch (SubSketch.java)
**Purpose**: Generates smaller sketches from input sketches with configurable processing parameters.
- **Core Function**: Reads input sketches, processes and resizes them based on target sketch size, and supports multiple output configuration modes
- **Key Features**:
  - Supports multi-file input processing with configurable output destination
  - Implements target sketch size calculation via `toSketchSize()` method
  - Handles blacklist filtering with `applyBlacklist()` method
  - Supports resizing of sketches larger than target size
  - Provides verbose mode for detailed processing information
- **Usage**: Reduce computational complexity by generating smaller, representative genomic sketches, filter out k-mers using blacklists for improved computational efficiency

## K-mer Processing and Utilities

### KmerLimit (KmerLimit.java)
**Purpose**: Manages k-mer count limit processing for genomic sequence data
- **Core Function**: Processes reads to generate and filter k-mers based on specified parameters, supporting early termination when target k-mer count is reached
- **Key Features**:
  - Configurable k-mer length (up to 32 bases)
  - Quality-aware k-mer processing with probability thresholds
  - Multithreaded read processing
  - Supports paired-end and single-end sequencing
  - Dynamic heap-based k-mer tracking
- **Usage**: Genomic analysis tools requiring controlled k-mer sampling or early processing termination

### KmerLimit2 (KmerLimit2.java)
**Purpose**: Command-line tool for subsampling genomic sequences based on k-mer occurrence and frequency
- **Core Function**: Computes k-mer frequency histogram and intelligently subsamples reads to meet target k-mer diversity
- **Key Features**:
  - Monte Carlo sampling rate calculation for precise read subsampling
  - Multi-threaded k-mer counting and histogram generation
  - Configurable k-mer length (default 32)
  - Advanced quality filtering with probability and base quality thresholds
  - Supports paired-end and single-end read processing
- **Usage**: Reduce computational complexity of large genomic datasets, generate representative subset of reads maintaining k-mer diversity

### InvertKey (InvertKey.java)
**Purpose**: Searches input sequences for specific k-mer hashes and returns their canonical representation
- **Core Function**: Performs targeted k-mer search in genomic sequences, finding occurrences of specified k-mer hashes and outputting their canonical sequence and location
- **Key Features**:
  - Supports single and multiple k-mer target searches
  - Handles both forward and reverse complement k-mer representations
  - Configurable k-mer length (primary and secondary)
  - Multiple input and output format support (FASTA, text)
  - Optional early termination mode
- **Usage**: Targeted k-mer location tracking in genomic sequences, finding specific k-mer occurrences in read datasets

## SSU (Small Subunit) Processing

### AddSSU (AddSSU.java)
**Purpose**: Adds or processes 16S and 18S ribosomal RNA sequences to genomic sketches.
- **Core Function**: Manages addition, filtering, and manipulation of SSU (Small Subunit) sequences for taxonomic sketches
- **Key Features**:
  - Supports loading 16S and 18S reference databases
  - Provides flexible taxonomic sequence filtering options
  - Handles eukaryotic and prokaryotic sequence processing
  - Supports selective SSU sequence addition or clearing
- **Usage**: Used for preparing and refining genomic sketches with specific ribosomal RNA sequence requirements, enabling detailed taxonomic classification and analysis

### CompareSSU (CompareSSU.java)
**Purpose**: Compares Sequence Similarity Units (SSUs) across multiple taxonomic levels with parallel processing.
- **Core Function**: Performs all-to-all or fractional matrix comparisons of SSU sequences with multi-threaded computational efficiency
- **Key Features**:
  - Supports multi-threaded SSU sequence comparison
  - Generates identity percentage metrics across taxonomic levels
  - Filters SSU sequences by length and N-base content
  - Can perform all-to-all or selective level comparisons
  - Utilizes tax tree for taxonomic level identification
- **Usage**: Used for comparative genomic analysis, specifically evaluating sequence similarity across different taxonomic classifications

### SSUMap (SSUMap.java)
**Purpose**: Manages loading and processing of 16S and 18S ribosomal RNA sequence maps for taxonomic reference.
- **Core Function**: Loads and maintains taxonomically-linked sequence mappings for ribosomal RNA reference sequences with efficient memory management and load-on-demand strategies.
- **Key Features**:
  - Supports loading of both 16S and 18S ribosomal RNA sequence maps
  - Provides synchronized map loading with thread-safe mechanisms
  - Handles automatic file selection with "auto" configuration option
  - Filters and selects longest representative sequences per taxonomic ID
  - Supports verbose logging and error state tracking
- **Usage**: Used in taxonomic classification and reference sequence management for genomic and microbial studies

## Data Management and Records

### Record (Record.java)
**Purpose**: Represents a genomic record with taxonomic and sequence comparison metadata.
- **Core Function**: Manages comparison data for taxonomic records, providing methods for SSU (Small Subunit) processing and comparison
- **Key Features**:
  - Implements `Cloneable` and `Comparable` interfaces for flexible record manipulation
  - Stores key taxonomic identifiers: query and reference Tax IDs
  - Tracks sequence bases, sizes, and sequence similarity metrics
  - Supports SSU (Small Subunit) sequence alignment and processing
  - Provides comparison mechanism based on Average Nucleotide Identity (ANI)
- **Usage**: Used in genomic sketch-based taxonomic comparisons, supports tracking and processing of taxonomic record metadata

### RecordSet (RecordSet.java)
**Purpose**: Manages and processes a collection of genomic Records with taxonomic level tracking and filtering.
- **Core Function**: Provides a container for genomic records with sorting, deduplication, and taxonomic level processing capabilities.
- **Key Features**:
  - Supports sorting and deduplicating records based on taxonomic levels
  - Implements a bit-mask based level tracking mechanism
  - Provides methods for SSU (Small Subunit) processing
  - Supports validation and testing of record taxonomic consistency
  - Allows adding and checking taxonomic levels
- **Usage**: Used in genomic sketch analysis for managing and processing sets of taxonomic records

## Parallel Processing Infrastructure

### AlignmentJob (AlignmentJob.java)
**Purpose**: Manages concurrent sequence alignment job processing with error resilience.
- **Core Function**: Executes sequence similarity calculations within a thread-safe, blocking queue-based concurrent processing framework
- **Key Features**:
  - Supports concurrent alignment job execution
  - Implements blocking queue mechanism for result delivery
  - Provides poison pill detection to prevent invalid job processing
  - Implements robust error handling with retry mechanisms
- **Usage**: Used in parallel processing scenarios for genomic sequence alignment tasks

### AlignmentThreadPool (AlignmentThreadPool.java)
**Purpose**: Manages thread pool for concurrent alignment comparisons with dynamic thread scaling.
- **Core Function**: Creates and manages a scalable thread pool for executing alignment jobs using a blocking work queue
- **Key Features**:
  - Dynamically spawns threads up to a maximum thread limit
  - Supports lazy thread creation based on workload
  - Implements a thread-safe, synchronous batch processing mechanism
  - Provides a poison pill shutdown mechanism for graceful thread termination
  - Uses atomic integer for tracking active jobs
- **Usage**: Handles parallel alignment comparisons by distributing work across multiple threads

## Utility Components

### DisplayParams (DisplayParams.java)
**Purpose**: Configuration and parameter management class for sketch output and comparison processing
- **Core Function**: Handles parsing, storing, and transforming display parameters for genomic sketch comparisons, supporting multiple output formats and filtering options
- **Key Features**:
  - Flexible parameter parsing from double-header lines
  - Supports multiple output formats (JSON, D3, etc.)
  - Comprehensive taxonomic filtering capabilities
  - Configurable result display options
  - Supports JSON and customizable result generation
- **Usage**: Configures display parameters for genomic sketch comparisons, controls what metadata and metrics are displayed in comparison results

### SketchIdComparator (SketchIdComparator.java)
**Purpose**: Provides a custom Comparator for sorting Sketch objects based on their unique sketch ID.
- **Core Function**: Compares two Sketch objects by their sketchID attribute, enabling sorted collection and retrieval of Sketch instances
- **Key Features**:
  - Implements java.util.Comparator<Sketch> interface
  - Provides a static final comparator instance for reuse
  - Uses simple integer subtraction for comparison
  - Singleton-like design with private constructor
- **Usage**: Used in sorting and ordering operations involving Sketch objects

## Package Overview

The sketch package provides a comprehensive MinHash sketching framework for:

1. **Genomic Sequence Analysis**: Rapid comparison of large genomic datasets using compact sketch representations
2. **Taxonomic Classification**: Efficient identification and classification of organisms based on genomic similarity
3. **Large-Scale Bioinformatics**: Scalable processing of genomic data with multi-threaded architecture
4. **Quality Control**: Advanced filtering mechanisms using blacklists, whitelists, and quality metrics
5. **Remote Processing**: Integration with remote servers for large-scale comparative genomics

The package supports multiple processing modes, from single-sequence analysis to large-scale metagenomics studies, with sophisticated parallel processing capabilities and comprehensive result analysis tools.