# Clump Package

The clump package provides tools for k-mer based sequence clumping, grouping similar reads together to improve compression efficiency and enable sophisticated sequence analysis and sorting operations.

## Classes

#### Clump (Clump.java)
**Purpose**: A specialized data structure for managing groups of reads that share a common k-mer hash.

- **Core Function**: Performs advanced read clustering, duplicate removal, and error correction for sequencing data
- **Key Features**:
  - Collects reads sharing a specific k-mer hash into a single cluster
  - Generates consensus reads by aggregating base and quality information
  - Implements sophisticated duplicate removal strategies
  - Supports error correction through statistical base calling
  - Handles optical and sequence-based duplicate detection
  - Allows configurable duplicate removal thresholds
  - Provides flexible consensus generation methods
- **Usage**: Used in bioinformatics pipelines for read preprocessing, deduplication, and quality improvement of sequencing data, particularly in DNA/RNA sequencing analysis

#### Clumpify (Clumpify.java)
**Purpose**: Advanced read clustering and sorting tool for sequencing data files, particularly FASTQ and FASTA formats.

**Core Function**: Reorganizes sequencing reads into optimized groups (clumps) using k-mer-based sorting and clustering algorithms to improve storage and processing efficiency.

**Key Features**:
- **Multi-Pass Clustering**: Supports multiple clustering passes to refine read organization
- **Automatic Memory Management**: Dynamically sets processing groups based on available system memory
- **Flexible Input Handling**: Supports single and paired-end sequencing files with various input formats
- **Advanced Sorting Options**: Provides name sorting, read repair, and quality quantization capabilities
- **Temporary File Management**: Intelligent temporary file creation with compression and naming options
- **Memory-Efficient Processing**: Uses k-mer comparison strategies to minimize memory consumption
- **Quality Control**: Optional quality-related transformations like ECCO and quality quantization

**Usage**: 
- Preprocessing sequencing data for more efficient storage and downstream analysis
- Reducing computational overhead in large sequencing datasets
- Preparing reads for alignment, assembly, or other bioinformatics workflows
- Optimizing read organization for improved I/O performance and reduced storage requirements

#### ClumpList (ClumpList.java)
**Purpose**: A specialized list for managing and processing groups of reads (clumps) with support for multi-threaded operations.

**Core Function**: 
- Manages collections of DNA read clumps, allowing dynamic addition, reordering, and processing of reads with k-mer based clustering
- Supports multi-threaded read processing and clump generation
- Provides methods for consensus read generation, error correction, and duplicate removal

**Key Features**:
- **Parallel Read Processing**: Can add reads using multi-threaded (`addReadsMT()`) or single-threaded (`addReadsST()`) methods
- **Clump Reordering**: Supports paired and unpaired read reordering with `reorderPaired()` and `reorder()` methods
- **Consensus Generation**: Ability to create consensus reads with optional simple consensus mode
- **Error Correction**: Supports read error correction and duplicate removal via `process()` method
- **K-mer Based Clustering**: Organizes reads into clumps based on k-mer similarity
- **Flexible Mode Operations**: Supports different processing modes like condensing, correcting, and deduplicating reads

**Usage**: 
- Used in DNA sequencing and read processing pipelines to organize, clean, and prepare sets of sequencing reads
- Particularly useful for preprocessing reads before further analysis, such as genome assembly or metagenomics studies
- Can handle large collections of reads with efficient multi-threaded processing
- Supports operations like generating consensus sequences, error correction, and removing duplicate reads

#### ClumpTools (ClumpTools.java)
**Purpose**: A utility class for managing k-mer count arrays and retrieving statistical tables for genomic read data.

**Core Function**: Provides static methods to create, retrieve, and manage a global k-mer count array using different input sources like read collections or file names.

**Key Features**:
- **Synchronized Table Retrieval**: Offers thread-safe methods to get k-mer count arrays from read collections or file names
- **Flexible Input Sources**: Supports creating k-mer count tables from:
  1. ArrayList of Read objects
  2. Input file names with configurable parameters
- **Dynamic Table Caching**: Maintains a static table that can be reused across method calls
- **Configurable k-mer Parameters**: Allows specifying k-mer size and minimum count threshold
- **Stream Management**: Automatically handles opening and closing input streams for read processing

**Usage**: Primarily used in genomic data analysis workflows for:
- Pre-processing read collections
- Generating k-mer frequency tables
- Supporting pivotal set computations in bioinformatics data processing

#### Hasher (Hasher.java)
**Purpose**: Provides robust hashing and comparison utilities for genetic sequence data, specifically designed for read-based comparisons and hash generation.

**Core Function**: Generates deterministic, low-collision hash codes for genetic sequences using a combination of length-based seeding and per-character randomized hash codes.

**Key Features**:
- **Deterministic Sequence Hashing**: Generates consistent long-integer hash values for byte-array genetic sequences
- **Multi-Mode Hash Generation**: Supports 32 different hashing modes using pre-computed randomized hash codes
- **Read Pair Hashing**: Provides specialized hashing for paired-end sequencing reads, combining individual read hashes
- **Case-Insensitive Character Handling**: Normalizes uppercase and lowercase characters during hash computation
- **Sequence Comparison**: Offers precise equality comparison for genetic sequences, handling null and empty reads
- **Rotation-Based Hash Mixing**: Uses bit rotation to enhance hash distribution and minimize collisions

**Usage**: 
- Used in bioinformatics tools for rapid sequence identification, deduplication, and clustering of genetic reads
- Supports comparing reads with different hashing strategies (single read, read pair)
- Critical for high-performance genomic data processing where quick sequence matching is required

#### KmerComparator (KmerComparator.java)
**Purpose**: A sophisticated k-mer-based read comparator for genomic sequence analysis with advanced hashing and comparison strategies.

**Core Function**: Provides a mechanism to compare genomic reads by generating canonical k-mer representations, handling various complexity scenarios like read orientation, quality, and length variations.

**Key Features**:
- **K-mer Canonicalization**: Generates canonical k-mers by selecting the lexicographically smaller representation between forward and reverse complement sequences
- **Adaptive Hashing**: Uses configurable multi-level hash functions (up to 8 hashes) for robust k-mer representation
- **Quality-Aware Comparison**: Incorporates base quality probabilities into k-mer selection and comparison process
- **Border Handling**: Supports configurable border regions to manage edge effects in sequence analysis
- **Paired-End Read Support**: Can merge and compare paired-end reads with intelligent overlap detection
- **Thread-Safe Design**: Implements thread-safe hashing through local thread keys and parallel processing mechanisms
- **Flexible Orientation Handling**: Option to reverse complement reads during comparison

**Usage**: Primarily used in genomic sequence clustering, deduplication, and comparative analysis workflows where precise read representation and comparison are critical. Particularly useful in:
- Metagenomic read classification
- Genome assembly preprocessing
- Sequence similarity clustering
- Error correction in sequencing data

#### KmerComparator2 (KmerComparator2.java)
**Purpose**: Lightweight comparative utility for k-mer (keyword-length DNA sequence) comparisons, providing a minimal implementation for sorting and comparing DNA reads.

- **Core Function**: Implements a flexible Comparator<Read> interface with advanced sequence comparison methods, allowing precise ordering of DNA reads based on multiple criteria
- **Key Features**:
  - Provides abstract comparison method for ReadKey objects, enabling custom comparison strategies
  - Supports multi-stage comparison including sequence, mate sequence, and read quality
  - Implements byte-level sequence comparison with null and length handling
  - Uses fallback to Read ID comparison for identical sequences
  - Offers static methods for comparing sequences and read qualities
  - Handles edge cases like null or different-length sequences
- **Usage**: Used in bioinformatics workflows for sorting, clustering, or organizing DNA reads based on sequence similarity, particularly in read preprocessing and clustering algorithms

#### KmerComparatorX (KmerComparatorX.java)
**Purpose**: Specialized k-mer comparator for advanced read key sorting in genomic sequencing data, extending base k-mer comparison logic with enhanced optical sorting capabilities.

- **Core Function**: Provides a custom comparison method for ReadKey objects, prioritizing k-mer size, strand orientation, and optional optical coordinates
- **Key Features**:
  - Sorts k-mers in descending order of k-mer value 
  - Handles both forward and reverse strand k-mer comparisons
  - Supports position-based secondary sorting
  - Implements optional optical coordinate sorting for sequencing data
  - Provides a static singleton comparator instance for efficient reuse
  - Allows conditional tile and lane sorting based on configuration flags
- **Usage**: Used in genomic sequencing data processing to deterministically order and group reads based on k-mer characteristics, particularly in optical mapping and read clustering workflows

#### KmerComparatorY (KmerComparatorY.java)
**Purpose**: Specialized k-mer comparator that provides a custom comparison strategy for read keys with a focus on Y-coordinate sorting for optical map generation.

**Core Function**: Implements a complex comparison method that prioritizes k-mer size, strand direction, position, and optionally Y-coordinate for sorting read keys in the clump package.

**Key Features**:
- **K-mer Size Prioritization**: Sorts read keys with larger k-mers first, enabling descending order comparison
- **Strand-Aware Comparison**: Differentiates between k-mers on positive and negative strands during sorting
- **Position-Based Secondary Sorting**: Uses read position as a secondary sorting criterion when k-mer and strand are identical
- **Optical Mapping Support**: Provides specialized sorting for optical mapping scenarios with optional lane and tile considerations
- **Y-Coordinate Refinement**: Offers granular sorting by Y-coordinate when optical-only mode is enabled
- **Static Comparator Instance**: Provides a pre-instantiated static comparator for efficient reuse
- **Extensible Design**: Extends KmerComparator2, allowing for inheritance and potential further customization

**Usage**: Used in bioinformatics workflows for sorting and organizing read keys during sequence analysis, particularly in optical mapping and k-mer-based clustering algorithms within the BBTools clump package.

#### KmerReduce (KmerReduce.java)
**Purpose**: Reduces reads to their feature k-mer, extracting a representative k-mer from input sequencing reads.

- **Core Function**: Processes DNA sequencing reads, converting entire reads into their most representative k-mer using a fixed-length k-mer hashing strategy.

- **Key Features**:
  - Supports concurrent read processing using multiple threads
  - Configurable k-mer length (default k=31, must be between 0-31)
  - Ability to perform optional read error correction (ECCO mode)
  - Flexible input handling for single and paired-end sequencing reads
  - Supports various file formats including FASTQ
  - Generates output with reduced read representation based on k-mer hashing
  - Handles complex read stream management with ConcurrentReadInputStream and ConcurrentReadOutputStream

- **Usage**: 
  - Primarily used in bioinformatics for reducing high-dimensional sequencing reads to compact, representative k-mer signatures
  - Useful in preprocessing steps for genomic analysis, clustering, or comparative genomics
  - Can be employed in workflows requiring dimensionality reduction of sequencing data
  - Command-line tool for extracting feature k-mers from large sequencing datasets

#### KmerSort (KmerSort.java)
**Purpose**: Abstract base class for k-mer based read sorting and processing in bioinformatics sequencing workflows.

**Core Function**: Provides a framework for k-mer based read manipulation, including preprocessing, sorting, clustering, and optional error correction or deduplication of sequencing reads.

**Key Features**:
- **K-mer Preprocessing**: Supports counting and filtering k-mers based on minimum occurrence threshold (via `preprocess()` method)
- **Parallel Processing**: Implements multi-threaded processing for hashing, sorting, and splitting reads using k-mer comparators
- **Read Transformation**: Supports multiple read processing modes including:
  - Error correction (`correct` flag)
  - Deduplication (`dedupe` flag)
  - Entry filtering
  - Read pair handling
- **Flexible Sorting**: Provides static methods for name-based and ID-based read sorting with pair-aware capabilities
- **Performance Tracking**: Comprehensive statistics tracking for reads processed, clumps formed, errors corrected, and duplicates removed
- **Configurable Parameters**: Supports various processing options like k-mer size, minimum count, quality quantization, and read name manipulation

**Usage**: Used as a base class in bioinformatics tools like Clumpify for preprocessing, filtering, and transforming sequencing reads based on k-mer characteristics.

#### KmerSort1 (KmerSort1.java)
**Purpose**: A command-line utility for sorting and processing k-mer based read data in bioinformatics sequencing files.

**Core Function**: Performs advanced read processing operations including sorting, deduplication, error correction, and read reordering using k-mer comparison strategies.

**Key Features**:
- **K-mer Sorting**: Uses KmerComparator to sort reads based on k-mer signatures with configurable comparison parameters
- **Deduplication**: Supports read deduplication across multiple passes with configurable seed and hash strategies
- **Read Reordering**: Offers multiple reordering modes including consensus-based and paired-end reordering
- **Multi-pass Processing**: Allows multiple processing passes with adaptive conservative strategies
- **Error Correction**: Provides optional read error correction mechanisms
- **Flexible Input Handling**: Supports single and paired-end FASTQ/SAM/BAM input formats with multi-group processing
- **Customizable Filtering**: Enables fine-tuned read filtering through various command-line parameters

**Usage**: Part of the BBTools suite for preprocessing and analyzing high-throughput sequencing data, particularly useful for removing duplicate reads, error correction, and preparing reads for downstream genomic analysis.

#### KmerSort2 (KmerSort2.java)
**Purpose**: Advanced k-mer sorting and processing utility for biological sequence data with multiple filtering and transformation capabilities.

**Core Function**: Implements a sophisticated k-mer based read processing pipeline that supports multiple operations including deduplication, condensing, error correction, and sorting of biological sequence reads.

**Key Features**:
- **Flexible Read Processing**: Supports multiple processing modes including deduplication, read condensing, and error correction
- **Multi-Pass Processing**: Allows multiple processing passes with configurable seed and comparison strategies
- **Parallel Processing**: Supports concurrent read input/output streams and multi-group processing
- **Adaptive Sorting**: Uses custom KmerComparator for advanced read sorting based on k-mer characteristics
- **Configurable Input/Output**: Handles various input file formats (FASTQ, SAM/BAM) with flexible output options
- **Advanced Filtering**: Supports prefiltering, entry filtering, and minimum count thresholds
- **Read Transformation**: Enables reverse complementing, read renaming, and consensus generation

**Usage**: Used in bioinformatics workflows for preprocessing, cleaning, and preparing sequencing reads for further analysis.

#### KmerSort3 (KmerSort3.java)
**Purpose**: Advanced k-mer sorting utility for processing and manipulating sequencing reads with multiple pass and clustering capabilities.

- **Core Function**: Perform multi-pass k-mer based sorting, clustering, and optional read processing (deduplication, condensing, error correction) for bioinformatics sequencing data.

- **Key Features**:
  - Multi-pass read processing with configurable sorting and clustering strategies
  - Concurrent read fetching and processing using multiple threads
  - Support for deduplication, read condensing, and error correction
  - Flexible input/output handling for FASTQ and paired-end sequencing data
  - Advanced memory management and file splitting for large datasets
  - Configurable k-mer comparison using custom KmerComparator
  - Entry filtering to reduce memory load and remove duplicate entries

- **Usage**: Used in bioinformatics workflows for preprocessing sequencing reads, including removing duplicate reads, condensing redundant reads, performing error correction on sequencing data, and preparing reads for downstream analysis.

#### KmerSplit (KmerSplit.java)
**Purpose**: Implements a parallel k-mer based read splitting and distribution utility for bioinformatics sequencing data processing.

**Core Function**: Splits input sequencing reads into multiple groups based on k-mer hash values, enabling parallel processing and data partitioning of large genomic datasets.

**Key Features**:
- **Parallel Processing**: Uses multi-threaded architecture to split reads across multiple output groups simultaneously
- **K-mer Hashing**: Employs k-mer comparator to distribute reads into specified number of groups based on hash values
- **Flexible Input Handling**: Supports both single and paired-end FASTQ input files
- **Configurable Splitting**: Allows customization of k-mer size, minimum k-mer count, and number of output groups
- **Quality Control**: Optional read validation, name shrinking, and quality quantization
- **Error Overlap Detection**: Supports optional error correction overlap (ECCO) detection for paired-end reads
- **Adaptive Threading**: Dynamically adjusts processing based on available system threads

**Usage**: Primarily used in bioinformatics workflows for distributing large sequencing datasets across multiple processing units and enabling parallel processing of large genomic sequencing files.

#### PivotSet (PivotSet.java)
**Purpose**: Reduces reads to their feature k-mer for computational analysis and filtering.

**Core Function**: Processes input sequencing reads to generate a k-mer count array (KCountArray) that captures the frequency and distribution of specific k-mer sequences across the input data.

**Key Features**:
- **K-mer Extraction**: Extracts k-mers (fixed-length DNA sequence subsequences) from input FASTQ reads using a configurable k-mer length (default k=31)
- **Parallel Processing**: Utilizes multi-threaded processing to extract k-mers efficiently, with thread count dynamically set based on system resources
- **Memory-Aware Filtering**: Dynamically calculates usable memory and allocates appropriate resources for k-mer counting using a fraction-based memory allocation strategy
- **Paired-End Read Support**: Handles both single-end and paired-end sequencing data, with optional error correction (ECCO) for read merging
- **Flexible Configuration**: Supports command-line configuration of k-mer length, minimum k-mer count, and other processing parameters
- **Pivot Estimation**: Provides estimated unique k-mer counts, including those appearing more than once (pivot points)

**Usage**: Primarily used in bioinformatics workflows for preprocessing sequencing reads, generating k-mer frequency distributions, and supporting read clustering and reduction algorithms.

#### ReadKey (ReadKey.java)
**Purpose**: A key representation for reads in the clump package, used for comparing, organizing, and tracking read characteristics.

**Core Function**: Provides a comprehensive key for reads that enables advanced comparison, positioning, and tracking of genomic reads across different dimensions (k-mer, strand, position, and optical properties).

**Key Features**:
- **K-mer Based Comparison**: Generates and compares read keys based on k-mer sequences, supporting strand-specific comparisons
- **Positional Tracking**: Maintains precise read position information, including ability to flip and reposition reads
- **Optical Coordinate Support**: Captures flowcell-specific coordinates (lane, tile, x, y) for optical duplicate detection
- **Distance Calculation**: Provides methods to calculate physical distances between read keys, supporting tile and coordinate-based proximity analysis
- **Unique Molecular Identifier (UMI) Matching**: Supports UMI comparison with configurable mismatch tolerance
- **Serializable Comparison**: Implements Comparable and Serializable interfaces for flexible read key manipulation
- **Memory Efficiency**: Includes methods for key creation and memory management, with built-in OutOfMemory error handling

**Usage**: Primarily used in genomic sequencing data processing to identify and group reads with similar characteristics, detect optical and PCR duplicates, and enable precise read positioning in bioinformatics workflows.

#### Splitter (Splitter.java)
**Purpose**: A utility class for splitting genomic clumps by allele variants with sophisticated variant detection and correlation strategies.

**Core Function**: Analyze genetic read clusters (clumps) and split them into subgroups based on genetic variations, helping to separate reads with different allelic compositions.

**Key Features**:
- **Variant Detection**: Identifies genetic variants by examining base counts and quality scores at different positions in genomic reads
- **Correlation-Based Splitting**: Uses advanced correlation algorithms to find meaningful pivot points for splitting clumps
- **Depth-Aware Splitting**: Implements minimum depth and fraction thresholds to ensure statistically valid splits
- **Quality-Sensitive Variant Tracking**: Considers both base counts and quality scores when determining variant significance
- **Multi-Variant Support**: Can handle complex scenarios with multiple variants in a single clump
- **Read Location Mapping**: Provides methods to convert between clump and read coordinate systems

**Usage**: Part of the genomic data processing pipeline, typically used for error correction, allele separation, and improving the accuracy of genetic sequencing data analysis.

#### StreamToOutput (StreamToOutput.java)
**Purpose**: Manages streaming input and output operations for processing read data with optional name-based sorting and multi-output distribution.

**Core Function**: Efficiently routes reads from input streams to one or multiple output streams, with optional name-based sorting and kmer-based group distribution.

**Key Features**:
- Supports processing reads from single or multiple concurrent input streams
- Provides flexible routing of reads to single or multiple output streams
- Implements optional name-based sorting of reads before processing
- Supports kmer-based read distribution across multiple output groups
- Tracks total reads and bases processed during streaming
- Handles error states and stream management
- Provides two constructors for different input stream initialization strategies
- Supports incremental kmer comparator generation

**Usage**: Used in bioinformatics data processing pipelines for routing and distributing read data, particularly useful for tasks requiring parallel processing of sequencing reads and complex read processing workflows.