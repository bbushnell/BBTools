# Cluster Package

The cluster package provides tools for clustering and analyzing sequence reads based on k-mer similarity, enabling efficient grouping and analysis of related genomic sequences.

## Classes

#### Cluster (Cluster.java)
**Purpose**: A thread-safe k-mer based clustering utility for organizing and scoring genomic reads using multiple k-mer lengths.

**Core Function**: Manages cluster-specific read classification through sophisticated k-mer analysis, supporting dual k-mer strategies (large and small k-mers) for precise genomic read categorization.

**Key Features**:
- **Dual K-mer Strategy**: Supports two k-mer lengths (k1 and k2) for flexible and nuanced read clustering
- **Thread-Safe Atomic Tracking**: Uses AtomicLongArray for concurrent k-mer counting and probabilistic tracking
- **Comprehensive Read Metadata**: Tracks read count, base count, GC content, and depth across multiple dimensions
- **Scoring Mechanisms**: Multiple scoring modes including difference, root mean square, AND count, multiplication, and Kolmogorov-Smirnov function
- **Dynamic Probability Calculation**: Recalculates k-mer probabilities with smoothing to handle low-frequency k-mers
- **Flexible Initialization**: Configurable cluster identification and array sizes for customized clustering

**Usage**: Utilized in genomic read clustering pipelines to categorize and analyze sequencing reads based on k-mer composition, particularly useful for metagenomics, assembly, and read classification tasks. Part of BBTools' advanced bioinformatics toolkit for processing complex sequencing data.

#### ClusterTools (ClusterTools.java)
**Purpose**: Provides utility methods for k-mer analysis and clustering operations in bioinformatics sequence processing.

**Core Function**: Implements k-mer manipulation and comparison techniques, including canonical k-mer generation, frequency analysis, and similarity scoring for sequence clustering.

**Key Features**:
- **K-mer Conversion**: Converts DNA sequences to k-mer representations with canonical encoding
  - `toKmers()`: Generates an array of canonical k-mers from a byte array of DNA bases
  - `toKmerCounts()`: Creates k-mer count arrays for frequency analysis
- **Canonical K-mer Identification**: 
  - `maxCanonicalKmer()`: Finds the maximum canonical k-mer for a given k-mer length
  - Handles reverse complement encoding to ensure consistent k-mer representation
- **Similarity Scoring Methods**: 
  - `andCount()`: Calculates overlap between k-mer sets using atomic long array
  - `innerProduct()`: Computes inner product of k-mer frequencies
  - `absDif()`: Calculates absolute difference between frequency arrays
  - `rmsDif()`: Computes root mean square difference between frequency arrays
  - `ksFunction()`: Implements a custom similarity scoring function using logarithmic differences
- **Flexible Input Handling**: 
  - Methods support various input types (byte arrays, int arrays)
  - Handles edge cases like short sequences and N-bases in DNA

**Usage**: 
- Used in bioinformatics sequence clustering and analysis pipelines
- Supports k-mer-based sequence comparison and classification
- Primarily utilized in genomic and metagenomic research for:
  - Sequence similarity assessment
  - Clustering DNA/RNA sequences
  - Frequency-based sequence comparison

#### MergeReadHeaders (MergeReadHeaders.java)
**Purpose**: A utility for replacing read headers in bioinformatics sequencing files using a separate header text file.

**Core Function**: Processes input sequencing files (FASTQ/FASTA) and replaces their original read headers with headers from a specified text file, supporting both single and paired-end read processing.

**Key Features**:
- Supports single and paired-end read file processing
- Reads headers from an external text file to replace existing read headers
- Handles multiple input/output file formats (FASTQ, FASTA)
- Configurable with various command-line options for input/output file specification
- Multithreaded processing with concurrent input/output streams
- Tracks and reports reads and bases processed
- Supports compressed and uncompressed file formats
- Flexible header replacement with optional interleaved output

**Usage**: Used in bioinformatics workflows to standardize or modify read headers across sequencing datasets, such as:
- Renaming reads from different sequencing runs
- Anonymizing read identifiers
- Preparing reads for downstream analysis that requires specific header formats
- Preprocessing sequencing data for comparative genomics or metagenomic studies

#### ReadTag (ReadTag.java)
**Purpose**: Manages metadata and clustering-related information for sequencing reads, computing GC content and supporting read classification.

**Core Function**: Creates a specialized tag object for a sequencing read that captures key genomic characteristics and supports clustering operations through k-mer analysis and strand tracking.

**Key Features**:
- **GC Content Computation**: Calculates the count of G and C bases in the read during instantiation
- **Strand-Aware Read Access**: Provides methods `r1()` and `r2()` to access reads based on strand orientation
- **K-mer Analysis**: Supports multiple k-mer representation methods:
  - `kmerArray1()`: Generates sorted long k-mer arrays
  - `kmerArray2()`: Creates canonically-ordered k-mer count arrays
  - `kmerFreq2()`: Computes normalized k-mer frequency distribution
- **Metadata Tracking**: Stores clustering-related metadata including:
  - Initial and final cluster assignments
  - Read depth
  - Strand information
- **Serializable Support**: Implements `Serializable` interface for potential object persistence

**Usage**: Used in bioinformatics clustering workflows to tag and analyze sequencing reads, particularly for grouping reads with similar genomic characteristics. Primarily utilized in computational genomics for read classification, clustering, and comparative analysis of sequencing data.

#### ReclusterByKmer (ReclusterByKmer.java)
**Purpose**: A specialized tool for advanced read clustering and reclustering using k-mer analysis techniques.

**Core Function**: Performs sophisticated read clustering by analyzing k-mer spectra across genomic reads, enabling dynamic read classification and reorganization based on sequence similarity.

**Key Features**:
- **Dual K-mer Analysis**: Supports two k-mer lengths (k1 and k2) for flexible and nuanced clustering strategies
- **Concurrent Processing**: Utilizes multi-threaded processing for efficient handling of large read datasets
- **Adaptive Clustering Modes**: Supports multiple clustering modes including create, recluster, and refine
- **Ambiguity Handling**: Provides configurable strategies for managing reads with ambiguous cluster assignments (best match, random assignment, toss, etc.)
- **Flexible Input Support**: Handles both single-end and paired-end read inputs in various file formats (FASTQ)
- **Configurable Cluster Assignment**: Allows fine-tuning of cluster assignment through command-line parameters

**Usage**: Used in genomic data analysis workflows to categorize and reorganize reads based on sequence composition, particularly useful in:
- Metagenomics research
- Genome assembly preprocessing
- Sequence classification and binning
- Reducing redundancy in sequencing datasets