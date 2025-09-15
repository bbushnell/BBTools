# Consensus Package

The consensus package provides advanced consensus sequence generation and gap filling algorithms for genomic assembly and sequence reconstruction workflows.

## Classes

#### BaseGraph (BaseGraph.java)
**Purpose**: Constructs a graph-based representation of sequence alignments to enable consensus sequence generation and variant analysis across reference sequences.

**Core Function**: Create a probabilistic graph structure that tracks base frequencies, quality scores, and alignment variants (insertions, deletions, substitutions) for generating high-confidence consensus sequences.

**Key Features**:
- **Dynamic Graph Construction**: Builds a flexible graph structure (`ref` and `del` arrays) representing reference sequence and possible variant paths
- **Quality-Weighted Alignment Processing**: Incorporates base quality scores and mapping quality into graph node weightings using sophisticated scoring algorithms
- **Variant Tracking**: Captures multiple alignment variants including substitutions, insertions, and deletions at each graph node
- **Consensus Generation**: Implements `traverse()` method to generate optimal consensus sequence by selecting most probable bases and accounting for depth and quality
- **Alignment Scoring**: Provides complex scoring methods (`score()` and `scoreOld()`) to evaluate alignment quality and confidence
- **Orientation Detection**: Supports finding the optimal read orientation through `findBestOrientation()` method
- **Flexible Padding**: Allows adding padding to reference sequences for robust edge case handling

**Usage**: 
- Used in bioinformatics pipelines for generating high-confidence consensus sequences from multiple sequence alignments
- Particularly useful in genomic variant calling, sequencing error correction, and de novo assembly workflows
- Supports processing of reads with complex alignment patterns, handling insertions, deletions, and substitutions
- Can be integrated into larger genomic analysis tools that require precise sequence reconstruction

#### BaseGraphPart (BaseGraphPart.java)
**Purpose**: Abstract base class for graph elements (nodes and edges) in a consensus graph representation, providing fundamental type classification and serialization support.

**Core Function**: Manages graph element type classification with strict validation, ensuring graph elements are categorized as Reference (REF), Insertion (INS), or Deletion (DEL) types.

**Key Features**:
- **Type Classification**: Supports three discrete graph element types (REF, INS, DEL) with runtime assertion validation
- **Serialization Support**: Implements Serializable interface for persistent storage and object transmission
- **Type String Conversion**: Provides `typeString()` method to convert numeric type to human-readable string representation
- **Abstract Part Identification**: Enforces `partString()` method in subclasses for specific element identification
- **Immutable Type Storage**: Stores graph element type as a final integer field for type consistency
- **Inheritance Framework**: Serves as a base class for BaseEdge and BaseNode, establishing common type management behavior

**Usage**: 
- Used as a parent class for graph structural elements in consensus-based graph algorithms
- Provides a standardized way to classify and manage graph element types across different graph-related operations
- Enables type-safe graph element creation and manipulation in consensus graph implementations

#### BaseNode (BaseNode.java)
**Purpose**: Represents a base position in a consensus graph, tracking statistical information about base observations and supporting consensus generation for genomic sequences.

**Core Function**: Maintain a statistical representation of base frequencies and qualities at a specific genomic position, enabling consensus base determination through sophisticated tracking and probabilistic analysis.

**Key Features**:
- **Base Frequency Tracking**: Maintains integer arrays `acgtCount` and `acgtWeight` to track occurrence counts and weighted qualities for A, C, G, and T bases
- **Consensus Calculation**: Implements advanced consensus determination method (`consensus()`) that considers reference base, base frequencies, and quality scores
- **Probabilistic Base Scoring**: Provides methods like `baseProb()` and `baseScore()` to calculate base probabilities and scoring
- **Position-Aware Representation**: Stores reference position (`rpos`) and reference base, enabling context-specific analysis
- **Edge Tracking**: Maintains optional references to reference (`refEdge`), insertion (`insEdge`), and deletion (`delEdge`) nodes for graph-based traversal
- **Allele Differentiation**: Includes `alleleDif()` method to calculate the difference between the most and second most common bases
- **Flexible Node Typing**: Supports different node types (reference, insertion, deletion) through type parameter in constructor

**Usage**: 
- Part of a consensus graph generation algorithm in genomic sequence analysis
- Used in bioinformatics tools for:
  1. Determining the most likely base at a specific genomic position
  2. Tracking base frequencies across multiple sequence reads
  3. Supporting variant calling and sequence alignment processes
  4. Providing probabilistic base quality assessments

#### ConsensusMaker (ConsensusMaker.java)
**Purpose**: Generates consensus sequences from aligned reads by analyzing and aggregating genomic sequence data across multiple input alignments.

**Core Function**: Multi-threaded reference sequence modification algorithm that processes SAM/BAM or FASTA alignments to generate consensus sequences with variant tracking and quality assessment.

**Key Features**:
- **Multi-threaded Alignment Processing**: Uses concurrent input/output streams and worker threads to process large genomic datasets efficiently
- **Flexible Input Handling**: Supports multiple input formats including SAM, BAM, and FASTA with configurable filtering and alignment parameters
- **Variant Tracking**: Records substitution, deletion, and insertion counts during consensus generation
- **Identity and Quality Scoring**: Calculates read alignment identity and optional model scoring for quality assessment
- **Comprehensive Statistics Reporting**: Generates detailed output including reads processed, bases processed, average identity, and variant counts
- **Configurable Consensus Generation**: Supports parameters like ploidy, minimum depth, and trimming options for refined consensus building
- **Optional Model-based Scoring**: Can incorporate pre-existing sequence models for advanced consensus generation

**Usage**: Used in bioinformatics workflows for:
- Generating consensus sequences from genomic alignment data
- Analyzing sequence variations across multiple read alignments
- Creating representative sequences from complex sequencing datasets
- Quality assessment of genomic alignments

#### ConsensusObject (ConsensusObject.java)
**Purpose**: Abstract base class for consensus-related data structures in the BBTools consensus package, providing core infrastructure for consensus object representation and manipulation.

**Core Function**: Serve as a foundational abstract class that defines key methods and constants for consensus-related objects, with a focus on text representation and type classification for genomic consensus processing.

**Key Features**:
- **Abstract Text Conversion**: Defines an abstract `toText()` method using `ByteBuilder` for converting consensus objects to text representation
- **Type Constants for Genomic Variations**: Defines static integer constants for genomic variation types:
  - `REF=2`: Match/Substitution or neutral-length node/edge
  - `INS=1`: Insertion node or edge
  - `DEL=0`: Deletion edge to non-adjacent node
- **Static Configuration Parameters**: Provides configurable static parameters for consensus processing:
  - Minimum depth threshold (`minDepth`)
  - Minor Allele Frequency (MAF) thresholds for substitutions, deletions, insertions
  - Trimming and conversion options (`trimDepthFraction`, `onlyConvertNs`)
- **Identity and Mapping Quality Flags**: Include configuration for identity processing:
  - `useMapq`: Toggle mapping quality usage
  - `invertIdentity`: Invert identity calculation
  - `identityCeiling`: Set maximum identity threshold (default 150)
- **Verbose Debugging**: Static `verbose` flag for enabling detailed logging
- **Consistent `toString()` Implementation**: Final method that leverages the abstract `toText()` method for string conversion

**Usage**: 
- Serves as a base class for more specific consensus object implementations in the BBTools consensus package
- Provides a standardized framework for representing and processing genomic consensus data
- Useful in bioinformatics workflows that require detailed tracking of genomic variations, such as variant calling, sequence alignment, and mutation analysis
- Allows subclasses to implement custom text representation while maintaining a consistent interface

#### FixScaffoldGaps (FixScaffoldGaps.java)
**Purpose**: Resizes scaffold gaps to represent the best estimate based on the insert size distribution of paired reads in genomic sequencing data.

**Core Function**: Dynamically adjusts the number of "N" bases (gap regions) in genomic scaffolds by analyzing read alignment depth, insert sizes, and coverage distribution.

**Key Features**:
- **Multi-Threaded Insert Size Analysis**: Uses concurrent processing to calculate insert size distribution across multiple read alignment threads
- **Dynamic Gap Modification**: Intelligently widens or narrows scaffold gaps based on statistical insert size proxies
- **Depth-Aware Gap Adjustment**: Considers read depth and local coverage when modifying scaffold gaps
- **Trim-Based Noise Reduction**: Applies configurable read edge trimming to reduce alignment noise (default 40% trim)
- **Comprehensive Scaffold Gap Tracking**: Maintains detailed statistics on gap modifications, including unchanged, widened, and narrowed gaps
- **Flexible Input Handling**: Supports various input formats (SAM/BAM) with configurable filtering and processing options

**Usage**: Used in genomic assembly and scaffolding workflows to improve scaffold representation by estimating and adjusting gap regions based on empirical read alignment data. Particularly useful for:
- Refining draft genome assemblies
- Improving scaffold contiguity estimates
- Preparing genomes for downstream comparative or functional analyses

#### Lilypad (Lilypad.java)
**Purpose**: Scaffolds contigs based on paired read mapping by using paired-end sequencing alignment data to connect and extend genomic contigs.

**Core Function**: Performs genomic scaffolding by analyzing SAM/BAM alignment files, identifying connections between contigs through paired-end read mappings, and constructing larger genomic sequences.

**Key Features**:
- **Multi-threaded Contig Scaffolding**: Uses parallel processing to analyze SAM file reads across multiple threads, enabling efficient large-scale genomic reconstruction
- **Sophisticated Edge Detection**: Employs complex edge-finding algorithms to determine optimal contig connections based on read mapping weight, strand consistency, and insert size
- **Adaptive Insert Size Estimation**: Calculates and tracks insert size distribution using percentile-based histogram analysis
- **Flexible Filtering Mechanism**: Supports configurable SAM line filtering with parameters like minimum mapping quality, strand consistency, and depth requirements
- **Dynamic Gap Filling**: Intelligently adds 'N' bases between contigs, with length inferred from insert size distribution
- **Strand-aware Contig Manipulation**: Supports strand-specific contig flipping and orientation preservation during scaffolding
- **Comprehensive Statistics Tracking**: Maintains detailed metrics about reads processed, scaffolds generated, gaps added, and insert size characteristics

**Usage**: 
- Genomic assembly pipeline for connecting short DNA sequence contigs into longer, more complete genomic sequences
- Particularly useful in metagenomic and de novo genome assembly workflows where paired-end sequencing data is available
- Command-line tool for bioinformaticians and genomic researchers working with complex sequencing datasets
