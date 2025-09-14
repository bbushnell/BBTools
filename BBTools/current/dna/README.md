# dna Package - DNA and Genomic Data Processing
*Comprehensive utilities for genomic sequence manipulation, amino acid translation, and chromosome data management*

## Package Overview
The `dna` package provides essential classes for genomic data processing, including DNA/RNA sequence manipulation, amino acid translation, chromosome array management, gene modeling, motif analysis, and scaffold handling.

---

## Core Classes

### Sequence Processing and Translation

#### AminoAcid (dna/AminoAcid.java)
**Purpose**: Comprehensive utility class for amino acid and nucleotide sequence conversion and manipulation in bioinformatics.

- **Core Function**: Provides advanced methods for translating between DNA, RNA, and amino acid representations, with robust encoding and decoding capabilities.

- **Key Features**:
  - **Nucleotide/Amino Acid Conversion**:
    - `stringToAAs()`: Converts DNA/RNA sequence to amino acid sequence
    - `toAAs()`: Translates byte arrays of nucleotides to amino acids
    - `toNTs()`: Converts amino acid sequences back to nucleotide sequences
    - `toAAQuality()`: Computes amino acid quality from nucleotide quality scores

  - **Encoding Methods**:
    - `kmerToString()`: Decodes bit-packed k-mer representations
    - `stringToKmer()`: Encodes DNA/RNA sequences into bit-packed representations
    - `kmerToStringAA()`: Decodes amino acid k-mer representations

  - **Complementary Sequence Operations**:
    - `reverseComplementBases()`: Generates reverse complement of DNA/RNA sequences
    - `complementBasesInPlace()`: In-place base complementation
    - `equalsReverseComp()`: Checks if sequences are reverse complements

  - **Biochemical Utility Methods**:
    - `codonToString()`: Converts numerical codon representations to nucleotide strings
    - `toAA()`: Translates codons to specific amino acids
    - `isFullyDefined()`: Validates base/amino acid definitions

- **Usage**: Primarily used in bioinformatics and genomic sequence processing, particularly for:
  - DNA/RNA sequence translation
  - Codon and amino acid manipulation
  - Sequence encoding and decoding
  - Nucleotide sequence validation and transformation

---

### Chromosome and Assembly Management

#### ChromArrayMaker (ChromArrayMaker.java)
**Purpose**: Creates chromosome arrays from genomic sequence files, supporting FASTA and FASTQ formats with advanced parsing and padding capabilities.

- **Core Function**: Converts genomic input files into structured chromosome arrays with configurable padding and scaffold management
  - Supports multiple input file formats (FASTA, FASTQ)
  - Generates chromosome-level data structures
  - Writes chromosome information to disk or memory

- **Key Features**:
  - Dynamic scaffold merging and padding strategies
  - Flexible genome build and chromosome numbering
  - Support for interchromosomal padding (start/mid/end padding)
  - Detailed chromosome and scaffold metadata generation
  - Configurable maximum chromosome length
  - Optional disk or in-memory chromosome array retention

- **Technical Implementation**:
  - Uses `CrisWrapper` for input file stream processing
  - Implements padding strategies with 'N' base insertions
  - Generates comprehensive chromosome metadata
  - Supports gzip compression for output files
  - Handles large genomic files with memory-efficient processing

- **Configuration Options**:
  ```java
  // Key configurable parameters
  public static int START_PADDING = 8000;   // Padding at chromosome start
  public static int MID_PADDING = 300;      // Padding between scaffolds
  public static int END_PADDING = 8000;     // Padding at chromosome end
  public static int MIN_SCAFFOLD = 1;       // Minimum scaffold length
  public static boolean MERGE_SCAFFOLDS = true;  // Allow scaffold merging
  ```

#### ChromosomeArray (ChromosomeArray.java)
**Purpose**: Central data structure for chromosome-level genomic information storage and retrieval with memory-efficient access patterns.

- **Key Features**:
  - Memory-efficient chromosome data structures
  - Supports large genome builds with optimized access patterns
  - Handles chromosome metadata and genomic coordinate management

#### FastaToChromArrays2 (FastaToChromArrays2.java)
**Purpose**: Converts FASTA files into chromosome arrays with efficient byte-level processing and configurable genome parsing strategies.

- **Core Function**: Parses FASTA input files, transforms genomic sequences into ChromosomeArray objects with advanced padding and scaffold management capabilities.

- **Key Features**:
  - **High-Performance Parsing**: Uses ByteFile1 for faster and lower-memory FASTA processing
  - **Flexible Scaffold Handling**: 
    - Supports merging scaffolds
    - Configurable padding strategies (start, mid, end padding)
    - Minimum scaffold length filtering
  - **Detailed Genome Metadata Generation**:
    - Generates chromosome information files
    - Tracks scaffold, contig, and base statistics
    - Supports writing chromosome arrays to disk

- **Command-Line Example**:
  ```bash
  java FastaToChromArrays2 ecoli_K12.fa 1 \
    writeinthread=false \
    genscaffoldinfo=true \
    retain \
    maxlen=536670912 \
    writechroms=true \
    minscaf=1 \
    midpad=300 \
    startpad=8000 \
    stoppad=8000
  ```

---

### Gene and Genomic Annotation

#### Gene (Gene.java)
**Purpose**: Comprehensive representation of a genomic gene with detailed coordinate and transcriptional information.

- **Core Function**: 
  - Modeling detailed genomic features including transcription start/stop, coding regions, exons, and strand information
  - Provides robust methods for comparing, intersecting, and analyzing gene structures

- **Key Features**:
  1. **Coordinate System**
     - Tracks transcription start/stop positions (txStart, txStop)
     - Tracks coding region start/stop positions (codeStart, codeStop)
     - Supports both plus and minus strand representations

  2. **Exon Management**
     - Stores exon array with precise start/end coordinates
     - Calculates exon lengths (exonLength, exonCodeLength)
     - Supports strand-compensated exon access

  3. **Genomic Metadata**
     - Captures gene symbol, mRNA accession, protein accession
     - Tracks gene status, completeness, and source
     - Supports pseudogene and untranslated gene identification

  4. **Intersection and Proximity Methods**
     - Methods to check intersections with transcription, coding, and exon regions
     - Calculates distances between genomic intervals
     - Finds nearest splice sites and exons

- **Notable Methods**:
  - `merge()`: Combines gene information from multiple sources
  - `intersectsTx()`, `intersectsCode()`: Check genomic region overlaps
  - `toGeneRelativeOffset()`: Converts absolute to gene-relative coordinates
  - `findClosestExon()`: Locates nearest exon to a given genomic region

#### Exon (dna/Exon.java)
**Purpose**: Represents a genomic exon with precise genomic coordinate tracking and gene structure modeling.

- **Core Function**: Models individual exons within a chromosome, capturing critical genomic annotation details including start/end positions, chromosome, strand, and exon type (UTR/CDS).

- **Key Features**:
  - **Precise Coordinate Representation**: 
    - Stores start (`a`) and end (`b`) positions as integers
    - Supports chromosome identification via byte representation
    - Tracks genomic strand (+/-) information

  - **Exon Type Classification**:
    - `utr` boolean flag for untranslated region exons
    - `cds` boolean flag for coding sequence exons
    - Enables fine-grained genomic annotation tracking

  - **Spatial Relationship Methods**:
    - `intersects()`: Checks if points or ranges overlap
    - `crosses()`: Determines partial intersection
    - `contains()`: Verifies complete containment
    - `distToSpliceSite()`: Calculates proximity to splice sites

  - **Exon Manipulation**:
    - `merge()`: Combines overlapping exons
    - `length()`: Calculates exon length
    - Supports multiple constructor variations for flexible instantiation

#### GeneSet (dna/GeneSet.java)
**Purpose**: Represents a collection of genes with shared characteristics, providing comprehensive genomic set management and annotation capabilities.

- **Core Function**: 
  - Manages a group of genes sharing a common identifier, chromosome, and potentially other genomic properties
  - Tracks transcript-level metadata such as genomic range, strand orientation, and gene type characteristics

- **Key Features**:
  - **Comprehensive Gene Collection Management**:
    - Stores multiple genes in an ArrayList
    - Tracks total number of transcripts
    - Calculates minimum and maximum genomic coordinates

  - **Genomic Annotation Handling**:
    - Captures chromosome information
    - Tracks strand orientation
    - Identifies special gene types like pseudogenes and untranslated genes

  - **Intersection and Range Detection**:
    - Supports point and range intersection checks
    - Enables precise genomic location queries

---

### Motif Analysis and Pattern Matching

#### Motif (Motif.java)
**Purpose**: Abstract base class for DNA motif representation and pattern matching algorithms with advanced sequence analysis capabilities.

- **Core Function**: Provides a flexible framework for identifying and counting DNA sequence motifs with exact and extended matching strategies.

- **Key Features**:
  - Abstract base class for motif matching implementations
  - Supports exact and extended DNA sequence pattern matching
  - Provides methods to count motif occurrences in byte arrays and strings
  - Includes built-in base probability and inversion calculations
  - Supports motif naming, length, and center position tracking

- **Key Methods**:
  - `countExact(String/byte[], int start, int end)`: Count exact motif matches
  - `countExtended(String/byte[], int start, int end)`: Count extended (fuzzy) motif matches
  - `matchStrength(byte[], int)`: Calculate match strength for a given position
  - `percentile(float strength)`: Convert match strength to percentile score

#### MotifMulti (dna/MotifMulti.java)
**Purpose**: Implements a multi-motif pattern matching strategy for DNA sequence analysis, allowing flexible matching across multiple motif patterns.

- **Core Function**: Extends the base `Motif` class to handle multiple motif patterns simultaneously, enabling complex sequence matching strategies
  - Supports exact and extended pattern matching across a collection of motifs
  - Dynamically selects the strongest matching motif in a given sequence region

- **Key Features**:
  - Multiple Motif Handling: Accepts an array of `Motif` objects during construction
  - Flexible Matching Methods:
    - `matchesExactly()`: Checks if any sub-motif matches exactly at a specific sequence position
    - `matchesExtended()`: Performs extended (potentially more flexible) matching across sub-motifs
    - `matchStrength()`: Determines the maximum normalized matching strength among sub-motifs

#### MotifProbsN (dna/MotifProbsN.java)
**Purpose**: Probabilistic motif analysis class for DNA/RNA sequence pattern matching with advanced N-base handling capabilities.

- **Core Function**: 
  - Implements probabilistic matching of sequence motifs with support for N-base variations
  - Calculates match strengths and normalizes probabilities for sequence patterns
  - Extends the base Motif class with enhanced probabilistic analysis

- **Key Features**:
  1. **Probabilistic Sequence Matching**
     - Computes match strength for DNA/RNA sequences using position-specific probability matrices
     - Handles variable-length motif patterns with configurable center positions
     - Supports normalization of match probabilities

  2. **N-Base Handling**
     - Configurable N-base parameter allows flexible handling of ambiguous bases
     - Uses `baseProb` array to manage base probabilities for N-base variations
     - Adapts probability calculations to account for base uncertainties

  3. **Position Importance Weighting**
     - Implements `positionImportance()` method to calculate position-specific significance
     - Dynamically adjusts motif probabilities based on position deviation from baseline

#### MotifSimple (dna/MotifSimple.java)
**Purpose**: Simple DNA motif matching class for exact and extended pattern recognition in byte-encoded sequences

- **Core Function**: Provides efficient methods for matching DNA sequence patterns with support for both exact and extended (degenerate) base matching
  - Extends the abstract `Motif` class
  - Handles DNA base matching with case-insensitive and extended matching capabilities

- **Key Features**:
  - **Exact Matching**: `matchesExactly()` method for precise sequence comparisons
  - **Extended Matching**: `matchesExtended()` method supporting degenerate base matching
  - **Case-Insensitive Search**: Compares both uppercase and lowercase base representations
  - **Flexible Initialization**: Constructor supports creating motifs with optional center positioning

---

### Utility and Support Classes

#### Data (Data.java)
**Purpose**: Central utility class for managing genomic data processing and genome information in BBTools.

- **Core Function**: Provides static methods and data structures for handling chromosome, gene, scaffold, and genome-related operations, including loading, unloading, and querying genomic information.

- **Key Features**:
  - **Genome Information Management**: 
    - Tracks chromosome lengths, defined bases, scaffolds, and contigs
    - Supports multiple genome builds (e.g., hg18, hg19)
    - Dynamically loads and unloads chromosome and gene data

  - **Chromosome and Scaffold Navigation**:
    - Methods to retrieve chromosome and scaffold information
    - Supports locating genes and gene sets on specific chromosomes
    - Handles scaffold name and location mapping

  - **Gene Processing**:
    - Load and manage gene data for different chromosomes
    - Generate gene range matrices for various genomic contexts
    - Support for gene set and transcript lookups

#### Matrix (Matrix.java)
**Purpose**: A flexible genomic matrix management class for storing, retrieving, and manipulating float-based grid data associated with genomic operations.

- **Core Function**: 
  - Manages and provides access to genomic scoring matrices from various build versions and genomic contexts
  - Supports dynamic loading of matrices from text files using a static HashMap for efficient retrieval
  - Enables extraction of matrix sub-grids with configurable prefix and length parameters

- **Key Features**:
  - **Static Matrix Repository**: Maintains a global table of matrices loaded from multiple file sources
  - **Flexible Matrix Retrieval**: 
    - Static `get(String name)` method for accessing pre-loaded matrices
    - Supports case-insensitive matrix name lookup
  - **Sub-Grid Extraction**: 
    - `subGrid(int prefixLength, int length)` method allows extracting specific sections of a matrix grid

#### ScafLoc (ScafLoc.java)
**Purpose**: Lightweight data structure for tracking genomic scaffold location with chromosome, name, and position.

- **Core Function**: Stores and manages genomic coordinate information for scaffolds
- **Key Features**:
  - Tracks scaffold name as a String
  - Stores chromosome number as an integer
  - Captures precise location within chromosome as an integer
- **Usage**: Used in genomic sequence mapping and coordinate tracking systems

#### Scaffold (Scaffold.java)
**Purpose**: Represents a genomic scaffold with metadata for sequence assembly and tracking.

- **Core Function**: Manages genomic scaffold information, parsing SAM (Sequence Alignment/Map) format entries and providing essential metadata about genomic sequences.

- **Key Features**:
  - Multiple constructor strategies for parsing scaffold information from different input formats
  - Captures scaffold name, length, and assembly information
  - Implements `Comparable` for ordered scaffold handling
  - Tracks various hit and base count statistics for genomic sequences

- **Scaffold Metadata Fields**:
  - `name`: Unique scaffold identifier
  - `assembly`: Assembly build information
  - `length`: Total scaffold sequence length
  - `basecount`: Base composition array (A,C,G,T,N)
  - `gc`: GC content percentage
  - `basehits`, `readhits`, `fraghits`: Sequence hit tracking metrics

---

## Package Usage
The dna package serves as the foundation for genomic data processing in BBTools, providing essential utilities for:
- DNA/RNA sequence manipulation and amino acid translation
- Chromosome and assembly data management
- Gene annotation and genomic feature modeling
- DNA motif analysis and pattern matching
- Genomic coordinate and scaffold tracking
- Matrix-based scoring and analysis operations

## Dependencies
- Relies on BBTools core utilities for file I/O and data processing
- Integrates with shared data structures and utility classes
- Supports multiple genome builds and annotation formats

---
*Documentation generated using evidence-based analysis of source code*