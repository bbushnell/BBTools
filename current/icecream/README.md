# icecream Package - PacBio Long-Read Processing Tools
*Specialized tools for Pacific Biosciences sequencing data analysis and quality control*

## Package Overview
The `icecream` package provides comprehensive tools for processing PacBio (Pacific Biosciences) long-read sequencing data, including alignment, quality control, chimera detection, and read simulation.

---

## Alignment Infrastructure

### IceCreamAligner (IceCreamAligner.java)
**Purpose**: Abstract base class for sequence alignment algorithms with configurable implementations
**Core Function**: Provides factory method to select alignment strategy (JNI or Java) and defines abstract alignment methods for forward and short alignments
**Key Features**: 
  - Factory method for aligner selection based on bit configuration
  - Abstract alignment methods with score and ratio thresholds
  - Supports different alignment implementations
**Usage**: `IceCreamAligner aligner = IceCreamAligner.makeAligner(32);`
**Dependencies**: aligner.AlignmentResult, shared.Shared

### IceCreamAlignerJava (IceCreamAlignerJava.java)
**Purpose**: Java implementation of a sequence alignment algorithm for DNA/protein sequences
**Core Function**: Performs forward and short-form alignments using dynamic programming, computing optimal alignment scores between query and reference sequences
**Key Features**: 
- Supports early exit for low-scoring alignments
- Computes alignment score with match/substitution/insertion/deletion penalties
- Calculates alignment ratio against query length
**Usage**: Align byte[] query against byte[] reference with configurable scoring parameters
**Dependencies**: aligner.AlignmentResult

### IceCreamAlignerJNI (IceCreamAlignerJNI.java)
**Purpose**: Java Native Interface (JNI) implementation for high-performance sequence alignment
**Core Function**: Performs fast byte array alignment using native methods, supporting different integer and short array sizes with early exit optimizations
**Key Features**: 
- Native method calls for performance-critical alignment
- Supports both int and short array alignment strategies
- Implements aggressive and safe early exit conditions
**Usage**: Extends IceCreamAligner for specialized JNI-based alignment
**Dependencies**: aligner.AlignmentResult, shared.KillSwitch, shared.Shared

---

## Quality Control and Analysis

### IceCreamFinder (IceCreamFinder.java)
**Purpose**: Detects inverted repeats in PacBio reads to identify chimeric reads with missing adapters
**Core Function**: Analyzes sequencing reads to find structural anomalies like inverted repeats, adapter issues, and low-complexity regions
**Key Features**:
  - Supports multi-threaded processing of DNA sequencing reads
  - Performs adapter alignment and trimming
  - Filters reads based on entropy and read quality
  - Generates detailed statistical output about read processing
**Usage**: Command-line tool for preprocessing and filtering PacBio sequencing data
**Dependencies**: Uses Java libraries for file I/O, alignment, and stream processing

### IceCreamGrader (IceCreamGrader.java)
**Purpose**: Read and classify sequencing reads based on potential "ice cream" contamination
**Core Function**: Processes input FASTQ files, identifies and tracks reads with "ice cream" markers, generating quality and contamination statistics
**Key Features**: 
- Handles paired/single-end read processing
- Counts good and bad reads/bases
- Generates detailed processing report
**Usage**: `java icecream.IceCreamGrader input.fastq`
**Dependencies**: fileIO, stream, shared, tracker packages

---

## Data Processing and Simulation

### IceCreamMaker (IceCreamMaker.java)
**Purpose**: Generates synthetic PacBio reads with chimeric inverted repeats
**Core Function**: Simulates PacBio sequencing by creating artificial genome fragments with configurable error rates, missing adapters, and molecular variations
**Key Features**:
  - Generates reads with configurable length, error rate, and genome characteristics
  - Supports simulating missing adapters and inverted repeats
  - Allows multiple read passes and base calling with error modeling
**Usage**: `java IceCreamMaker in=input.fasta out=output.fastq`
**Dependencies**: fileIO, shared, stream, structures packages

### ReformatPacBio (ReformatPacBio.java)
**Purpose**: PacBio data reformatting and processing tool
**Core Function**: Processes and filters PacBio sequencing reads, supporting ZMW-aware subsampling, trimming, and quality control
**Key Features**:
  - Supports exact read/base/ZMW sampling
  - Entropy-based read filtering
  - Consensus sequence generation
  - Trimming of low-quality read ends
**Usage**: Command-line tool for preprocessing PacBio sequencing data
**Dependencies**: stream, fileIO, shared packages

---

## Data Structures and Utilities

### ZMW (ZMW.java)
**Purpose**: Container for reads from a single PacBio Zero-Mode Waveguide (ZMW)
**Core Function**: Extends ArrayList<Read> to manage and analyze PacBio sequencing reads, providing methods to calculate read lengths, median reads, and estimate sequencing passes
**Key Features**: 
- Count and analyze read bases within a ZMW
- Find median and longest reads
- Estimate number of sequencing passes
**Usage**: Tracking and processing PacBio sequencing reads in bioinformatics pipelines
**Dependencies**: java.util.ArrayList, stream.Read, structures.IntList

### ZMWStreamer (ZMWStreamer.java)
**Purpose**: Wrapper for read input streams that processes Zero-Mode Waveguide (ZMW) reads
**Core Function**: Streams and organizes PacBio reads into ZMW-specific lists, processing reads sequentially from either a ConcurrentReadInputStream or SamReadStreamer
**Key Features**: 
- Supports limiting total reads or ZMWs processed
- Thread-safe queue-based read management
- Handles both BAM/SAM and other file formats
**Usage**: Used in bioinformatics pipelines for processing PacBio sequencing data
**Dependencies**: stream.*, fileIO.*, shared.*

### PBHeader (PBHeader.java)
**Purpose**: Parses and extracts structured metadata from Pacific Biosciences sequencing run headers
**Core Function**: Breaks down complex sequencing header strings into component parts like run ID, ZMW ID, and movie coordinates
**Key Features**: 
- Parses headers with `/` and `_` delimiters
- Extracts run ID, ZMW ID, start/stop coordinates
- Handles headers with potential trailing whitespace
**Usage**: `PBHeader header = new PBHeader("m64021_190821_100154/102/32038_35649")`
**Dependencies**: shared.Parse

### ReadBuilder (ReadBuilder.java)
**Purpose**: Constructs and manipulates PacBio sequencing read metadata with multi-subread support
**Core Function**: Parses and reconstructs complex sequencing read information, handling movie start/stop coordinates, number of passes, subreads, and error rates
**Key Features**: 
  - Static parsing method from Read objects
  - Generates read headers with detailed metadata
  - Supports appending multiple read builders
**Usage**: Create read metadata for PacBio sequencing data processing
**Dependencies**: stream.Read, structures.ByteBuilder

### PolymerTrimmer (PolymerTrimmer.java)
**Purpose**: Trims byte arrays by detecting homogeneous symbol sequences with dynamic scoring
**Core Function**: Implements left and right side trimming of byte arrays by calculating a dynamic score for symbol matching, allowing flexible polymer sequence removal
**Key Features**: 
- Dynamic scoring mechanism for symbol trimming
- Configurable minimum polymer length and fraction
- Supports both left and right side trimming
**Usage**: Used in bioinformatics for sequence preprocessing and adapter removal
**Dependencies**: None (standalone utility class)

---

## Package Usage

The icecream package enables comprehensive PacBio sequencing data processing including:

- **Alignment and Scoring**: High-performance sequence alignment with JNI and Java implementations
- **Quality Control**: Chimera detection, adapter trimming, and contamination analysis  
- **Data Simulation**: Synthetic read generation for testing and validation
- **ZMW Processing**: Zero-Mode Waveguide read organization and analysis
- **Format Conversion**: PacBio-specific data reformatting and preprocessing

---
*Documentation generated using evidence-based analysis of source code*