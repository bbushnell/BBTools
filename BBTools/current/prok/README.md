# Prok Package - Prokaryotic Genome Annotation and Gene Prediction

BBTools prokaryotic package providing comprehensive tools for bacterial and archaeal genome annotation, gene prediction, ribosomal sequence processing, and statistical modeling of genomic features.

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