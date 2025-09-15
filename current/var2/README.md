# BBTools var2 Package

Advanced variant calling and genomic variation analysis framework providing comprehensive tools for detecting, filtering, analyzing, and manipulating genetic variants from high-throughput sequencing data, including neural network-based quality assessment and multi-sample variant merging capabilities.

---

## AnalyzeVars (AnalyzeVars.java)
**Purpose**: Utility class for analyzing and manipulating variants in reads
**Core Function**: Processes genomic reads to detect, fix, and filter known genetic variants
**Key Features**:
- Identifies substitutions, insertions, and deletions in reads
- Converts variant markers to indicate known variant status
- Filters variants based on coverage and allele depth
- Loads and processes variants from VCF files
**Usage**: Supports variant calling and analysis in genomic sequencing data

## ApplyVariants (ApplyVariants.java)
**Purpose**: Applies variants to genomic sequences based on VCF files and optional depth information
**Core Function**: Processes input reads, applies genetic variants from VCF files with depth-based filtering
**Key Features**:
- Supports filtering variants by minimum coverage depth
- Handles insertions, deletions, and substitutions
- Masks low-coverage regions with 'N'
- Configurable variant application rules
**Usage**: Modifies genomic sequences by applying genetic variants with optional depth-based quality control

## CallVariants (CallVariants.java)
**Purpose**: Calls variants from one or more SAM or BAM files using multithreaded processing.
**Core Function**: Performs comprehensive genetic variant discovery and analysis across input sequencing data files.
**Key Features**:
- Multithreaded variant calling with configurable thread pools
- Optional prefiltering using Bloom filter-like structures
- Read realignment and quality-based trimming
- Comprehensive variant statistics and filtering
- Support for forced variants from input VCF files
- Multiple output formats (VAR, VCF, GFF)
**Usage**: Genetic variant identification and analysis for genomic sequencing data, enabling detailed characterization of genetic variations across entire datasets.

## CallVariants2 (CallVariants2.java)
**Purpose**: Calls variants from one or more SAM or BAM files for multi-sample analysis.
**Core Function**: Implements a comprehensive multi-threaded variant calling pipeline that processes aligned reads, discovers variants, and generates VCF output.
**Key Features**:
- Supports multisample variant discovery with parallel processing
- Optional k-mer prefiltering to reduce memory usage for large cohorts
- Neural network-based variant quality scoring
- Configurable read realignment and quality trimming
- Detailed variant statistics and filtering
**Usage**: Used for large-scale genomic variant identification across multiple samples, supporting complex variant calling workflows with flexible configuration options.

## CompareVCF (CompareVCF.java)
**Purpose**: Perform set operations on multiple VCF (Variant Call Format) files
**Core Function**: Compares and processes VCF files using set operations like difference, union, and intersection
**Key Features**:
- Supports difference, union, and intersection set operations on VCF files
- Can split complex variants into components
- Filters variants by quality score
- Handles multiple input VCF files
**Usage**: Compare genetic variants across multiple samples or experiments

## CVOutputWriter (CVOutputWriter.java)
**Purpose**: Utility class for writing CallVariants output files with static methods for generating variant output formats.
**Core Function**: Manages writing variant-related output files including Var, VCF, and GFF formats, and generates histograms for variant statistics.
**Key Features**:
- Writes output files for variant analysis
- Generates histograms for score, zygosity, and quality distributions
- Supports multiple output file formats (Var, VCF, GFF)
- Calculates statistical summaries like mean, median, and mode
**Usage**: Used in variant calling pipelines to output and analyze genetic variation data

## FeatureVectorMaker (FeatureVectorMaker.java)
**Purpose**: Converts Var objects into feature vectors for neural network analysis.
**Core Function**: Implements a factory pattern for generating specialized feature vectors using intern-developed neural network models.
**Key Features**:
- Three feature extraction modes: ELBA (32 features), LAWRENCE (8 features), DONOVAN (16 features)
- Supports quality score prediction and genotype calling
- Implements min-max scaling and normalization techniques
- Extensible framework for variant feature extraction
**Usage**: Generates standardized feature vectors for variant analysis in neural network models, allowing flexible feature vector creation for different analysis tasks

## FilterSam (FilterSam.java)
**Purpose**: Removes lines with unsupported variations from a sam file.
**Core Function**: Processes SAM/BAM files, filtering reads based on variant criteria and quality metrics.
**Key Features**:
- Supports filtering reads with substitutions, insertions, and deletions
- Configurable maximum variant thresholds
- Parallel processing with multithreading support
- Separate output streams for good and bad reads
**Usage**: Used for quality control and variant filtering in genomic sequencing data processing

## FilterVCF (FilterVCF.java)
**Purpose**: Filters VCF files based on variant quality, type, position, and statistical criteria.
**Core Function**: Process VCF files through comprehensive filtering, supporting both single-threaded and multithreaded operations.
**Key Features**:
- Statistical filtering using VarFilter and SamFilter criteria
- Position-based and variant type filtering
- Multi-allelic variant splitting capabilities
- Quality score histogram generation
**Usage**: Used for preprocessing and cleaning variant call files (VCF) before downstream genomic analysis

## MergeSamples (MergeSamples.java)
**Purpose**: Merges VCF files from multiple samples into a unified variant call set.
**Core Function**: Synchronously reads and combines variant calls across multiple samples, preserving individual sample information and aggregating statistical evidence.
**Key Features**:
- Processes VCF files with identical genomic positions
- Supports multithreaded merging using producer-consumer architecture
- Aggregates statistical metrics like read counts, quality scores, and coverage
- Optionally applies neural network-based variant scoring
**Usage**: Used in multi-sample variant calling workflows to create comprehensive merged variant lists with per-sample details preserved

## Realign (Realign.java)
**Purpose**: Realigns samlines to a reference.
**Core Function**: Processes SAM/BAM files to realign reads against a reference genome, applying quality and mapping filters.
**Key Features**:
- Supports multi-threaded read processing
- Applies quality trimming and border trimming
- Filters reads based on mapping quality and pair status
- Attempts realignment of reads to improve mapping accuracy
**Usage**: Used for improving read alignments in genomic sequencing data, particularly for enhancing mapping precision of sequenced reads.

## Realigner (Realigner.java)
**Purpose**: Realigns reads using Multiple Sequence Alignment (MSA) to improve alignments from non-affine-gap aligners.
**Core Function**: Performs glocal alignment with padding around the original alignment region, retaining only realignments that improve alignment score.
**Key Features**:
- Handles reads with indels longer than 1bp
- Analyzes alignment quality before realignment
- Supports unclipping of terminal indels
- Tracks realignment metrics (attempts, successes, improvements)
**Usage**: Improves alignment accuracy for genomic reads with complex alignment patterns, especially those poorly aligned by simpler algorithms.

## SamFilter (SamFilter.java)
**Purpose**: Filters SAM/BAM alignments, VCF variants, and genomic data based on configurable criteria.
**Core Function**: Applies multi-dimensional filtering to genomic data using mapping quality, position ranges, alignment identity, and SAM flags.
**Key Features**:
- Supports flexible filtering of mapped/unmapped reads
- Configurable mapping quality thresholds
- Alignment identity range filtering
- Contig name whitelist filtering
**Usage**: Used to selectively process genomic alignment data during bioinformatics analysis and preprocessing.

## Scaffold (Scaffold.java)
**Purpose**: Represents a single scaffold (chromosome/contig) in a reference genome, handling coverage tracking and sequence storage.
**Core Function**: Manages scaffold metadata, tracks read coverage, and supports parsing SAM header lines for genome reference information.
**Key Features**:
- Supports lazy initialization of coverage arrays
- Optional strand-specific coverage tracking
- Handles coordinate conversion between SAM and internal representations
- Calculates coverage for different variant types
**Usage**: Used in genomic variant analysis and reference genome coverage calculations, particularly in bioinformatics sequencing pipelines.

## ScafMap (ScafMap.java)
**Purpose**: Maps scaffold (chromosome/contig) names to Scaffold objects with efficient lookup capabilities.
**Core Function**: Provides a mapping system for tracking scaffolds from different genomic file formats, supporting SAM, VCF, and FASTA references.
**Key Features**:
- Supports loading scaffold information from multiple file headers
- Handles alternative scaffold name mappings with whitespace variations
- Tracks scaffold metadata like name, length, and coverage
- Provides methods for retrieving scaffolds by name or numeric ID
**Usage**: Used for managing and querying reference genome scaffold information during bioinformatics data processing.

## SoftClipper (SoftClipper.java)
**Purpose**: Performs soft-clipping operations on sequence alignments by identifying poorly aligned terminal regions.
**Core Function**: Uses dynamic scoring to find optimal alignment regions and clip suboptimal terminal bases.
**Key Features**:
- Computes alignment quality with scoring system for matches, substitutions, insertions, and deletions
- Dynamically finds highest-scoring contiguous alignment region
- Handles complex cases with deletions and insertions
- Adjusts alignment coordinates during clipping
**Usage**: Used in genomic sequence alignment to improve alignment quality by removing unreliable terminal bases.

## Var (Var.java)
**Purpose**: Represents a single genomic variant with comprehensive quality metrics and statistical evidence.
**Core Function**: Core variant calling class that accumulates evidence from multiple reads, calculates statistical quality scores, and outputs results in standard formats (VAR, VCF).
**Key Features**:
- Multi-threaded variant accumulation from read alignments
- Detailed statistical scoring for variant confidence
- Bias detection (strand, read, positional)
- Support for multiple variant types: substitutions, insertions, deletions
**Usage**: Used for processing and analyzing genomic variants from sequencing data, providing detailed variant calling and quality assessment.

## VarFilter (VarFilter.java)
**Purpose**: Comprehensive filtering system for genetic variants based on multiple quality metrics
**Core Function**: Evaluates variant quality through multi-tier statistical filtering using depth, quality, proximity, and integrated scoring
**Key Features**:
- Implements 7-tier filtering approach with early rejection strategies
- Supports configurable thresholds for read depth, quality scores, and variant frequency
- Handles proximity-based filtering for nearby variants
- Integrates optional neural network scoring for advanced variant assessment
**Usage**: Used in genetic variant analysis to distinguish genuine variants from sequencing artifacts or low-confidence calls

## VarHelper (VarHelper.java)
**Purpose**: Utility class providing helper methods for variant processing and output formatting.
**Core Function**: Generates headers, calculates scores, analyzes homopolymers, and processes junction variants for variant calling.
**Key Features**:
- Generates VAR and VCF format headers with detailed sequencing metadata
- Converts variant quality scores to Phred scale
- Counts homopolymer lengths around substitution positions
- Detects junction variants from clipped read alignments
**Usage**: Supporting variant analysis in BBTools by providing utility functions for variant detection and reporting

## VarKey (VarKey.java)
**Purpose**: Simplified key representation of genetic variants for efficient hashing and comparison
**Core Function**: Creates lightweight variant identifiers with minimal information storage
**Key Features**:
- Converts full Var objects to compact variant keys
- Supports different variant types (SUB, INS, DEL)
- Implements custom hashCode and comparison methods
- Efficiently stores scaffold, position, length, type, and first allele nucleotide
**Usage**: Used for variant identification and storage in hash-based collections when full variant details are unnecessary

## VarMap (VarMap.java)
**Purpose**: Thread-safe container for managing large collections of genomic variants using sharded concurrent hash maps.
**Core Function**: Provides efficient, multithreaded storage, retrieval, and processing of genomic variants with minimal lock contention.
**Key Features**:
- Sharded ConcurrentHashMap storage for scalable variant management
- Multithreaded variant processing with parallel statistical analysis
- Custom iterator for seamless traversal across variant shards
- Advanced filtering and quality assessment pipeline
**Usage**: Used in genomic variant calling and analysis, supporting complex multi-sample variant detection workflows.

## VarProb (VarProb.java)
**Purpose**: Statistical probability calculations for variant calling quality assessment.
**Core Function**: Provides binomial probability analysis to evaluate the likelihood of observed variant patterns occurring by chance.
**Key Features**:
- Calculates adjusted probability of binomial events with sequencing bias
- Precomputes factorial and binomial coefficient matrices for efficiency
- Handles large value calculations with iterative precision management
- Supports cumulative probability matrix generation
**Usage**: Used in variant calling to distinguish real genetic variants from sequencing artifacts or random noise.

## VCFFile (VCFFile.java)
**Purpose**: Loads and parses VCF (Variant Call Format) files into memory
**Core Function**: Processes VCF files by parsing headers, extracting sample names, and loading variant lines into a LinkedHashMap
**Key Features**:
- Handles VCF header parsing and sample name extraction
- Supports complex variant splitting 
- Creates scaffold mapping from contig headers
- Provides methods for variant line retrieval
**Usage**: Used for VCF file comparison and processing operations, particularly in genomic variant analysis

## VCFLine (VCFLine.java)
**Purpose**: Represents a single line from a VCF (Variant Call Format) file, handling complex parsing and manipulation of genetic variant data.
**Core Function**: Parses and manages variant information including scaffold, position, reference, alternate alleles, quality scores, and sample details.
**Key Features**:
- Supports parsing multiple variant types: substitutions, insertions, deletions, multi-allelic variants
- Implements advanced splitting methods for complex variant representations
- Handles canonical trimming of variant representations
- Provides type detection and classification for genetic variants
**Usage**: Used in genomic variant analysis to process and manipulate variant information from VCF files, supporting detailed variant examination and comparison.

## VcfLoader (VcfLoader.java)
**Purpose**: Multithreaded loader for VCF and VAR format files using producer-consumer pattern.
**Core Function**: Parses variant files with parallel processing to maximize I/O and CPU efficiency.
**Key Features**:
- Supports both VAR and VCF file formats
- Uses producer-consumer threading model
- Handles header and data line parsing
- Supports optional coverage and extended information parsing
**Usage**: Load variant files with automatic format detection and parallel processing.

## VcfToVar (VcfToVar.java)
**Purpose**: Utility class for parsing VCF format data into Var objects with complex coordinate and allele conversion.
**Core Function**: Converts VCF text lines into Var objects, handling chromosome mapping, coordinate normalization, and statistical field extraction.
**Key Features**:
- Supports basic and extended VCF parsing modes
- Handles indel coordinate normalization
- Extracts coverage and statistical information
- Resolves scaffold names to internal numbering
**Usage**: Converts VCF variant data for genomic variant tracking and analysis, used in variant calling and genome comparison workflows.

## VcfWriter (VcfWriter.java)
**Purpose**: Multithreaded writer for variant data in VCF, VAR, or GFF formats
**Core Function**: Uses producer-consumer pattern to format and write variant data across multiple threads while maintaining file order
**Key Features**:
- Supports writing in VCF, VAR, and GFF output formats
- Handles variant filtering and statistical metadata inclusion
- Uses ArrayBlockingQueue for parallel variant processing
- Provides format-specific variant output methods
**Usage**: Writing variant data files from genomic analysis, supporting different file format outputs

## VectorDonovan (VectorDonovan.java)
**Purpose**: Implements feature vector generation and normalization for variant quality score prediction
**Core Function**: Generates a 16-element normalized feature vector for neural network classification of genetic variants
**Key Features**:
- Quantile transformation of variant features
- Log and power transformations for feature scaling
- Caches transformers for performance optimization
- Handles 16 distinct variant characteristics
**Usage**: Preprocessing genetic variant data for machine learning models, specifically neural network classification

## VectorElba (VectorElba.java)
**Purpose**: Feature vector generation for quality score prediction of genetic variants
**Core Function**: Creates a normalized 11-dimensional feature vector for machine learning quality prediction
**Key Features**:
- Transforms input variant metrics using log and power functions
- Normalizes features to standard scale using min-max normalization
- Encodes variant type, allele fraction, and read characteristics
- Handles variant metrics like strand bias, read length, and mapping quality
**Usage**: Used in variant quality assessment and machine learning models for genetic variant analysis