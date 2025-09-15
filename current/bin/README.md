# Bin Package

This package provides comprehensive genomic binning, taxonomic classification, and metagenomic analysis tools. It includes advanced clustering algorithms, refinement strategies, and similarity metrics for organizing and analyzing genomic contigs.

## AbstractRefiner (AbstractRefiner.java)
**Purpose**: Abstract base class for bin refinement strategies that analyze and potentially split impure computational bins.
**Core Function**: Provides a framework for bin refinement by offering abstract methods to refine input bins and validate split quality.
**Key Features**:
- Factory method for creating specialized refiners (CRYSTAL, GRAPH, EVIDENCE, ENSEMBLE)
- Abstract refinement method returning null or split bins
- Split validation logic checking bin size and composition
- Supports multiple refinement parameter configurations
**Usage**: Used in genomic binning to improve cluster purity by detecting and splitting incorrectly grouped contigs.

## AdjustEntropy (AdjustEntropy.java)
**Purpose**: Calculates and compensates for entropy variations across different GC content levels in genomic sequences.
**Core Function**: Generates random sequences with specified GC content and measures their k-mer entropy.
**Key Features**:
- Supports dynamic GC content range sampling
- Computes sliding window entropy using k-mer frequency
- Generates statistically robust entropy estimates
- Provides strandedness and entropy compensation metrics
**Usage**: Used for analyzing genomic sequence entropy and generating normalized entropy references.

## AllToAllVectorMaker (AllToAllVectorMaker.java)
**Purpose**: Generate vector representations of biological contigs for comparative analysis
**Core Function**: Generates stratified vector samples across contigs using configurable parameters
**Key Features**:
- Creates vector representations from contigs and contig clusters
- Supports configurable sampling rates and filtering criteria
- Handles multi-cluster and single-contig vector generation
- Tracks k-mer differences and genome similarity metrics
**Usage**: Machine learning preprocessing for genomic sequence comparison and taxonomic analysis

## Bin (Bin.java)
**Purpose**: Core data structure representing a genomic bin with metadata and quality metrics
**Core Function**: Stores bin composition, statistics, and provides methods for bin manipulation and comparison
**Key Features**:
- Tracks bin contamination, completeness, and quality metrics
- Supports bin merging and contig management
- Calculates bin statistics and taxonomic assignments
- Implements comparable interface for sorting
**Usage**: Central data structure for metagenomic binning workflows

## BinComparator (BinComparator.java)
**Purpose**: Provides a custom comparison strategy for Bin objects with multi-tier ranking.
**Core Function**: Compares Bin objects using contamination level, size, and ID to establish a consistent sorting order.
**Key Features**:
- Prioritizes bins with lower contamination levels
- Uses size as secondary sorting criterion
- Ensures deterministic ordering through ID comparison
**Usage**: Used for sorting and ranking Bin objects in collections, prioritizing cleaner and larger bins.

## BinMap (BinMap.java)
**Purpose**: Maps genomic features to bins for taxonomic and functional classification
**Core Function**: Provides mapping infrastructure between genomic elements and their assigned bins
**Key Features**:
- Efficient bin-to-feature mapping
- Supports multiple mapping strategies
- Handles taxonomic ID resolution
- Optimized for large-scale genomic datasets
**Usage**: Core mapping infrastructure for genomic binning workflows

## Binner (Binner.java)
**Purpose**: Implements advanced binning algorithms for clustering genomic contigs
**Core Function**: Performs multi-pass refinement of bin/cluster assignments using sophisticated similarity metrics
**Key Features**:
- Multiple comparison modes (refine, residue, purify)
- Advanced binning parameter configuration
- Multithreaded parallel processing
- Size-based clustering optimization
**Usage**: Genomic sequence binning and metagenomic contig classification

## BinObject (BinObject.java)
**Purpose**: Superclass for binner classes that provides core k-mer and taxonomic processing utilities
**Core Function**: Manages k-mer remapping, canonical k-mer generation, and statistical calculations for genomic binning
**Key Features**:
- K-mer quantization and canonical representation
- GC content mapping for genomic sequences
- Taxonomic ID parsing and resolution
- Multiple distance/difference calculation methods
**Usage**: Provides utility methods for genomic sequence analysis, binning, and taxonomic classification

## BinSketcher (BinSketcher.java)
**Purpose**: Handles sketches and taxonomic assignments for contigs and clusters.
**Core Function**: Manages multi-threaded sketching process for biological sequence elements with configurable processing modes.
**Key Features**:
- Supports bulk and single-element sketching modes
- Configurable thread scaling based on workload
- Taxonomy assignment at genus level
- Dynamically adjusts sketch creation based on object size
**Usage**: Used for creating and managing taxonomic sketches in bioinformatics sequence analysis

## BinStats (BinStats.java)
**Purpose**: Represents statistical metadata and quality classification for genomic bins
**Core Function**: Captures bin characteristics like size, contamination, completeness, and RNA marker gene counts
**Key Features**:
- Classifies bin quality into tiers (UHQ, VHQ, HQ, MQ, VLQ, LQ)
- Tracks bin metadata: size, contigs, taxonomy ID, GC content, depth
- Implements Comparable for sorting bins by size and name
- Validates RNA marker gene completeness for quality assessment
**Usage**: Used in metagenomic binning to evaluate and categorize genome assemblies

## BinStatsComparator (BinStatsComparator.java)
**Purpose**: Custom Comparator for comparing BinStats objects with specific sorting criteria
**Core Function**: Sorts BinStats objects by contamination level, size, and ID
**Key Features**:
- Prioritizes lower contamination first
- Sorts by descending bin size when contamination is equal
- Provides stable sorting by ID for equal contamination and size
**Usage**: Used for sorting and ranking BinStats collections based on contamination and size

## ChartMaker (ChartMaker.java)
**Purpose**: Generates statistical charts and plots for bin-related data analysis
**Core Function**: Writes chart files with bin statistics using ByteStreamWriter
**Key Features**:
- Creates charts from bin statistics with sorting
- Generates contamination histogram
- Supports multiple chart types (bin stats, completeness/contamination plot)
- Handles floating-point precision for metrics
**Usage**: Generates visualization data for genomic bin analysis and quality assessment

## Clade (Clade.java)
**Purpose**: Represents a taxonomic clade with k-mer frequency signatures and genome comparison statistics.
**Core Function**: Stores and analyzes taxonomic data using 1-5 mer frequency counts and genome composition metrics.
**Key Features**:
- Calculates canonical k-mer distributions from 1 to 5 mers
- Tracks taxonomic metadata like ID, level, name, and lineage
- Computes GC content, entropy, and sequence strandedness
- Supports merging multiple taxonomic entries
**Usage**: Used in taxonomic classification and comparative genomics for tracking sequence characteristics across different taxonomic groups

## CladeIndex (CladeIndex.java)
**Purpose**: Indexes taxonomic clades for efficient lookup and comparison operations
**Core Function**: Provides fast access to clade data through optimized indexing structures
**Key Features**:
- Efficient clade lookup by taxonomic ID
- Supports range queries and similarity searches
- Optimized for large taxonomic databases
- Thread-safe access patterns
**Usage**: Core indexing infrastructure for taxonomic classification workflows

## CladeLoader (CladeLoader.java)
**Purpose**: Loads fasta files with TID-labeled contigs to produce Clade record output with kmer frequencies
**Core Function**: Processes input files, extracts taxonomic information from reads, and generates Clade records with metadata
**Key Features**:
- Supports concurrent processing of input files
- Can load clades from sequence or clade-format files
- Handles 16S and 18S ribosomal sequence additions
- Flexible taxonomic ID resolution and mapping
**Usage**: Used for taxonomic classification and analysis of genomic sequences, particularly in metagenomics and microbial research

## CladeLoaderMF (CladeLoaderMF.java)
**Purpose**: Designed to load one clade per file in metagenomic bins using parallel processing.
**Core Function**: Loads and processes taxonomic clades from multiple files simultaneously, with optional per-contig clade generation.
**Key Features**:
- Parallel processing of multiple files using one thread per file
- Supports loading clades from different file formats (FASTA, Clade)
- Optional per-contig clade generation
- Entropy tracking during sequence processing
**Usage**: Used in metagenomic classification and taxonomic analysis of genomic sequences

## CladeSearcher (CladeSearcher.java)
**Purpose**: Searches taxonomic databases for clade matches and similarity comparisons
**Core Function**: Performs efficient similarity searches across taxonomic clade databases
**Key Features**:
- Fast similarity search algorithms
- Supports multiple comparison metrics
- Handles large-scale taxonomic databases
- Optimized for real-time classification
**Usage**: Core search engine for taxonomic classification and clade identification workflows

## Cluster (Cluster.java)
**Purpose**: K-mer-based clustering and scoring utility for DNA reads
**Core Function**: Thread-safe cluster representation with atomic k-mer counting and probabilistic scoring
**Key Features**:
- Supports two k-mer lengths (k1, k2) for matching
- Atomic arrays for thread-safe k-mer tracking
- Dynamic read scoring based on k-mer probabilities
- Calculates read cluster metrics like GC content and depth
**Usage**: Used in bioinformatics for read classification and clustering operations

## Comparison (Comparison.java)
**Purpose**: Compares two Clade objects to determine their genomic similarity.
**Core Function**: Multi-method comparison of genomic profiles using k-mer frequency analysis.
**Key Features**:
- Supports 5 comparison methods: absolute difference, cosine, Hellinger, Euclidean, GC-compensated
- Uses 3, 4, and 5-mer k-mer frequency profiles
- Performs early exit optimizations to reduce computation
- Calculates similarity based on GC content, strandedness, entropy
**Usage**: Taxonomic classification and genomic profile comparison in bioinformatics research.

## ComparisonHeap (ComparisonHeap.java)
**Purpose**: Maintains a heap of the top N Comparisons encountered.
**Core Function**: Manages a priority queue of comparisons, keeping the N best matches efficiently.
**Key Features**:
- Uses custom WorstFirstComparator to prioritize worst comparisons
- Defensive copying to prevent direct modification of comparisons
- Automatic replacement of worst comparisons when heap is full
- Supports sorting and retrieving top comparisons
**Usage**: Efficiently tracking and managing the top N comparisons in sequence alignment or matching tasks.

## ConservationModel (ConservationModel.java)
**Purpose**: Models realistic genomic conservation variation using summed sine waves
**Core Function**: Calculates position-dependent mutation probabilities in sequence simulation
**Key Features**:
- Generates mutation rates using multiple sine waves
- Combines sine waves with configurable amplitude and period
- Supports random phase offsets for wave variation
- Provides methods to get mutation probability or determine mutation status
**Usage**: Simulate sequence mutations with position-dependent rates during genomic analysis

## Contig (Contig.java)
**Purpose**: Represents a genomic contig with sequence and metadata for bioinformatics analysis.
**Core Function**: Manages contig sequence data, provides iterator, and supports sketch generation.
**Key Features**:
- Stores contig name, ID, and base sequence
- Generates k-mer frequency profiles (tetramers, trimers, pentamers)
- Supports conversion to FASTA format
- Creates MinHash sketches for sequence comparison
**Usage**: Used in genome assembly and taxonomic classification workflows to represent individual sequence contigs.

## ContigRenamer (ContigRenamer.java)
**Purpose**: Renames contigs based on a SAM file, appending coverage and optional taxonomy information.
**Core Function**: Processes SAM alignments to rename contigs with depth and taxonomy details.
**Key Features**:
- Calculates per-contig coverage from SAM alignment data
- Optionally adds taxonomy ID to unnamed contigs
- Supports SAM file filtering for precise alignment processing
- Multithreaded processing of input files
**Usage**: Used in bioinformatics workflows to annotate genomic contigs with mapping statistics

## CrystalChamber (CrystalChamber.java)
**Purpose**: Recrystallization-based bin refinement using centroid clustering.
**Core Function**: Dissolves clusters and rebuilds them using iterative centroid assignment to find natural partitions.
**Key Features**:
- Uses deterministic random seeding for reproducible results
- Implements iterative k-means like clustering
- Performs binary cluster splitting with quality validation
- Uses Oracle for contig similarity measurements
**Usage**: Refines genomic bin clusters by finding more optimal partitions through iterative clustering

## CrystalChamber2 (CrystalChamber2.java)
**Purpose**: Enhanced recrystallization-based bin refinement using iterative centroid clustering.
**Core Function**: Performs binary splitting of contigs into clusters using Oracle similarity calculations instead of Euclidean distance.
**Key Features**:
- Uses farthest-first centroid initialization for maximum initial separation
- Implements k=2 binary splitting with assignment stability detection
- Prevents over-fragmentation by enforcing minimum split improvement thresholds
- Handles biological sequence clustering with reproducible random seeding
**Usage**: Refines contig clusters by iteratively separating them into more homogeneous groups using similarity-based assignments.

## DataLoader (DataLoader.java)
**Purpose**: Manages loading and processing of genomic contigs, reads, and associated metadata
**Core Function**: Handles data input, depth calculation, and contig loading from various file formats
**Key Features**:
- Supports loading contigs from FASTA, SAM, and BAM files
- Calculates contig depth using multiple methods (coverage stats, Bloom filter)
- Parses taxonomy IDs and contig metadata
- Handles multi-sample depth tracking
**Usage**: Primary data loader for genomic assembly and binning processes in BBTools

## EnsembleRefiner (EnsembleRefiner.java)
**Purpose**: Ensemble-based bin refinement combining multiple clustering strategies using consensus voting.
**Core Function**: Applies CrystalChamber, GraphRefiner, and EvidenceRefiner to cluster contigs, selecting results with multi-method agreement.
**Key Features**:
- Parallel refinement using three different clustering methods
- Builds co-occurrence matrix to determine robust cluster boundaries
- Applies consensus threshold to filter splitting decisions
- Tracks split attempts and successful refinements
**Usage**: Used for advanced contig clustering with multiple algorithmic perspectives, improving clustering accuracy in genomic assemblies.

## EvidenceRefiner (EvidenceRefiner.java)
**Purpose**: Evidence-based bin refinement using DBSCAN-style density clustering.
**Core Function**: Identifies dense regions of similar contigs separated by sparse boundaries, automatically determining cluster count and identifying outliers.
**Key Features**:
- Uses Oracle similarity metric for contig clustering
- Implements DBSCAN algorithm for density-based clustering
- Validates cluster quality through internal cohesion and external separation
- Handles noise points and small clusters dynamically
**Usage**: Refines genomic bin clusters by finding natural community boundaries in contig datasets.

## FileRenamer (FileRenamer.java)
**Purpose**: Renames files with top sketch hit taxonomic ID
**Core Function**: Batch renames input files by prefixing their filename with the top taxonomic hit from a RefSeq sketch query
**Key Features**:
- Multi-threaded sketch generation for input files
- Queries RefSeq database to find top taxonomic matches
- Generates new filename using taxonomic ID prefix
- Performs atomic file renaming operation
**Usage**: Used for batch renaming files based on taxonomic identification of genomic sketches

## GeneTools (GeneTools.java)
**Purpose**: Static utility class for configuring and managing gene detection models in prokaryotic and eukaryotic genomes
**Core Function**: Provides static methods to initialize, configure, and retrieve gene callers using pre-defined Markov models
**Key Features**:
- Lazy initialization of gene detection models
- Configurable detection for 16S, 18S, 5S, 23S rRNAs, tRNA, and coding sequences
- Automatic file path resolution for gene models
- Synchronized thread-safe model loading
**Usage**: Used for automated gene detection and annotation in genomic sequence analysis

## GradeBins (GradeBins.java)
**Purpose**: Grades and evaluates genomic bins for quality and taxonomic composition
**Core Function**: Processes genomic bins, calculates completeness, contamination, and taxonomic metrics
**Key Features**:
- Supports multi-threaded bin processing
- Calculates bin quality metrics like completeness and contamination
- Supports multiple quality assessment methods (CheckM, EukCC)
- Generates taxonomic classification and lineage reporting
**Usage**: Used for evaluating metagenomic assembly quality, bin classification, and taxonomic profiling

## GraphRefiner (GraphRefiner.java)
**Purpose**: Graph-based bin refinement using modularity maximization for community detection in contigs.
**Core Function**: Constructs similarity graph between contigs and applies Louvain-style community detection to identify natural cluster boundaries.
**Key Features**:
- Builds weighted similarity graph from contigs
- Applies modularity-maximizing community detection
- Iteratively moves nodes to maximize community structure
- Validates split using modularity improvement threshold
**Usage**: Refines genome binning by detecting more natural contig groupings beyond centroid methods

## GTDBLine (GTDBLine.java)
**Purpose**: Parses and manages Genome Taxonomy Database (GTDB) taxonomic lineage information
**Core Function**: Extracts and organizes taxonomic classification from delimited strings
**Key Features**:
- Parses taxonomy hierarchy into domain, phylum, class, order, family, genus, species
- Supports retrieval of taxonomic levels through hierarchical method chaining
- Handles multi-level classification parsing with prefix-based identification
- Provides method to get specific taxonomic level by numeric index
**Usage**: Processing and organizing genomic taxonomy data from tab and semicolon-delimited input

## IDComparator (IDComparator.java)
**Purpose**: Comparator for Bin objects to enable sorting based on their ID
**Core Function**: Provides a singleton comparator that compares Bin objects by their integer ID
**Key Features**:
- Implements Comparator<Bin> interface
- Uses simple integer subtraction for comparison
- Private constructor prevents multiple instances
- Provides a static final singleton comparator
**Usage**: Used for sorting collections of Bin objects in ascending order by ID

## Key (Key.java)
**Purpose**: Quantization utility for binning genomic data by GC content and depth
**Core Function**: Converts continuous genomic measurements into discrete quantization levels
**Key Features**:
- Logarithmic depth quantization with low-depth offset
- Linear GC content quantization 
- Bit-packed hash code generation
- Static parsing and configuration methods
**Usage**: Used for binning and comparing genomic sequence characteristics across different scales

## KeyValue (KeyValue.java)
**Purpose**: Represents a key-value pair that supports custom sorting for IntHashMap conversion
**Core Function**: Converts IntHashMap to a sorted ArrayList of key-value pairs
**Key Features**:
- Sorts entries by value in descending order
- Provides secondary sorting by key in ascending order
- Supports conversion from IntHashMap to list
- Implements Comparable for custom sorting
**Usage**: Used for transforming and sorting integer-based hash map data

## KmerProb (KmerProb.java)
**Purpose**: Calculates k-mer occurrence probabilities for genomic sequence analysis
**Core Function**: Provides statistical modeling of k-mer frequencies for sequence classification
**Key Features**:
- Models k-mer probability distributions
- Supports multiple k-mer length analysis
- Handles frequency normalization and scoring
- Optimized for large-scale genomic datasets
**Usage**: Core component for k-mer based sequence classification and binning workflows

## Oracle (Oracle.java)
**Purpose**: Machine learning-based similarity comparison tool for genomic binning
**Core Function**: Compares genomic bins using multi-level k-mer, depth, and network-based similarity metrics
**Key Features**:
- Uses multiple k-mer comparison techniques (3-mer, 4-mer, 5-mer)
- Applies machine learning network for advanced bin similarity assessment
- Supports depth ratio, GC content, and edge weight comparisons
- Handles taxonomic compatibility checks
**Usage**: Assists in genomic bin classification and merging during metagenomic assembly processes

## QuickBin (QuickBin.java)
**Purpose**: Prototype for metagenome contig binning.
**Core Function**: Processes genomic contigs by clustering, refining, and organizing them into taxonomic or functional groups.
**Key Features**:
- Supports multiple clustering strategies (by taxonomy, tetramer frequency)
- Configurable cluster refinement and purification
- Multi-pass edge following for bin merging
- Detailed performance and quality metrics reporting
**Usage**: Used for organizing and classifying genomic contigs in metagenomic datasets, enabling taxonomic and functional analysis.

## SamLoader (SamLoader.java)
**Purpose**: Loads and processes SAM/BAM alignment files for genomic analysis
**Core Function**: Handles SAM/BAM file parsing and extracts alignment information for binning workflows
**Key Features**:
- Efficient SAM/BAM file parsing
- Extracts alignment statistics and coverage information
- Supports multi-threaded processing
- Handles large-scale alignment datasets
**Usage**: Core component for processing alignment data in genomic binning and analysis workflows

## ScoreComparator (ScoreComparator.java)
**Purpose**: Implements a custom Comparator for Bin objects with multi-level sorting logic
**Core Function**: Provides a consistent comparison method for sorting Bin objects based on score, size, and ID
**Key Features**:
- Sorts primarily by bin score (ascending order)
- Secondary sorting by bin size 
- Final tiebreaker using bin ID
**Usage**: Used for ordering and sorting collections of Bin objects with precise comparison rules

## SimilarityMeasures (SimilarityMeasures.java)
**Purpose**: Calculates similarity metrics between numerical arrays using various statistical measures
**Core Function**: Computes distance and similarity metrics like cosine difference, Euclidean distance, Jensen-Shannon divergence
**Key Features**:
- Supports multiple distance calculation methods: cosine, Euclidean, absolute difference
- Handles float, int, and long array comparisons
- Implements advanced statistical divergence measures like Jensen-Shannon
- Provides flexible normalization and compensation techniques
**Usage**: Used for comparing frequency distributions, especially in bioinformatics and machine learning contexts

## Sketchable (Sketchable.java)
**Purpose**: Interface defining sketch generation capabilities for genomic objects
**Core Function**: Provides contract for objects that can generate MinHash sketches for similarity comparison
**Key Features**:
- Standardized sketch generation interface
- Supports multiple sketch types and parameters
- Enables polymorphic sketch creation
- Optimized for large-scale similarity analysis
**Usage**: Core interface for sketch-based genomic similarity analysis and classification workflows

## SketchRecord (SketchRecord.java)
**Purpose**: Represents a genomic sketch record with metadata for similarity analysis
**Core Function**: Stores sketch data and associated metadata for genomic sequence comparison
**Key Features**:
- Encapsulates sketch data with genomic metadata
- Supports multiple sketch formats and types
- Handles taxonomic and sequence information
- Optimized for efficient storage and retrieval
**Usage**: Core data structure for sketch-based genomic analysis and classification workflows

## SpectraCounter (SpectraCounter.java)
**Purpose**: Concurrent processing and analysis of genomic contigs with multi-threaded feature extraction
**Core Function**: Manages parallel processing of contigs, extracting genomic metrics like entropy, depth, and potential gene markers
**Key Features**:
- Concurrent contig processing with thread-based load balancing
- Optional entropy and strandedness calculation
- 16S/18S ribosomal gene detection support
- Depth and taxonomic ID parsing capabilities
**Usage**: Used for detailed genomic contig characterization in bioinformatics pipelines