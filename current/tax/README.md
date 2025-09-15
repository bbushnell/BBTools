# Tax Package

This package provides comprehensive taxonomic infrastructure for biological sequence classification, identifier mapping, and phylogenetic analysis. It includes tools for NCBI taxonomy integration, accession number processing, and hierarchical tree management.

## AccessionToTaxid (AccessionToTaxid.java)
**Purpose**: Tool for mapping sequence accession numbers to taxonomy IDs
**Core Function**: Processes input files to build a mapping between sequence accession numbers and NCBI taxonomy identifiers
**Key Features**:
- Supports parsing input files with accession and taxid information
- Uses hash-based and map-based indexing for efficient lookup
- Handles multiple input file formats and compression types
- Provides methods to validate and extract accession numbers
**Usage**: Used in bioinformatics pipelines to map sequence identifiers to taxonomic classifications

## AccessionToTaxid_old (AccessionToTaxid_old.java)
**Purpose**: Legacy tool for mapping sequence accession numbers to taxonomy IDs
**Core Function**: Processes input files to build a mapping between sequence accession numbers and NCBI taxonomy identifiers
**Key Features**:
- Supports parsing input files with accession and taxid information
- Uses hash-based and map-based indexing for efficient lookup
- Handles multiple input file formats and compression types
- Provides methods to validate and extract accession numbers
**Usage**: Used in bioinformatics pipelines to map sequence identifiers to taxonomic classifications

## AnalyzeAccession (AnalyzeAccession.java)
**Purpose**: Counts patterns in Accessions and handles hashing for Accession to TaxID lookups
**Core Function**: Analyzes input files to extract and count unique accession number patterns through multi-threaded processing
**Key Features**:
- Supports parallel processing of multiple input files
- Generates combinatorial complexity metrics for accession patterns
- Maps characters to abstract patterns (digits, letters, delimiters)
- Encodes accession numbers with pattern-based compression
**Usage**: Used for analyzing and encoding biological sequence accession numbers, particularly in taxonomic or genomic datasets

## AnalyzeAccession_ST (AnalyzeAccession_ST.java)
**Purpose**: Single-threaded version for counting patterns in Accessions and handling hashing for Accession to TaxID lookups
**Core Function**: Processes input files to analyze and count unique accession patterns using character remapping
**Key Features**:
- Pattern detection for accession codes (digits, letters)
- Supports multiple input file processing
- Generates pattern count statistics
- Calculates potential pattern combinations
**Usage**: Analyzes taxonomic accession codes, counting and classifying their structural patterns

## CanonicalLineage (CanonicalLineage.java)
**Purpose**: Represents canonical taxonomic lineage with ordered hierarchical levels
**Core Function**: Manages taxonomic hierarchy from domain to species level with standardized lineage representation
**Key Features**:
- Maintains ordered taxonomic levels (domain, phylum, class, order, family, genus, species)
- Supports lineage comparison and validation
- Handles missing or incomplete taxonomic information
- Provides standardized string representation of lineages
**Usage**: Used for taxonomic classification and lineage comparison in phylogenetic analysis

## ExplodeTree (ExplodeTree.java)
**Purpose**: Constructs a directory and file tree of sequences corresponding to a taxonomic tree
**Core Function**: Processes input FASTA files, categorizes sequences by taxonomic node, and writes them to corresponding directory structures
**Key Features**:
- Creates directory tree based on taxonomic hierarchy
- Writes name files for each taxonomic node
- Processes FASTA files and tracks read/line/base counts
- Supports optional results file generation
**Usage**: Organizing genomic sequences into a structured taxonomy-based file system

## FilterByTaxa (FilterByTaxa.java)
**Purpose**: Filters sequences according to their taxonomy, as determined by the sequence name
**Core Function**: Processes input read files and filters reads based on taxonomic inclusion/exclusion criteria
**Key Features**:
- Supports filtering sequences labeled with gi number or NCBI taxID
- Can process single and paired-end read files
- Supports best-effort taxonomy matching
- Generates optional results file with matched taxonomy nodes
**Usage**: Used to filter genomic sequencing data by specific taxonomic constraints

## FindAncestor (FindAncestor.java)
**Purpose**: Finds and classifies taxonomy ancestors for given identifiers in biological data
**Core Function**: Processes input lines with GI numbers, converts to taxonomy IDs, and determines common ancestors
**Key Features**:
- Converts GI numbers to NCBI taxonomic IDs
- Finds lowest common ancestor for a set of taxa
- Computes majority vote taxonomy classification
- Supports extensive taxonomy tree traversal
**Usage**: Taxonomic classification and ancestral analysis in biological sequence datasets

## GiToTaxid (GiToTaxid.java)
**Purpose**: Utility for mapping GenInfo (GI) numbers to NCBI Taxonomy IDs
**Core Function**: Parses and converts biological sequence identifiers to taxonomy numbers
**Key Features**:
- Supports parsing GI numbers from multiple sequence formats (EMBL, GenBank)
- Handles both string and byte array input parsing
- Validates and extracts taxonomy IDs from complex sequence headers
- Provides thread-safe initialization of GI-to-Taxid mapping tables
**Usage**: Used in taxonomic classification of biological sequences by converting identifier numbers

## GiToTaxidInt (GiToTaxidInt.java)
**Purpose**: Utility class for converting GenBank identifier (GI) numbers to NCBI taxonomy IDs
**Core Function**: Parses and maps GI numbers or NCBI taxonomy identifiers to corresponding taxonomy IDs
**Key Features**:
- Supports parsing GI numbers from multiple string/byte formats
- Handles various delimiter types ('|', '~', '_')
- Can initialize mapping arrays from files
- Robust error handling for invalid inputs
**Usage**: Used in bioinformatics workflows to translate sequence identifiers to taxonomic classifications

## IDNode (IDNode.java)
**Purpose**: Support class for creating hierarchical ID trees with maximum value tracking
**Core Function**: Builds hierarchical tree structure using priority queue and node merging algorithm
**Key Features**:
- Creates tree from array of nodes using priority-based merging
- Implements custom compareTo for node priority
- Supports Newick tree format generation
- Tracks node positions and maximum values using BitSet
**Usage**: Used in taxonomic tree construction and hierarchical data representation

## IDTree (IDTree.java)
**Purpose**: Hierarchical tree structure for managing taxonomic identifiers and relationships
**Core Function**: Builds and maintains tree-based data structure for efficient taxonomic lookup and traversal
**Key Features**:
- Supports tree construction from taxonomic data
- Provides efficient node lookup and traversal algorithms
- Handles hierarchical relationships between taxonomic levels
- Optimized for large-scale taxonomic databases
**Usage**: Core data structure for taxonomic tree management and phylogenetic analysis

## ImgRecord (ImgRecord.java)
**Purpose**: Represents and manages IMG (Integrated Microbial Genomes) record data
**Core Function**: Stores and processes microbial genome metadata from IMG database entries
**Key Features**:
- Stores IMG identifier, taxonomy information, and genome metadata
- Supports parsing and validation of IMG record formats
- Handles genome classification and annotation data
- Provides access methods for microbial genome properties
**Usage**: Used for processing and organizing microbial genome metadata from IMG database

## ImgRecord2 (ImgRecord2.java)
**Purpose**: Enhanced version for representing and parsing IMG (Integrated Microbial Genomes) record identifiers from text files
**Core Function**: Converts IMG record text files into HashMap or array of ImgRecord2 objects
**Key Features**:
- Parses tab-delimited IMG record files with image, taxonomy, and name information
- Converts records to HashMap for efficient lookup by image ID
- Handles parsing of IMG header formats with flexible delimiter detection
- Optional name storage controlled by static boolean flag
**Usage**: Used for processing and organizing microbial genome metadata from IMG file formats

## Lineage (Lineage.java)
**Purpose**: Represents complete taxonomic lineage information with hierarchical structure
**Core Function**: Manages taxonomic classification hierarchy from kingdom to species level
**Key Features**:
- Stores complete taxonomic path (kingdom, phylum, class, order, family, genus, species)
- Supports lineage comparison and similarity calculations
- Handles partial and complete lineage information
- Provides methods for lineage validation and formatting
**Usage**: Used for taxonomic classification, lineage analysis, and phylogenetic studies

## PrintTaxonomy (PrintTaxonomy.java)
**Purpose**: Utility for displaying and formatting taxonomic hierarchy information
**Core Function**: Processes and prints taxonomic trees and lineage information in various formats
**Key Features**:
- Supports multiple output formats for taxonomic data
- Handles tree traversal and hierarchical display
- Provides formatted output for taxonomic lineages
- Supports filtering and selection of taxonomic levels
**Usage**: Used for visualizing and reporting taxonomic classifications and tree structures

## Query (Query.java)
**Purpose**: Provides query interface for taxonomic database searches and lookups
**Core Function**: Handles taxonomic database queries and identifier resolution
**Key Features**:
- Supports multiple query types for taxonomic data
- Handles identifier lookup and conversion
- Provides efficient database search algorithms
- Supports batch query processing
**Usage**: Core query engine for taxonomic database operations and identifier resolution

## RenameGiToTaxid (RenameGiToTaxid.java)
**Purpose**: Renames sequence files by converting GI numbers to taxonomy IDs in filenames
**Core Function**: Processes files to replace GI-based naming with taxonomy-based naming conventions
**Key Features**:
- Batch processing of file renaming operations
- Converts GI numbers to corresponding taxonomy IDs
- Supports various filename formats and conventions
- Handles error checking and validation during renaming
**Usage**: Used for organizing sequence files using taxonomic identifiers instead of GI numbers

## RenameIMG (RenameIMG.java)
**Purpose**: Renames files using IMG (Integrated Microbial Genomes) identifier conventions
**Core Function**: Processes and renames files according to IMG database naming standards
**Key Features**:
- Supports IMG-specific file naming conventions
- Handles batch renaming operations
- Validates IMG identifiers during renaming process
- Provides error handling for invalid identifiers
**Usage**: Used for organizing microbial genome files according to IMG database standards

## ShrinkAccession (ShrinkAccession.java)
**Purpose**: Compresses and optimizes accession number storage and processing
**Core Function**: Reduces memory footprint of accession number databases through compression techniques
**Key Features**:
- Implements compression algorithms for accession databases
- Reduces memory usage for large-scale accession processing
- Maintains lookup efficiency while reducing storage requirements
- Supports various compression strategies
**Usage**: Used for optimizing large-scale accession number databases and improving memory efficiency

## SplitByTaxa (SplitByTaxa.java)
**Purpose**: Splits sequence datasets based on taxonomic classification
**Core Function**: Processes input sequences and divides them into separate files based on taxonomic groups
**Key Features**:
- Supports taxonomic-based file splitting operations
- Handles multiple taxonomic levels for splitting criteria
- Processes large sequence datasets efficiently
- Maintains sequence integrity during splitting process
**Usage**: Used for organizing sequence datasets by taxonomic classification for downstream analysis

## TaxApp (TaxApp.java)
**Purpose**: Main application interface for taxonomic analysis and processing workflows
**Core Function**: Provides command-line interface for executing taxonomic analysis pipelines
**Key Features**:
- Integrates multiple taxonomic processing tools
- Supports various input/output formats
- Provides workflow management for taxonomic analysis
- Handles parameter configuration and validation
**Usage**: Primary entry point for taxonomic analysis workflows and batch processing operations

## TaxClient (TaxClient.java)
**Purpose**: Client interface for taxonomic database services and remote queries
**Core Function**: Handles communication with taxonomic database servers and remote services
**Key Features**:
- Supports remote taxonomic database queries
- Handles network communication and error handling
- Provides caching mechanisms for efficient queries
- Supports various database connection protocols
**Usage**: Used for accessing remote taxonomic databases and distributed taxonomic services

## TaxFilter (TaxFilter.java)
**Purpose**: Filters and processes taxonomic data based on specified criteria
**Core Function**: Applies filtering rules to taxonomic datasets based on user-defined parameters
**Key Features**:
- Supports multiple filtering criteria for taxonomic data
- Handles complex filtering logic and rule combinations
- Processes large taxonomic datasets efficiently
- Provides validation and error checking for filter parameters
**Usage**: Used for selective processing and analysis of taxonomic datasets

## TaxNode (TaxNode.java)
**Purpose**: Represents individual nodes in taxonomic tree structures
**Core Function**: Stores taxonomic information and maintains relationships between taxonomic levels
**Key Features**:
- Stores taxonomic ID, name, rank, and parent relationships
- Supports tree navigation and traversal operations
- Handles taxonomic hierarchy validation
- Provides efficient node lookup and comparison methods
**Usage**: Core data structure for building and managing taxonomic trees and hierarchies

## TaxServer (TaxServer.java)
**Purpose**: Server component for hosting taxonomic database services
**Core Function**: Provides server infrastructure for taxonomic database queries and services
**Key Features**:
- Handles multiple concurrent client connections
- Provides taxonomic database query services
- Supports various query types and response formats
- Implements caching and optimization for improved performance
**Usage**: Used for hosting centralized taxonomic database services and supporting distributed queries

## TaxSize (TaxSize.java)
**Purpose**: Calculates and manages size metrics for taxonomic datasets
**Core Function**: Computes statistics and size information for taxonomic databases and trees
**Key Features**:
- Calculates taxonomic tree size and depth statistics
- Provides memory usage estimates for taxonomic data
- Supports performance optimization based on size metrics
- Handles large-scale taxonomic dataset analysis
**Usage**: Used for taxonomic database optimization and performance analysis

## TaxTree (TaxTree.java)
**Purpose**: Core implementation of taxonomic tree data structure and operations
**Core Function**: Manages complete taxonomic hierarchy with efficient tree operations and queries
**Key Features**:
- Implements full taxonomic tree with NCBI taxonomy support
- Provides efficient ancestor lookup and lineage computation
- Supports tree traversal and node relationship queries
- Handles dynamic tree updates and modifications
**Usage**: Primary data structure for taxonomic classification, lineage analysis, and phylogenetic operations