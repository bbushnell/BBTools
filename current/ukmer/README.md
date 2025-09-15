# UKmer Package

This package provides specialized k-mer processing utilities with advanced data structures for efficient k-mer storage, counting, and analysis. It includes optimized hash tables, buffering systems, and multi-dimensional storage implementations.

## AbstractKmerTableU (AbstractKmerTableU.java)
**Purpose**: Abstract base class for k-mer hash table implementations in ukmer package
**Core Function**: Provides abstract methods for k-mer table operations like increment, set, and retrieve values
**Key Features**:
- Abstract methods for k-mer value manipulation
- Support for k-mer ownership tracking
- Static utility methods for k-mer text conversion
- Thread-safe allocation methods
**Usage**: Base class for different k-mer table data structures in genomic sequence analysis

## DumpThreadU (DumpThreadU.java)
**Purpose**: Multi-threaded k-mer table dumping utility for efficient serialization of k-mer data
**Core Function**: Handles concurrent dumping of k-mer tables to output streams with load balancing
**Key Features**:
- Supports parallel dumping of multiple k-mer tables
- Thread-safe output stream management
- Configurable thread count and load balancing
- Handles various output formats and compression
**Usage**: Used for serializing large k-mer datasets to files with optimal performance

## HashArrayU (HashArrayU.java)
**Purpose**: Stores kmers in a long[] and values in an int[][], with a victim cache for collision handling
**Core Function**: Multi-dimensional hash table implementation for k-mer storage with linear probing and secondary caching
**Key Features**:
- Supports variable-length k-mer storage using multi-dimensional long arrays
- Implements linear probing with configurable extra search slots
- Uses a victim cache (HashForestU) to handle hash collisions
- Supports single and multi-value k-mer storage modes
**Usage**: Used in k-mer processing and genomic data analysis for efficient k-mer storage and retrieval

## HashArrayU1D (HashArrayU1D.java)
**Purpose**: Single-dimensional hash array for k-mer storage with optimized memory layout
**Core Function**: Provides efficient k-mer storage using single-dimensional arrays for improved cache performance
**Key Features**:
- Single-dimensional array storage for reduced memory overhead
- Optimized for cache locality and memory efficiency
- Supports standard k-mer operations (increment, set, get)
- Linear probing with configurable collision handling
**Usage**: Used for memory-efficient k-mer storage in applications with limited memory requirements

## HashArrayU2D (HashArrayU2D.java)
**Purpose**: Two-dimensional hash array for k-mer storage with enhanced collision management
**Core Function**: Implements k-mer storage using two-dimensional arrays for better collision distribution
**Key Features**:
- Two-dimensional array structure for improved collision handling
- Enhanced load balancing across array dimensions
- Supports multi-value k-mer storage and retrieval
- Configurable array sizing and growth strategies
**Usage**: Used for k-mer storage in high-collision scenarios requiring better distribution

## HashArrayUHybrid (HashArrayUHybrid.java)
**Purpose**: Stores kmers in a long[] and counts in an int[], with a victim cache
**Core Function**: Hybrid hash array implementation for k-mer storage with flexible value management
**Key Features**:
- Supports single and multi-value k-mer count storage
- Implements increment and value insertion methods
- Uses IntList2 for secondary storage of multi-value arrays
- Handles hash collisions with a victim cache
**Usage**: Advanced k-mer counting and storage for genomic sequence analysis

## HashBufferU (HashBufferU.java)
**Purpose**: Multi-way buffered k-mer hash table for efficient k-mer tracking and management
**Core Function**: Distributes k-mers across multiple backend hash tables using hash-based routing
**Key Features**:
- Smart buffer flushing with force and try-lock mechanisms
- Supports count-only and value-based k-mer insertion modes
- Handles k-mer ownership and routing across multiple tables
- Efficient k-mer reconstruction and incremental storage
**Usage**: Used in k-mer processing pipelines for distributed, memory-efficient k-mer tracking

## HashForestU (HashForestU.java)
**Purpose**: Forest-based secondary storage system for handling hash table overflows
**Core Function**: Provides tree-based storage for k-mers that don't fit in primary hash tables
**Key Features**:
- Tree-based data structure for overflow k-mer storage
- Supports dynamic growth and efficient k-mer lookup
- Handles high-collision scenarios with balanced tree operations
- Integrates with primary hash tables as victim cache
**Usage**: Used as secondary storage for k-mer hash tables when primary storage reaches capacity

## HistogramMakerU (HistogramMakerU.java)
**Purpose**: Creates frequency histograms for k-mer tables with multi-threaded and single-threaded processing options
**Core Function**: Generates population frequency distribution for AbstractKmerTableU collections by counting occurrences
**Key Features**:
- Supports parallel processing for large table collections
- Dynamic thread allocation based on system thread count
- Load-balanced thread processing with atomic work allocation
- Fallback to single-threaded mode for limited thread environments
**Usage**: Generates frequency distributions for k-mer tables in genomic analysis workflows

## Kmer (Kmer.java)
**Purpose**: Core k-mer representation and manipulation utilities for genomic sequence analysis
**Core Function**: Provides k-mer creation, comparison, and conversion operations with support for various k-mer lengths
**Key Features**:
- Supports k-mer creation from sequences and conversion to strings
- Implements k-mer comparison and equality operations
- Handles canonical k-mer representation and reverse complements
- Provides k-mer encoding and decoding functionality
**Usage**: Fundamental data structure for k-mer-based genomic sequence analysis and processing

## KmerBufferU (KmerBufferU.java)
**Purpose**: Buffered k-mer storage system for efficient batch processing of k-mers
**Core Function**: Provides buffered insertion and batch processing of k-mers with automatic flushing
**Key Features**:
- Automatic buffer management with configurable flush thresholds
- Supports batch k-mer insertion for improved performance
- Handles buffer overflow and memory management
- Integrates with k-mer hash tables for efficient storage
**Usage**: Used in k-mer processing pipelines for batch insertion and memory-efficient k-mer handling

## KmerNodeU (KmerNodeU.java)
**Purpose**: Node representation for k-mer data structures with support for linked storage
**Core Function**: Provides node-based storage for k-mers in linked data structures and trees
**Key Features**:
- Supports k-mer storage with associated values and metadata
- Implements node linking and traversal operations
- Handles k-mer comparison and sorting within node structures
- Supports dynamic node allocation and memory management
**Usage**: Used in tree-based and linked k-mer data structures for complex k-mer organization

## KmerNodeU1D (KmerNodeU1D.java)
**Purpose**: Single-dimensional k-mer node for optimized linear storage structures
**Core Function**: Provides k-mer node implementation optimized for single-dimensional storage arrays
**Key Features**:
- Single-dimensional array-based k-mer node storage
- Optimized for linear access patterns and cache efficiency
- Supports basic k-mer operations with minimal overhead
- Handles k-mer value storage and retrieval
**Usage**: Used in linear k-mer storage structures requiring minimal memory overhead

## KmerNodeU2D (KmerNodeU2D.java)
**Purpose**: Two-dimensional k-mer node for enhanced storage organization
**Core Function**: Provides k-mer node implementation with two-dimensional organization for better data locality
**Key Features**:
- Two-dimensional node organization for improved access patterns
- Enhanced k-mer storage with multi-dimensional indexing
- Supports complex k-mer relationships and hierarchical storage
- Optimized for two-dimensional k-mer access patterns
**Usage**: Used in advanced k-mer data structures requiring two-dimensional organization

## KmerTableSetU (KmerTableSetU.java)
**Purpose**: Collection manager for multiple k-mer tables with unified operations
**Core Function**: Manages collections of k-mer tables and provides unified access to distributed k-mer data
**Key Features**:
- Manages multiple AbstractKmerTableU instances as a unified collection
- Provides load balancing and distribution across table instances
- Supports parallel operations across multiple k-mer tables
- Handles table creation, sizing, and memory management
**Usage**: Used for managing large-scale k-mer datasets distributed across multiple hash tables

## OwnershipThread (OwnershipThread.java)
**Purpose**: Thread management utility for k-mer table ownership and concurrent access control
**Core Function**: Manages thread ownership and synchronization for concurrent k-mer table operations
**Key Features**:
- Provides thread-safe ownership tracking for k-mer tables
- Handles concurrent access coordination and lock management
- Supports thread-local storage and ownership transfer
- Manages resource allocation and cleanup for threaded operations
**Usage**: Used in multi-threaded k-mer processing for safe concurrent table access

## WalkerU (WalkerU.java)
**Purpose**: Iterator utility for traversing k-mer tables and extracting k-mer data
**Core Function**: Provides iteration and traversal functionality for k-mer hash tables and data structures
**Key Features**:
- Supports sequential and random access k-mer traversal
- Handles k-mer extraction and enumeration from hash tables
- Provides filtering and selection capabilities during traversal
- Optimized for large-scale k-mer table iteration
**Usage**: Used for k-mer table traversal, data extraction, and batch processing operations