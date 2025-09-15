# kmer Package - K-mer Processing and Storage Infrastructure
*High-performance k-mer counting, storage, and analysis tools for genomic sequence processing*

## Package Overview
The `kmer` package provides comprehensive infrastructure for k-mer (k-length nucleotide subsequences) processing, including efficient hash tables, counting algorithms, and parallel processing capabilities essential for genomic analysis.

---

## Abstract Base Classes

### AbstractKmerTable (AbstractKmerTable.java)
**Purpose**: Abstract base class for k-mer counting and storage strategies in bioinformatics
**Core Function**: Provides uniform interface for k-mer table implementations with flexible storage methods
**Key Features**:
  - Generic k-mer increment/retrieval methods
  - Thread-safe ownership tracking
  - Multi-format k-mer export (FASTA/text)
  - Memory-safe allocation utilities
**Usage**: Base class for specialized k-mer table implementations in genomic analysis

### AbstractKmerTableSet (AbstractKmerTableSet.java)
**Purpose**: Manages k-mer loading and processing for Tadpole's computational workflows
**Core Function**: Handles multi-threaded k-mer table construction from input files
**Key Features**:
  - Supports prefiltering of low-depth k-mers
  - Generates k-mer count histograms
  - Tracks input reads, bases, and loaded k-mers
**Usage**: Base class for k-mer table implementations in bioinformatics sequence analysis

---

## Hash Table Implementations

### HashArray (HashArray.java)
**Purpose**: Abstract k-mer storage implementation with linear probing and victim cache
**Core Function**: Efficient k-mer hash table with automatic resizing and thread-safe access
**Key Features**:
  - Linear probing for collision resolution
  - Automatic resizing with configurable load factors
  - Thread-safe ownership and insertion mechanisms
  - Supports both 1D and 2D value storage
**Usage**: Base class for specialized k-mer hash table implementations in genomic data processing

### HashArray1D (HashArray1D.java)
**Purpose**: Stores k-mers in a long[] with integer counts using linear probing and victim cache
**Core Function**: Manages k-mer hash table with dynamic resizing and overflow handling
**Key Features**:
  - Linear probing for k-mer insertion with victim cache fallback
  - Automatic resizing when load factor exceeded
  - Integer overflow protection by capping at MAX_VALUE
**Usage**: Efficient k-mer counting in bioinformatics sequence analysis

### HashArray2D (HashArray2D.java)
**Purpose**: Stores k-mers in a long[] with values in a 2D int[][] array, supporting efficient value management
**Core Function**: Manages k-mer hash storage with automatic array resizing and value insertion
**Key Features**:
  - Supports multiple values per k-mer using 2D integer arrays
  - Handles dynamic memory growth with prime number scheduling
  - Prevents duplicate value insertions
**Usage**: Specialized k-mer hash table for storing and retrieving k-mer sequences with associated values

### HashArrayHybrid (HashArrayHybrid.java)
**Purpose**: Stores k-mers in a long[] and counts in an int[], with a victim cache for efficient storage.
**Core Function**: Hybrid hash array that supports incremental k-mer count tracking with overflow handling
**Key Features**:
  - Automatic resizing when storage limit reached
  - Supports single and multiple value insertions
  - Efficient memory management with victim cache
**Usage**: K-mer frequency tracking in bioinformatics data processing

### HashArrayHybridFast (HashArrayHybridFast.java)
**Purpose**: Stores k-mers in a hybrid hash array with efficient count tracking and victim cache
**Core Function**: Incrementally store and count k-mer occurrences using a primary hash table and overflow mechanism
**Key Features**:
  - Linear probing with extra space for fast primary table lookups
  - Integrated victim cache for handling hash collisions
  - Automatic resizing when size limits are exceeded
**Usage**: Efficient k-mer counting in bioinformatics data processing

### HashForest (HashForest.java)
**Purpose**: Specialized k-mer hash table with tree-based storage and dynamic resizing
**Core Function**: Implements a prime-sized hash array of binary search trees for k-mer frequency tracking
**Key Features**:
  - Supports 1D and 2D k-mer node storage
  - Automatic resizing with prime number capacity
  - Efficient incremental updates and retrievals
**Usage**: Tracks k-mer frequencies in bioinformatics sequence analysis

---

## Node and Link Structures

### KmerNode (KmerNode.java)
**Purpose**: Abstract base class for k-mer navigation and storage using a binary search tree structure.
**Core Function**: Manages k-mer counts with binary tree traversal and incremental counting
**Key Features**:
  - Binary tree-based k-mer storage with pivot-based navigation
  - Supports incremental counting with overflow protection
  - Abstract methods for node creation and value management
**Usage**: Used in k-mer frequency tracking and genome analysis algorithms

### KmerNode1D (KmerNode1D.java)
**Purpose**: One-dimensional k-mer node implementation for storing a single integer value associated with a k-mer pivot.
**Core Function**: Represents a single-value k-mer node in a binary tree structure
**Key Features**:
  - Supports single integer value storage
  - Implements k-mer dumping methods for bytes and text
  - Part of k-mer data structure for bioinformatics processing
**Usage**: Used in k-mer frequency analysis and genomic data processing

### KmerNode2D (KmerNode2D.java)
**Purpose**: Enables multiple values per k-mer in a specialized node structure
**Core Function**: Two-dimensional k-mer node supporting multiple integer values per node
**Key Features**:
  - Dynamically resizable value storage
  - Supports inserting multiple values per k-mer
  - Implements binary tree-like node structure
**Usage**: Used in k-mer frequency and sequence analysis algorithms

### KmerLink (KmerLink.java)
**Purpose**: Linked-list implementation for storing k-mer counts with dynamic chaining mechanism
**Core Function**: Manages k-mer count tracking using a recursive linked chain
**Key Features**:
  - Increments k-mer counts recursively
  - Supports dynamic node insertion
  - Provides thread-safe ownership tracking
**Usage**: Used in k-mer frequency analysis and counting in genomic sequence processing

---

## Table Management and Processing

### KmerTable (KmerTable.java)
**Purpose**: Hash table implementation for storing and managing k-mer counts with automatic resizing
**Core Function**: Provides efficient storage and retrieval of k-mer frequencies using linked list hash table
**Key Features**:
  - Automatic hash table resizing when load factor exceeded
  - Thread-safe k-mer count tracking
  - Supports incremental counting with integer overflow protection
**Usage**: Tracking and analyzing k-mer frequencies in genomic sequence analysis

### KmerTableSet (KmerTableSet.java)
**Purpose**: Manages k-mer tables for processing DNA sequence data in Tadpole algorithm
**Core Function**: Loads, stores, and manipulates k-mers across multiple hash tables
**Key Features**:
  - Handles forward and reverse complement k-mer storage
  - Supports multi-threaded k-mer loading
  - Enables quality trimming and filtering of reads
  - Provides methods for k-mer counting and ownership tracking
**Usage**: Used in DNA sequence analysis for tasks like assembly and error correction

### SimpleKmerTable (SimpleKmerTable.java)
**Purpose**: Lightweight k-mer table for storing and managing k-mer sequences
**Core Function**: Provides basic infrastructure for k-mer sequence storage
**Key Features**:
  - Supports custom initial table size via constructor
  - Minimal implementation, likely intended for extension
**Usage**: Base class for k-mer sequence management in genomic analysis

---

## Buffer and Utility Classes

### KmerBuffer (KmerBuffer.java)
**Purpose**: Efficient buffer for storing k-mer sequences with optional associated values
**Core Function**: Manages a dynamic list of long-encoded k-mers with optional integer values
**Key Features**:
  - Supports adding single or multiple k-mers
  - Optional value tracking alongside k-mers
  - Constant-time clear and size operations
**Usage**: Temporary storage during k-mer processing and analysis in bioinformatics algorithms

### HashBuffer (HashBuffer.java)
**Purpose**: Buffered k-mer hash table with multi-way distribution for efficient k-mer tracking
**Core Function**: Routes k-mers across multiple tables using consistent hash distribution
**Key Features**:
  - Distributes k-mers across multiple AbstractKmerTables
  - Supports incremental k-mer addition with buffer management
  - Handles k-mer routing via `kmerToWay()` method
**Usage**: Intermediate storage and routing for k-mer counting in genomic analysis

### Walker (Walker.java)
**Purpose**: Abstract base class for iterating through k-mer hash map data structures
**Core Function**: Provides abstract methods for traversing k-mer collections
**Key Features**:
  - `next()`: Advances iterator
  - `kmer()`: Retrieves current k-mer key
  - `value()`: Retrieves current value
**Usage**: Used as base for k-mer iteration strategies in BBTools k-mer processing packages

---

## Analysis and Processing Tools

### HistogramMaker (HistogramMaker.java)
**Purpose**: Creates k-mer frequency histograms from hash tables with multi-threaded and single-threaded processing strategies.
**Core Function**: Generates frequency distribution of k-mer occurrences across multiple hash tables
**Key Features**:
  - Automatic thread selection based on CPU cores
  - Single-threaded processing for small datasets
  - Multi-threaded processing for high-performance computing
**Usage**: Analyze k-mer frequencies in genomic sequence data processing

### TableReader (TableReader.java)
**Purpose**: Reads and processes k-mer tables for sequence matching and masking.
**Core Function**: Performs k-mer based read processing, including matching, counting, and masking
**Key Features**:
  - Supports forward and reverse complement k-mer matching
  - Configurable k-mer length and Hamming distance search
  - Finds best sequence matches based on k-mer hits
**Usage**: Used in bioinformatics for sequence alignment and filtering reads

### TableLoaderLockFree (TableLoaderLockFree.java)
**Purpose**: Parallel k-mer table loader for genomic reference sequences with thread-safe, lock-free processing
**Core Function**: Load k-mers from reference sequences into concurrent hash tables using multiple threads
**Key Features**:
  - Supports variable-length k-mer processing
  - Handles reverse complement and canonical k-mer representation
  - Supports error correction and k-mer mutations
  - Configurable k-mer storage modes (set/increment)
**Usage**: Load and index genomic k-mers for efficient sequence comparison and analysis

---

## Threading and Coordination

### DumpThread (DumpThread.java)
**Purpose**: Multi-threaded worker for efficiently dumping k-mer tables across multiple processing threads.
**Core Function**: Coordinates parallel k-mer table dumping with atomic work distribution
**Key Features**:
  - Dynamically allocates threads based on table count
  - Supports minimum and maximum k-mer count filtering
  - Thread-safe output using synchronized ByteStreamWriter
**Usage**: Parallel processing of large k-mer tables during genomic analysis

### OwnershipThread (OwnershipThread.java)
**Purpose**: Parallel thread for initializing or clearing ownership of AbstractKmerTable instances
**Core Function**: Processes AbstractKmerTable ownership using multi-threaded approach
**Key Features**:
  - Supports concurrent initialization and clearing of table ownership
  - Uses AtomicInteger for thread-safe index tracking
  - Dynamically adjusts thread count based on available resources
**Usage**: Manages concurrent ownership operations for k-mer tables in parallel processing scenarios

### ScheduleMaker (ScheduleMaker.java)
**Purpose**: Calculates optimal memory allocation and resizing strategy for k-mer hash tables
**Core Function**: Generates prime-sized hash table schedule with intelligent memory management
**Key Features**:
  - Dynamically calculates hash table sizes based on available memory
  - Converts table sizes to prime numbers for optimal hash distribution
  - Supports multiple hash table shards for parallel processing
**Usage**: Memory-efficient k-mer hash table initialization in computational genomics

---

## Package Usage

The kmer package provides the foundational infrastructure for:

- **K-mer Counting**: Efficient storage and tracking of k-mer frequencies
- **Memory Management**: Intelligent hash table sizing and memory allocation
- **Parallel Processing**: Multi-threaded k-mer loading and analysis
- **Data Structures**: Optimized hash tables, trees, and linked structures
- **Genomic Analysis**: Support for forward/reverse complement processing
- **Performance Optimization**: Multiple storage strategies for different use cases

---
*Documentation generated using evidence-based analysis of source code*