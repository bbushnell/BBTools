# Bloom Package - Bloom Filters and K-mer Counting

## Overview
The `bloom` package provides comprehensive implementations of Bloom filters and k-mer counting data structures for efficient sequence analysis. Bloom filters are probabilistic data structures that test set membership with space-efficient storage and fast queries, while k-mer counting provides exact frequency tracking for substrings in biological sequences. This package includes multi-threaded implementations, error correction algorithms, and specialized data structures optimized for genomic applications with large datasets.

## Key Components

### Core Data Structures

#### BloomFilter (BloomFilter.java)
**Purpose**: Probabilistic data structure for efficient set membership testing
- **Core Function**: Wraps KCountArray implementations to provide Bloom filter functionality with multithreaded reference loading
- **Key Features**:
  - Configurable hash functions and bit arrays for optimal performance
  - Memory-efficient storage with adjustable false positive rates
  - Canonical k-mer representation for consistent queries
  - Multithreaded reference sequence loading and processing
  - Serialization support for persistent storage
  - Parameter optimization for different data types
  - Integration with k-mer counting systems
- **Usage**: Fast membership testing for large k-mer sets and sequence filtering

#### KCountArray (KCountArray.java)
**Purpose**: Abstract base class for k-mer counting data structures
- **Core Function**: Provides interface for exact and approximate k-mer frequency counting
- **Key Features**:
  - Multiple implementation variants (KCountArray2-8MT) for different use cases
  - Configurable cell size and hash functions
  - Thread-safe operations with atomic counters
  - Memory-efficient packing of count values
  - Canonical k-mer support for strand-independent counting
  - Prefiltering capabilities for improved performance
  - Dynamic array sizing based on memory constraints
- **Usage**: Foundation for all k-mer counting operations

#### KCountArray7MTA (KCountArray7MTA.java)
**Purpose**: Multi-threaded k-mer counting array with atomic operations
- **Core Function**: High-performance k-mer counting with thread-safe increment operations
- **Key Features**:
  - Atomic integer arrays for thread-safe counting
  - Optimized memory layout for cache efficiency
  - Configurable bit-packing for different count ranges
  - Prefilter integration for two-stage counting
  - Parallel processing support for large datasets
  - Memory usage optimization and monitoring
- **Usage**: High-throughput k-mer counting in multithreaded environments

### Error Correction

#### ErrorCorrect (ErrorCorrect.java)
**Purpose**: Sequence error correction using k-mer frequency information
- **Core Function**: Identifies and corrects sequencing errors based on k-mer abundance patterns
- **Key Features**:
  - Multi-pass error correction algorithm
  - Configurable error detection thresholds
  - Quality score integration for correction decisions
  - Parallel processing for improved performance
  - Canonical k-mer analysis for consistent results
  - Stepwise refinement with multiple passes
  - Integration with Bloom filters for memory efficiency
- **Usage**: Preprocessing step for improving sequence quality before analysis

#### BloomFilterCorrector (BloomFilterCorrector.java)
**Purpose**: Bloom filter-based error correction system
- **Core Function**: Uses Bloom filters to efficiently identify and correct sequencing errors
- **Key Features**:
  - Memory-efficient error detection using Bloom filters
  - Fast k-mer lookup for error identification
  - Configurable correction parameters
  - Integration with quality scores
  - Batch processing capabilities
- **Usage**: Memory-efficient error correction for large datasets

### Specialized Counting

#### KmerCount3-6 (KmerCount3.java, KmerCount4.java, etc.)
**Purpose**: Specialized k-mer counting implementations for different k-mer sizes
- **Core Function**: Optimized counting for specific k-mer lengths with direct encoding
- **Key Features**:
  - Direct k-mer encoding without hashing for small k
  - Optimized data structures for specific k-mer lengths
  - Memory-efficient storage for exact counts
  - Fast lookup and increment operations
  - Integration with larger counting systems
- **Usage**: Optimized counting for small k-mer sizes (k=3 to k=6)

#### LargeKmerCount (LargeKmerCount.java)
**Purpose**: K-mer counting for large k-mer sizes using hash-based storage
- **Core Function**: Handles k-mer counting for k>31 using hash tables
- **Key Features**:
  - Hash-based storage for large k-mers
  - Dynamic memory allocation
  - Collision handling and resolution
  - Memory usage optimization
  - Integration with standard counting interfaces
- **Usage**: K-mer counting for large k-mer sizes where direct encoding is infeasible

### Filtering and Processing

#### PolyFilter (PolyFilter.java)
**Purpose**: Filters sequences based on polynucleotide content
- **Core Function**: Identifies and removes low-complexity sequences with repetitive content
- **Key Features**:
  - Polynucleotide run detection
  - Configurable filtering thresholds
  - Quality-aware filtering decisions
  - Parallel processing support
  - Integration with other filtering systems
- **Usage**: Quality control preprocessing to remove low-complexity sequences

#### ReadCounter (ReadCounter.java)
**Purpose**: Counts reads and analyzes sequence statistics
- **Core Function**: Provides comprehensive read counting and statistical analysis
- **Key Features**:
  - Read counting with filtering options
  - Length distribution analysis
  - Quality score statistics
  - Parallel processing for large datasets
  - Integration with other analysis tools
- **Usage**: Quality control and dataset characterization

### Wrapper Classes

#### BloomFilterWrapper (BloomFilterWrapper.java)
**Purpose**: High-level interface for Bloom filter operations
- **Core Function**: Provides simplified interface for common Bloom filter operations
- **Key Features**:
  - Streamlined API for common use cases
  - Parameter optimization and tuning
  - Error handling and validation
  - Integration with file I/O systems
  - Batch processing capabilities
- **Usage**: Simplified access to Bloom filter functionality

#### BloomFilterCorrectorWrapper (BloomFilterCorrectorWrapper.java)
**Purpose**: Wrapper for error correction using Bloom filters
- **Core Function**: Provides high-level interface for Bloom filter-based error correction
- **Key Features**:
  - Simplified error correction workflow
  - Parameter management and optimization
  - Integration with quality control systems
  - Batch processing support
  - Progress monitoring and reporting
- **Usage**: Streamlined error correction for production workflows

## Architecture

### Design Patterns
- **Template Method**: KCountArray provides common interface with specialized implementations
- **Strategy Pattern**: Multiple counting strategies for different data characteristics
- **Factory Pattern**: Automatic selection of optimal implementation based on parameters
- **Decorator Pattern**: Wrapper classes enhance functionality of core components
- **Observer Pattern**: Progress monitoring and statistics collection
- **Builder Pattern**: Complex parameter configuration for optimal performance

### Dependencies
- `dna.AminoAcid` - DNA sequence encoding and manipulation
- `shared.Tools` - Utility functions and data manipulation
- `structures.*` - Efficient data structures for storage and processing
- `stream.*` - Sequence I/O and processing streams
- `fileIO.*` - File handling and serialization
- `tracker.*` - Statistics collection and monitoring

### Memory Management
- **Bit Packing**: Efficient storage of count values using configurable bit widths
- **Array Segmentation**: Large arrays split into segments for improved cache performance
- **Memory Monitoring**: Dynamic memory usage tracking and optimization
- **Garbage Collection**: Minimized object allocation for improved performance

## Common Usage Examples

### Basic Bloom Filter
```java
// Create Bloom filter with specified parameters
String[] args = {"k=31", "bits=2", "hashes=3", "ref=reference.fasta"};
BloomFilter bf = new BloomFilter(args);

// Test membership
boolean contains = bf.contains("ATCGATCGATCGATCGATCGATCGATCGATCG");
```

### K-mer Counting
```java
// Create k-mer counter
long cells = 1000000000L;  // 1 billion cells
int cellBits = 8;          // 8 bits per cell (max count 255)
int hashes = 1;            // Single hash function
KCountArray kca = KCountArray.makeNew(cells, cellBits, hashes);

// Count k-mers in sequence
String sequence = "ATCGATCGATCGATCG";
int k = 21;
for (int i = 0; i <= sequence.length() - k; i++) {
    String kmer = sequence.substring(i, i + k);
    long encoded = encodeKmer(kmer);
    kca.increment(encoded);
}
```

### Error Correction
```java
// Set up error correction parameters
String[] args = {
    "in=reads.fastq",
    "out=corrected.fastq",
    "k=23",
    "cbits=4",
    "thresh1=1",
    "thresh2=2",
    "maxerrors=3"
};

// Run error correction
ErrorCorrect ec = new ErrorCorrect();
ec.process(args);
```

### Large K-mer Counting
```java
// For k-mers larger than 31
LargeKmerCount lkc = new LargeKmerCount();
lkc.setK(51);  // 51-mer counting
lkc.loadReads("input.fastq");
lkc.process();

// Get count for specific k-mer
long count = lkc.getCount("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG");
```

### Filtering Low-Complexity Sequences
```java
// Filter polynucleotide sequences
PolyFilter pf = new PolyFilter();
pf.setParameters("minLength=10", "maxRatio=0.8");
boolean isLowComplexity = pf.isPolynucleotide("AAAAAAAAAAAAAAAA");
```

## Performance Considerations
- **Memory Usage**: Choose appropriate cell sizes and hash functions based on available memory
- **Thread Safety**: Use atomic operations for multithreaded environments
- **Cache Efficiency**: Segment large arrays for better cache performance
- **Hash Functions**: Balance between speed and collision resistance
- **Bit Packing**: Optimize storage density based on expected count ranges
- **Prefiltering**: Use two-stage counting for improved memory efficiency

## Error Handling
- **Parameter Validation**: Comprehensive validation of input parameters
- **Memory Limits**: Graceful handling of memory constraints
- **Overflow Protection**: Detection and handling of count overflow
- **Data Integrity**: Validation of data structures and operations
- **Progress Monitoring**: Real-time monitoring of processing progress

## Quality Control Features
- **False Positive Rates**: Configurable trade-offs between memory and accuracy
- **Count Verification**: Validation of count accuracy and consistency
- **Memory Monitoring**: Real-time memory usage tracking
- **Performance Metrics**: Detailed performance statistics and optimization guidance
- **Error Detection**: Automatic detection of potential issues and errors

## Applications
- **Sequence Assembly**: K-mer counting for genome assembly and validation
- **Error Correction**: Preprocessing step for improving sequence quality
- **Contamination Detection**: Identification of foreign sequences in datasets
- **Repeat Analysis**: Analysis of repetitive elements and copy number
- **Quality Control**: Assessment of sequence quality and complexity
- **Metagenomics**: Analysis of complex microbial communities
- **Comparative Genomics**: Comparison of k-mer profiles between samples