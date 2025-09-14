# Cardinality Package

The cardinality package provides probabilistic data structures for estimating the cardinality (number of unique elements) in large datasets, primarily using LogLog-family algorithms for memory-efficient cardinality tracking.

## Classes

#### BBLog (BBLog.java)
**Purpose**: Implements a probabilistic cardinality estimation algorithm using the LogLog method for counting unique elements in large datasets.

- **Core Function**: Provides an efficient probabilistic technique for estimating the number of unique elements using a fixed-size array of buckets and hash-based tracking
- **Key Features**:
  - Supports both atomic and non-atomic tracking modes for thread-safe and non-thread-safe scenarios
  - Configurable number of buckets (default 2048) and k-mer length (default 31)
  - Uses 64-bit hash function (Tools.hash64shift) for distributing elements across buckets
  - Implements multiple cardinality estimation strategies (total, estimated sum, median estimation)
  - Optional count tracking for each bucket to provide additional statistical information
  - Supports merging cardinality trackers with maximum value selection
- **Usage**: Primarily used in bioinformatics and data processing applications for quickly estimating the number of unique elements in large datasets without storing all elements, such as unique k-mer counting in genomic sequences

#### BBLog_simple (BBLog_simple.java)
**Purpose**: A simplified LogLog-based cardinality estimation tracker for efficiently estimating the number of unique elements in a set.

- **Core Function**: Implements a probabilistic data structure that uses a hashing technique to estimate the number of unique elements with low memory overhead
- **Key Features**:
  - Configurable number of buckets (default 2048) for tracking unique elements
  - Supports kmer length specification (default k=31)
  - Allows random seed initialization for reproducibility
  - Efficient hash-based cardinality estimation algorithm
  - Optional count tracking for additional analysis
  - Supports merging of multiple cardinality trackers
- **Usage**: Used in bioinformatics and data processing scenarios where approximate unique element counting is needed, such as:
  - Estimating unique k-mers in genomic sequences
  - Determining approximate set sizes with minimal memory consumption
  - Providing fast, probabilistic unique element counting for large datasets

#### CardinalityTracker (CardinalityTracker.java)
**Purpose**: Abstract superclass for cardinality-tracking data structures that estimate the number of unique elements in a dataset.

- **Core Function**: Implements probabilistic counting algorithms to estimate the number of distinct elements using minimal memory, supporting various LogLog-based tracking methods.

- **Key Features**:
  - Factory method support for creating different types of cardinality trackers (BBLog, LogLog, LogLog2, LogLog16, LogLog8)
  - Flexible hashing mechanisms for tracking both small and large k-mers 
  - Quality-aware probabilistic tracking with minimum probability thresholds
  - Supports both DNA sequence and numeric element tracking
  - Thread-safe design with thread-local k-mer reuse
  - Configurable bucket count for trade-off between accuracy and memory usage
  - Ability to generate k-mer frequency histograms

- **Usage**: Used in bioinformatics applications to estimate unique k-mer counts in DNA sequences, genomic data analysis, and other scenarios requiring memory-efficient distinct element counting.

#### LogLog (LogLog.java)
**Purpose**: Implements the LogLog cardinality estimation algorithm for efficient unique element counting using probabilistic data structures.

- **Core Function**: Estimates the number of unique elements (cardinality) using a space-efficient probabilistic counting technique based on tracking the maximum number of leading zeros in hash values across multiple buckets.

- **Key Features**:
  - Supports both thread-safe (atomic) and non-atomic implementations of cardinality tracking
  - Utilizes a configurable number of buckets (default 2048) for tracking unique element estimates
  - Implements multiple cardinality estimation methods: standard cardinality(), cardinalityH(), with different statistical approaches
  - Uses 64-bit hash function (Tools.hash64shift) for distributing elements across buckets
  - Supports adding (merging) multiple LogLog trackers with maximum value selection
  - Handles quality filtering of elements based on minimum probability threshold
  - Provides compensation factors for log-scaled bucket corrections

- **Usage**: Typically used in bioinformatics and big data applications to efficiently estimate the number of unique elements (such as unique k-mers) without storing all elements, making it ideal for memory-constrained environments that require approximate unique element counting.

#### LogLog16 (LogLog16.java)
**Purpose**: Implements a 16-bit LogLog cardinality estimation algorithm for efficiently counting unique elements with low memory overhead.

- **Core Function**: Uses a probabilistic counting method to estimate the number of unique elements in a dataset by tracking the maximum number of leading zeros in hash values across multiple buckets
- **Key Features**:
  - Supports configurable number of buckets (default 2048) for cardinality tracking
  - Uses 64-bit hash function (Tools.hash64shift) for generating hash keys
  - Implements adaptive cardinality estimation with multiple statistical methods (mean, median, harmonic mean)
  - Allows dynamic configuration of kmer length and minimum probability filtering
  - Supports adding and merging LogLog16 trackers for distributed counting
  - Provides constant-time O(1) space complexity with low memory requirements
- **Usage**: Primarily used in bioinformatics for estimating unique k-mer counts in genomic data, providing a memory-efficient alternative to exact counting methods for large datasets

#### LogLog2 (LogLog2.java)
**Purpose**: Second-generation LogLog cardinality estimation algorithm for efficient unique element counting in large datasets.

**Core Function**: Implements a probabilistic data structure to estimate the number of unique elements using a logarithmic counting technique with advanced hashing and statistical estimation strategies.

**Key Features**:
- **Configurable Bucket System**: Supports custom number of buckets (default 2048) for flexible cardinality tracking
- **Advanced Hashing Mechanism**: Uses 64-bit hash function (Tools.hash64shift) to distribute elements across buckets
- **Multiple Estimation Strategies**: Supports different cardinality estimation methods including mean, median, harmonic mean, and geometric mean
- **Atomic and Non-Atomic Modes**: Provides thread-safe implementation with AtomicIntegerArray for concurrent environments
- **Precision Control**: Configurable mantissa bits (default 20) for fine-tuning estimation accuracy
- **Quality Filtering**: Option to ignore elements below a specified probability threshold

**Usage**: Efficiently estimate the number of unique elements in large datasets, such as genomic k-mer counting, stream processing, and cardinality estimation in big data analytics without storing all unique elements.

#### LogLog8 (LogLog8.java)
**Purpose**: An 8-bit implementation of the LogLog cardinality estimation algorithm for efficient, approximate counting of unique elements.

- **Core Function**: Provides a memory-efficient probabilistic data structure for estimating the number of unique elements in a stream using logarithmic counting and minimal memory overhead.

- **Key Features**:
  - Supports configurable number of buckets (default 2048) for cardinality estimation
  - Uses 8-bit byte array for tracking maximum leading zero counts of hashed elements
  - Implements custom hash function (Tools.hash64shift) for element distribution
  - Supports adding multiple LogLog8 trackers together by taking maximum of corresponding buckets
  - Provides probabilistic cardinality estimation with configurable Kmer length and random seed
  - Low memory footprint with approximately O(log log n) space complexity
  - Configurable minimum probability threshold for filtering low-confidence elements

- **Usage**: 
  - Ideal for estimating unique element count in large datasets where exact counting is memory-prohibitive
  - Commonly used in bioinformatics for k-mer counting and unique sequence estimation
  - Supports quick approximate cardinality estimation with minimal computational overhead
  - Can be used in streaming data scenarios where memory is constrained and approximate counts are acceptable

#### LogLog8_simple (LogLog8_simple.java)
**Purpose**: A lightweight, simplified 8-bit implementation of the LogLog cardinality estimation algorithm for efficiently tracking unique element counts.

- **Core Function**: Uses an 8-bit probabilistic data structure to estimate the number of unique elements (cardinality) in a dataset with minimal memory overhead
- **Key Features**:
  - Supports configurable number of buckets (default 2048)
  - Tracks unique elements using a compact byte array representation
  - Implements hash-based counting with leading zero tracking
  - Provides efficient cardinality estimation with low memory requirements
  - Supports merging multiple LogLog trackers
  - Uses 64-bit hash function for uniform element distribution
  - Allows optional minimum probability filtering for element inclusion
- **Usage**: Ideal for approximate counting of unique elements in large datasets, such as unique k-mer tracking in genomic sequences, web analytics, network flow analysis, and other scenarios requiring memory-efficient unique element estimation

#### LogLogWrapper (LogLogWrapper.java)
**Purpose**: Implements a flexible cardinality estimation wrapper using the LogLog algorithm for unique k-mer counting in biological sequence data.

- **Core Function**: Estimates the number of unique k-mers in DNA/RNA sequence datasets using a probabilistic counting technique with configurable parameters
- **Key Features**:
  - Supports multiple input file processing (single or paired-end reads)
  - Configurable number of buckets (default 2048) for cardinality estimation
  - Supports both real sequence file and synthetic read generation modes
  - Multi-threaded k-mer hashing and cardinality tracking
  - Configurable k-mer length (default 31)
  - Optional count tracking and detailed statistical reporting
  - Handles various input formats (FASTQ supported)

- **Usage**: Used for rapidly estimating the number of unique k-mers in large genomic or transcriptomic sequencing datasets, particularly useful for assessing sequence diversity and complexity without full set enumeration

#### LogLog_old (LogLog_old.java)
**Purpose**: A legacy implementation of the LogLog cardinality estimation algorithm for probabilistic counting of unique elements in data streams.

- **Core Function**: Estimate the number of unique elements (cardinality) in a dataset using a probabilistic hashing technique that uses logarithmic space complexity

- **Key Features**:
  - Uses multiple hash tables with configurable parameters (buckets, bits, k-mer length)
  - Supports quality-aware probabilistic hashing for biological sequence data
  - Handles both small (k<32) and large (k>=32) k-mer based probabilistic counting
  - Supports concurrent processing of reads using thread-local storage
  - Implements both atomic and non-atomic counting modes
  - Allows filtering of low-probability k-mers based on base quality scores

- **Usage**: 
  - Primarily used in bioinformatics for estimating unique k-mer counts in sequencing data
  - Can estimate cardinality of large datasets with minimal memory overhead
  - Useful for genomic sequence analysis, metagenomic studies, and sequence deduplication tasks

#### MultiLogLog (MultiLogLog.java)
**Purpose**: Manages multiple LogLog cardinality trackers for simultaneous k-mer size tracking in genomic sequencing data.

**Core Function**: Creates and manages an array of CardinalityTracker instances to estimate unique k-mer counts across different k-mer lengths simultaneously.

**Key Features**:
- Supports multiple k-mer lengths in a single instance
- Configurable via Parser or direct parameter specification
- Validates and filters k-mer lengths dynamically
- Allows custom number of buckets for cardinality estimation
- Supports optional random seed for hash function
- Provides minimum probability threshold for k-mer quality filtering
- Automatically sorts and deduplicates k-mer lengths

**Usage**: Used in genomic sequencing analysis to estimate the number of unique k-mers across multiple k-mer sizes, helping assess sequence diversity and complexity. Typical workflow involves creating a MultiLogLog instance with desired k-mer lengths and then hashing reads to track cardinality.