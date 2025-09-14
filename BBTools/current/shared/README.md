# BBTools shared Package

Shared utilities and core infrastructure components providing foundational functionality across BBTools suite, including random number generation, text parsing, SIMD operations, synchronization primitives, and system utilities.

---

## Colors (Colors.java)
**Purpose**: Utility class for terminal text color formatting and color array generation
**Core Function**: Provides ANSI color codes and methods to format text with colors and styles
**Key Features**:
- Supports standard and bright color variations (red, green, blue, etc.)
- Handles text formatting like underline and color reset
- Platform-aware color support (disables on unsupported systems)
- Generates predefined color arrays for easy use
**Usage**: Used for colorizing terminal output across BBTools applications

## FastRandom (FastRandom.java)
**Purpose**: A fast, seedable random number generator for non-cryptographic purposes.
**Core Function**: Implements XorShift128+ algorithm for generating high-quality pseudorandom numbers with low computational overhead.
**Key Features**:
- Uses XorShift128+ algorithm for fast random number generation
- Supports seeding with custom or system time-based seeds
- Generates pseudorandom values for long, int, float, double, and boolean types
- Implements efficient byte array population method
**Usage**: Recommended for performance-critical applications requiring fast, non-cryptographic random number generation.

## FastRandomAES (FastRandomAES.java)
**Purpose**: Cryptographically secure pseudorandom number generator using AES encryption
**Core Function**: Generates high-quality random numbers using AES cipher in counter mode
**Key Features**:
- Implements efficient random number generation methods (nextInt, nextLong, etc.)
- Uses AES encryption for cryptographic randomness
- Supports seeding with custom or system-generated seeds
- Optimized buffer-based random number generation
**Usage**: Used for generating secure, high-performance random numbers in cryptographic and computational applications

## FastRandomSIMD (FastRandomSIMD.java)
**Purpose**: SIMD-based pseudorandom number generator extending java.util.Random
**Core Function**: Generates random numbers using vectorized operations with SIMD instructions
**Key Features**:
- Uses 256-bit vector processing for random number generation
- Implements efficient random number generation methods
- Supports multiple random number generation strategies
- Includes optimized methods for nextInt(), nextLong(), etc.
**Usage**: Generating high-performance pseudorandom numbers in performance-critical applications

## FastRandomXoshiro (FastRandomXoshiro.java)
**Purpose**: A fast, seedable random number generator based on xoshiro256+ algorithm
**Core Function**: Generates high-quality pseudorandom numbers with improved statistical properties
**Key Features**:
- Extends java.util.Random for standard random generation
- Uses xoshiro256+ algorithm for efficient random number generation
- Supports seeding with system time or custom seed
- Provides methods for generating random int, long, float, double, and boolean values
**Usage**: Used in performance-critical applications requiring fast, high-quality random number generation

## KillSwitch (KillSwitch.java)
**Purpose**: Monitors CPU utilization to determine if the program has crashed.
**Core Function**: Provides runtime monitoring, forced VM shutdown, and safe memory allocation mechanisms.
**Key Features**:
- CPU utilization tracking with customizable thresholds
- Forced process termination when system load exceeds limits
- Safe memory allocation methods with OutOfMemory error handling
- Static utility methods for array allocation and copying
**Usage**: Ensures program stability by detecting and preventing hung or unresponsive processes

## LineParser (LineParser.java)
**Purpose**: Interface for parsing delimited byte[] lines with flexible term extraction methods.
**Core Function**: Provides methods to parse and manipulate delimited text lines across various data types.
**Key Features**:
- Parse terms as int, long, float, double, byte, and String
- Set and reset parsing boundaries
- Check term characteristics (length, start, equality)
- Convert parsed line to ArrayList of Strings
**Usage**: Used for processing delimited text data with fine-grained parsing control.

## LineParser1 (LineParser1.java)
**Purpose**: Finds delimiters of a text line efficiently to allow for parsing.
**Core Function**: Splits text lines into terms using a specified delimiter, enabling fast parsing of individual terms.
**Key Features**:
- Parses various data types (int, long, float, double, byte)
- Supports parsing specific terms by index
- Memory-efficient delimiter tracking
- Flexible term extraction methods
**Usage**: Used for parsing delimited text lines with high performance, particularly in data processing scenarios

## LineParser2 (LineParser2.java)
**Purpose**: Line parsing utility with bounded memory footprint for processing delimited text lines
**Core Function**: Parses byte arrays into individual terms using a specified delimiter
**Key Features**:
- Supports parsing various primitive types (int, long, float, double)
- Allows advancing through line terms incrementally
- Handles line parsing with minimal memory overhead
- Implements LineParser interface for flexible text processing
**Usage**: Parsing delimited text files or data streams with controlled memory consumption

## LineParser4 (LineParser4.java)
**Purpose**: Uses multiple ordered delimiters for parsing lines of text.
**Core Function**: Parses byte arrays into terms using specified delimiters, allowing flexible text splitting.
**Key Features**:
- Supports multiple ordered delimiters (e.g., ",. ,,")
- Parse methods for various data types (int, long, float, double)
- Tracks term boundaries using IntList
- Supports parsing specific terms or current field
**Usage**: Parsing structured text files with complex delimiter patterns

## LineParserS (LineParserS.java)
**Purpose**: Interface for parsing delimited string lines with flexible parsing capabilities.
**Core Function**: Provides methods to set, parse, and manipulate text lines with controlled parsing.
**Key Features**:
- Set parsing target line with optional term limit
- Parse specific characters from designated fields
- Reset parsing state
- Chainable parsing methods
**Usage**: Used in text processing and parsing scenarios requiring granular line manipulation

## LineParserS1 (LineParserS1.java)
**Purpose**: Finds delimiters of a text line efficiently to enable parsing
**Core Function**: Splits text lines using a specified delimiter, allowing term-based parsing
**Key Features**:
- Parses integers, longs, floats, doubles from specific line terms
- Supports parsing characters and byte arrays from line segments
- Allows checking term contents and term length
- Memory proportional to number of delimiters per line
**Usage**: Used for efficiently parsing structured text lines by breaking them into terms using a delimiter

## LineParserS2 (LineParserS2.java)
**Purpose**: Line parsing utility with bounded memory footprint for parsing delimited text files
**Core Function**: Parses lines into terms using a specified delimiter, supporting various type conversions
**Key Features**:
- Supports parsing integers, longs, floats, doubles, and strings
- Allows advancing through line terms dynamically
- Implements methods for comparing and manipulating line terms
- Bounded memory usage for processing very long lines
**Usage**: Used for efficient parsing of delimited text files with controlled memory consumption

## LineParserS3 (LineParserS3.java)
**Purpose**: Line parsing utility that implicitly uses tab and space as delimiters
**Core Function**: Parse and process text lines with flexible delimiter handling
**Key Features**:
- Supports custom delimiter character
- Parses lines into terms with multiple parsing methods
- Handles integer, long, float, double, and string parsing
- Provides bounds tracking and term manipulation
**Usage**: Used for parsing structured text files with variable delimiters

## LineParserS4 (LineParserS4.java)
**Purpose**: Uses multiple ordered delimiters for parsing lines into terms
**Core Function**: Parses lines using a sequence of delimiters to split text into discrete terms
**Key Features**:
- Supports parsing multiple terms from a single line
- Provides type-specific parsing (int, long, float, double)
- Handles string and byte array extraction
- Flexible delimiter-based line parsing
**Usage**: Used for tokenizing and extracting structured data from text lines with complex delimiter patterns

## LineParserS4Reverse (LineParserS4Reverse.java)
**Purpose**: Parses lines right-to-left using multiple ordered delimiters
**Core Function**: Breaks down lines from end to start while maintaining left-to-right delimiter interpretation
**Key Features**:
- Supports multiple delimiter parsing
- Reverses delimiter order for right-to-left processing
- Parses various data types (int, long, float, string)
- Handles dynamic line parsing with flexible bounds
**Usage**: Used for parsing lines where end structure is known but prefix is uncertain

## LineParserSimple (LineParserSimple.java)
**Purpose**: Simple byte array line parsing utility that segments lines by a specified delimiter
**Core Function**: Iterates through byte arrays, breaking lines into segments at a predefined delimiter
**Key Features**:
- Supports custom delimiter byte definition
- Tracks parsing state with internal index variables
- Provides advance() and advanceInner() methods for line segmentation
- Handles byte array line parsing without string conversion
**Usage**: Used for efficient, low-overhead parsing of byte-based text lines in performance-critical applications

## MetadataWriter (MetadataWriter.java)
**Purpose**: Writes execution metadata for BBTools workflows in TSV or JSON format
**Core Function**: Captures runtime information including host, version, command, and read statistics
**Key Features**:
- Generates metadata with time, host, BBTools version, and Java version
- Supports both TSV and JSON output modes
- Captures input/output read and base counts
- Translates command-line tools to shell script equivalents
**Usage**: Logging metadata for BBTools computational workflows and tracking job execution details

## Parse (Parse.java)
**Purpose**: Utility class for parsing various data types from strings and byte arrays
**Core Function**: Provides flexible parsing methods for numbers, booleans, arrays, and special string formats
**Key Features**:
- Supports parsing integers, longs, floats, doubles with multiple formats
- Handles numeric suffixes like K, M, G for scale (e.g., 5K = 5000)
- Parses boolean values with flexible input rules
- Provides methods for parsing arrays and special string formats
**Usage**: Used for converting string representations to numeric types in BBTools parsing operations

## Parser (Parser.java)
**Purpose**: Utility class for parsing command-line arguments and configuring BBTools processing parameters
**Core Function**: Provides extensive parsing methods for various input flags and options across multiple BBTools modules
**Key Features**:
- Supports parsing of quality trimming parameters
- Handles input/output file configurations
- Manages read filtering and processing options
- Configures interleaved read processing
- Supports SAM/BAM file parameter settings
**Usage**: Used as a central configuration parser for most BBTools command-line tools, enabling flexible read processing and analysis workflows

## PreParser (PreParser.java)
**Purpose**: Pre-processes command line arguments for BBTools applications
**Core Function**: Handles command line argument parsing, version printing, output stream redirection, and help flag processing
**Key Features**:
- Strips leading hyphens from command line arguments
- Supports Java flag parsing
- Handles output stream redirection
- Generates optional JSON metadata for command execution
**Usage**: Used as a preprocessing step for parsing command line arguments in BBTools applications

## Primes (Primes.java)
**Purpose**: Utility class for finding prime numbers around a given value
**Core Function**: Provides methods to find primes at least or at most a specified number
**Key Features**:
- Finds nearest prime using binary search on pre-computed prime list
- Supports finding primes for various integer ranges
- Handles large prime number searches recursively
- Loads pre-computed prime list from compressed text file
**Usage**: Used for mathematical computations requiring prime number lookup or generation

## Shared (Shared.java)
**Purpose**: Centralized utility functions and constants used across BBTools.
**Core Function**: Provides global configuration, environment detection, resource management, and utility methods for BBTools applications.
**Key Features**:
- Manages thread, buffer, and memory configurations dynamically
- Detects runtime environment (OS, hardware, cluster systems)
- Provides flexible sorting algorithms with parallel/sequential options
- Handles JNI library loading and server URL management
**Usage**: Used as a central configuration and utility class for managing runtime parameters across the entire BBTools suite.

## SIMD (SIMD.java)
**Purpose**: Holds SIMD (Single Instruction, Multiple Data) vectorization methods for performance optimization.
**Core Function**: Provides high-performance vector operations using Java's Vector API for various data types.
**Key Features**:
- Vectorized mathematical operations (sum, max, add, multiply)
- Performance-optimized methods for float, int, long, byte, and double arrays
- SIMD loop implementations with residual scalar processing
- Supports vector species for 256-bit vector operations
**Usage**: Used to accelerate computational kernels in BBTools, particularly in machine learning and data processing components.

## SIMDAlign (SIMDAlign.java)
**Purpose**: Holds SIMD (Single Instruction, Multiple Data) methods for high-performance alignment calculations
**Core Function**: Implements vectorized alignment algorithms using Java Vector API for parallel score computation
**Key Features**:
- Supports vectorized band alignment with SIMD operations
- Handles match, substitution, insertion, and deletion scoring
- Computes maximum alignment scores using vector-based computations
- Provides specialized methods for different data types (long, int)
**Usage**: Used in high-throughput bioinformatics alignment algorithms requiring efficient score computations

## SIMDAlignByte (SIMDAlignByte.java)
**Purpose**: Holds SIMD methods for sequence alignment using vector operations
**Core Function**: Implements high-performance byte-level sequence alignment using Java's Vector API
**Key Features**:
- Supports vectorized substitution counting across multiple diagonals
- Implements band-based sequence alignment using SIMD instructions
- Handles vector and scalar processing modes dynamically
- Supports alignment with clipping and substitution limits
**Usage**: Used in bioinformatics sequence alignment algorithms requiring fast, parallel processing of genomic data

## SIMDAlignOld (SIMDAlignOld.java)
**Purpose**: Holds SIMD methods for efficient sequence alignment using vector operations.
**Core Function**: Implements vectorized alignment algorithms using Java's Vector API for long and integer data types.
**Key Features**:
- Supports vector-based band alignment with match, substitution, and insertion scoring
- Handles special cases for N-character matches and small vector ranges
- Provides methods for processing deletion tails and cross-cut diagonals
- Implements both vector and scalar fallback processing strategies
**Usage**: Used in genomic sequence alignment algorithms requiring high-performance vector computations.

## SyncHeart (SyncHeart.java)
**Purpose**: Thread-safe synchronization utilities for mutable configuration values in BBTools.
**Core Function**: Provides safe concurrent access to shared state using ReadWriteLock with optimized read performance.
**Key Features**:
- Manages thread-safe getters and setters for multiple configuration parameters
- Uses ReentrantReadWriteLock to allow multiple concurrent readers
- Supports configuration for threading, MPI, file I/O, and buffer settings
- Provides explicit read and write locking mechanisms
**Usage**: Used to safely modify and access shared configuration values across multi-threaded BBTools applications.

## Timer (Timer.java)
**Purpose**: Lightweight timing utility for measuring code execution duration
**Core Function**: Tracks time elapsed between start and stop calls using nanosecond precision
**Key Features**:
- Captures start and stop timestamps using System.nanoTime()
- Converts elapsed time to seconds with configurable decimal precision
- Optional output stream for logging timing results
- Supports chained timing with stopAndStart() method
**Usage**: Performance measurement and code execution time tracking in development and testing

## Tools (Tools.java)
**Purpose**: Utility class providing static methods for file handling, data processing, and string formatting
**Core Function**: Implements generic helper methods for input/output operations, type conversions, and system-independent utilities
**Key Features**:
- File input/output validation methods (testInputFiles, testOutputFiles)
- String and number formatting utilities (format, padLeft, padRight)
- Linear regression calculation for data analysis
- Integer list and set loading from various sources
**Usage**: Core utility class used across BBTools for cross-cutting concerns like file processing and data manipulation

## TrimRead (TrimRead.java)
**Purpose**: Helper class for inline quality trimming of genetic sequencing reads
**Core Function**: Dynamically removes low-quality bases from the start and end of DNA/RNA reads
**Key Features**:
- Quality-aware base trimming with multiple modes (optimal, window, standard)
- Supports both left and right-side read trimming
- Handles read match string modification during trimming
- Configurable minimum quality and length thresholds
**Usage**: Used in bioinformatics pipelines to preprocess sequencing data by removing unreliable base calls

## Vector (Vector.java)
**Purpose**: Protects normal classes from seeing SIMD implementation details
**Core Function**: Provides vector-based computational methods with optional SIMD acceleration
**Key Features**:
- Implements array manipulation methods like sum, max, add
- Supports both scalar and vector (SIMD) computation paths
- Handles multiple numeric array types (int, float, long, etc.)
- Provides performance-optimized array operations
**Usage**: Used for high-performance mathematical and machine learning computations