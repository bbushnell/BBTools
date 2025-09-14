# fun Package - Utility Tools and Experimental Functions
*Collection of utility classes, experimental algorithms, and testing tools for various computational tasks*

## Package Overview
The `fun` package contains diverse utility classes for mathematical calculations, disk performance analysis, pathfinding algorithms, genetic simulations, palindrome detection, and experimental data processing tools.

---

## Mathematical and Statistical Utilities

### Calc (Calc.java)
**Purpose**: Command-line statistical calculation tool that generates cumulative distribution statistics for bitwise-encoded integer combinations.

- **Key Features**:
  - Efficient bitwise statistical aggregation
  - Flexible command-line configuration
  - Configurable verbosity and output streaming
  - Performance timing and comprehensive error tracking

- **Usage**:
  ```bash
  java fun.Calc numstats=6 verbose=t
  ```

- **Technical Details**:
  - Generates distribution by iterating 2^(numStats*5) possible bit combinations
  - Supports custom statistical bucket configurations
  - Outputs cumulative percentage distribution to standard output

### Chance (Chance.java)
**Purpose**: Statistical probability calculator for calculating the likelihood of achieving a minimum number of successes in a given number of draws.

- **Core Function**: Performs Monte Carlo simulations to estimate the probability of achieving a specified number of successes across multiple experimental rounds

- **Key Features**:
  - Calculates probability using random sampling
  - Supports configurable parameters for draws, minimum successes, and probability
  - Uses thread-local random number generation for improved performance
  - Provides precise probability estimation through multiple simulation rounds

- **Technical Implementation**:
  - Utilizes Monte Carlo simulation method
  - Aggregates results across multiple simulation rounds
  - Outputs probability as a percentage with 6 decimal places

- **Example Invocation**:
  ```bash
  # Calculate probability of at least 3 successes in 10 draws with 0.5 probability
  java fun.Chance 10 3 0.5 100000
  ```

### ParseDouble (ParseDouble.java)
**Purpose**: Provides custom, high-performance implementations for parsing byte arrays into double-precision floating-point numbers with advanced handling of numeric conversions.

- **Key Features**:
  - **Efficient Numeric Parsing**: Two custom parsing methods (`parseDouble` and `parseDouble2`) that directly convert byte arrays to double values
  - **Performance Optimization**: Includes built-in warmup and benchmarking mechanism
  - **Special Value Handling**: Supports parsing of NaN, Infinity, and negative numbers
  - **Flexible Input Processing**: Handles decimal points and complex number formats

- **Technical Implementation Details**:
  - Uses pre-computed `DECIMAL_INV_MULT` lookup table for efficient decimal place calculations
  - Implements custom numeric map (`NUMERIC_MAP`) for rapid character validation
  - Falls back to `Double.parseDouble()` for complex number formats

---

## Probability and Simulation Tools

### ProbShared (ProbShared.java)
**Purpose**: Calculates probabilistic k-mer intersection and cardinality for genomic sequence analysis.

- **Core Function**: Computes statistical probabilities of k-mer overlaps between two DNA/RNA sequences of different lengths
  - Uses combinatorial probability to estimate unique k-mer counts
  - Calculates intersection probability between sequence k-mer sets

- **Key Features**:
  - `cardinality(int k, int seqLength)`: Estimates unique k-mer count in a sequence
  - `probIntersect(int k, int len1, int len2)`: Determines probability of k-mer set intersection

- **Technical Implementation**:
  - Uses `Math.pow(4, k)` to calculate k-mer sequence space (4 DNA bases)
  - Employs iterative probability reduction algorithm

### ProbShared2 (ProbShared2.java)
**Purpose**: Probabilistic simulation of k-mer shared sequence occurrence between two random DNA sequences.

- **Core Function**: Calculates the probability of finding a common k-mer between two randomly generated sequences of specified lengths.

- **Key Features**:
  - Monte Carlo simulation-based probability estimation
  - Random DNA sequence generation
  - K-mer frequency analysis
  - Configurable sequence lengths and simulation rounds

- **Technical Implementation**:
  - Uses bit manipulation for efficient k-mer encoding
  - Employs HashSet for unique k-mer tracking
  - Supports variable k-mer sizes and sequence lengths

### ProbShared3 (ProbShared3.java)
**Purpose**: Probabilistic simulation tool for k-mer occurrence probability across random sequences.

- **Core Function**: Estimates the probability of finding shared k-mers between two randomly generated sequences

- **Key Features**:
  - Monte Carlo simulation of k-mer overlap
  - Configurable parameters for sequence length and k-mer size
  - Statistical sampling using random generation
  - Precise probability calculation through multiple simulation rounds

---

## System Performance and Analysis

### DiskBench (DiskBench.java)
**Purpose**: A sophisticated disk I/O performance benchmarking tool for measuring read, write, and read-write performance across multiple file I/O methods.

- **Key Features**:
  - **Flexible Benchmarking Modes**:
    - Read-only mode
    - Write-only mode
    - Simultaneous read-write mode
  - **Multiple I/O Method Support**:
    - ByteFile, QuickFile, TextFile
    - BufferedInputStream, FileInputStream
  - **Configurable Benchmark Parameters**:
    - Adjustable data size
    - Configurable number of passes
    - Multi-threaded performance testing
  - **Performance Reporting**:
    - Detailed timing information
    - Throughput calculations in MB/s

- **Sample Invocation**:
  ```bash
  java fun.DiskBench mode=rw data=1G passes=3 method=bytefile
  ```

### DiskSpaceAnalyzer (DiskSpaceAnalyzer.java)
**Purpose**: Recursive disk space analysis tool that generates an HTML report of directory sizes and storage utilization.

- **Core Function**: 
  - Recursively scans a root directory to calculate total file and subdirectory sizes
  - Generates a detailed HTML report showing disk usage across different directories

- **Key Features**:
  - **Recursive Directory Traversal**: Analyzes entire directory tree
  - **Size Calculation**: Computes total bytes for each directory and individual files
  - **HTML Reporting**: Creates visually formatted disk usage report
  - **Flexible Root Path**: Allows specifying custom root directory

- **Example Usage**:
  ```bash
  java fun.DiskSpaceAnalyzer /path/to/analyze
  # Generates disk_usage.html in current directory
  ```

---

## Algorithms and Data Structures

### FindPath (FindPath.java)
**Purpose**: A lightweight graph pathfinding utility that computes the shortest path between nodes in an undirected graph.

- **Core Function**: 
  - Compute the shortest path between two nodes using a breadth-first graph traversal strategy
  - Loads graph connectivity from a tab-separated text file
  - Supports finding paths between arbitrary nodes with minimal computational overhead

- **Key Features**:
  - Dynamic graph construction from input file
  - Shortest path computation with distance tracking
  - Bidirectional edge support (graph is undirected)
  - Provides detailed path output with node sequence and total distance

- **Command-Line Usage**:
  ```bash
  java fun.FindPath StartNode EndNode graph.txt
  ```

### MakeAdjacencyList (MakeAdjacencyList.java)
**Purpose**: Generates a randomized adjacency matrix representing a graph with configurable parameters, writing the result to a text file.

- **Core Function**: Create a probabilistic graph with randomly generated edge weights between nodes

- **Key Features**:
  - Configurable number of nodes via `nodes` parameter (default 10)
  - Random edge generation with probability control (`prob` parameter, default 0.3)
  - Flexible edge weight range using `minlen` and `maxlen` parameters (default 5-25)
  - Deterministic graph generation using optional random seed

- **Example Command-Line Usage**:
  ```bash
  java fun.MakeAdjacencyList nodes=15 prob=0.4 minlen=10 maxlen=50 out=graph.txt
  ```

---

## Genetic Algorithms and Simulations

### Genetic (Genetic.java)
**Purpose**: Genetic algorithm implementation for solving optimization problems using evolutionary computation techniques.

- **Core Function**: Uses genetic algorithm principles to evolve solutions over multiple generations

- **Key Features**:
  - Population-based evolutionary optimization
  - Configurable mutation and crossover rates
  - Fitness-based selection mechanisms
  - Multi-generational solution evolution

### Life (Life.java)
**Purpose**: Implements Conway's Game of Life cellular automata simulation on a configurable 2D grid with random initial state.

- **Core Function**: 
  - Creates a 2D grid of cells that evolve through discrete generations based on specific neighborhood rules
  - Simulates cell life/death cycles using neighbor-based state transitions

- **Key Features**:
  - **Random Initialization**: Populates initial grid with live/dead cells based on provided probability
  - **Toroidal Grid**: Uses modular arithmetic to create wrap-around grid boundaries
  - **State Tracking**: Detects repeating grid states to prevent infinite loops
  - **Terminal Visualization**: Prints grid state using '@' and space characters

- **Example Execution**:
  ```bash
  java fun.Life 40 20 100 0.25
  # Creates 40x20 grid, runs 100 rounds, 25% initial live probability
  ```

---

## String Processing and Analysis

### Palindrome (Palindrome.java)
**Purpose**: Advanced string processing utility for finding the longest palindrome in sequences with configurable mismatch tolerance and optional reverse complement support.

- **Core Function**: Locates the longest palindromic substring within input strings, with advanced matching capabilities including mismatch allowance and reverse complement detection.

- **Key Features**:
  - Finds longest palindrome with support for both even and odd-length palindromes
  - Configurable mismatch tolerance (maxMismatches parameter)
  - Optional reverse complement (rcomp) mode for biological sequence analysis
  - Supports multiple input methods: direct string, file input, command-line arguments

- **Command-Line Interface**:
  ```bash
  # Basic usage
  java Palindrome "racecar"
  
  # With mismatch tolerance
  java Palindrome 2 "abcdefg"
  
  # Reverse complement mode
  java Palindrome rcomp "ATCG"
  ```

### Palindrome2 (Palindrome2.java)
**Purpose**: Advanced palindrome detection utility for identifying the longest palindromic substring with configurable mismatch tolerance and sequence processing.

- **Core Function**: Finds the longest palindromic substring with advanced features like mismatch allowance and complementary sequence detection.

- **Key Features**:
  - **Flexible Palindrome Detection**: Supports both odd and even-length palindrome detection
  - **Multiple Input Methods**: Direct string input, file-based sequence reading
  - **Advanced Palindrome Length Calculation**: Separate methods for odd and even-centered sequences
  - **Sequence Preprocessing**: Handles multi-sequence files with FASTA-like headers

### PalSim (PalSim.java)
**Purpose**: Monte Carlo simulation for analyzing product distribution and placement of random float triples.

- **Core Function**: Performs a probabilistic simulation to calculate placement and distribution of three-dimensional product spaces

- **Key Features**:
  - Generates random float triples within a constrained range (1.0 to 1.3)
  - Runs multiple simulations (default 1,000,000 iterations)
  - Calculates multiple placement metrics

- **Execution Example**:
  ```bash
  java fun.PalSim 0.1 0.2 0.3 500000
  ```

---

## Data Processing Utilities

### Crunch (Crunch.java)
**Purpose**: File size aggregation utility for processing file metadata from a pipe-delimited input file.

- **Core Function**: Calculates total size of files matching specific criteria by parsing a metadata file

- **Key Features**:
  - Reads file metadata from pipe-delimited text file
  - Filters files containing 'F' in the 7th column
  - Aggregates total file size across matching files
  - Converts final size to terabytes (TB)

### Dongle (Dongle.java)
**Purpose**: Utility class for time-based encoding, decoding, and validation of long integer values using bitwise manipulation and rotation techniques.

- **Core Function**: Provides specialized encoding and decoding mechanisms for timestamps with built-in validation and parsing capabilities

- **Key Features**:
  - **Time Encoding**: Transforms timestamps using complex bitwise XOR and bit rotation operations
  - **Time Decoding**: Reverses the encoding process to retrieve original timestamp values
  - **Date Parsing**: Converts string representations of dates to millisecond timestamps
  - **Timestamp Validation**: Checks if a timestamp falls within a specified time range

### Merced (Merced.java)
**Purpose**: Demonstration and exploration of Java HashMap usage and basic object-oriented programming concepts in a fun, educational context.

- **Core Function**: Provides a lightweight example of key-value data structures and nested class design

- **Key Features**:
  - Nested `Person` class with demographic attributes
  - HashMap demonstrations with type-specific mappings
  - Simple computational method for calculating personal "attractiveness"

---

## Experimental File Processing (Foo Series)

### Foo (Foo.java)
**Purpose**: A utility class for reading, processing, and analyzing text files with flexible parsing and filtering capabilities.

- **Core Function**: 
  - Reads a text file line by line
  - Parses and filters lines based on specific criteria
  - Generates statistical analysis of file contents, particularly file sizes

### Foo2 (Foo2.java)
**Purpose**: Efficient file processing utility for parsing delimited text files with size-based filtering.

- **Core Function**: Processes pipe-delimited text files, extracting and aggregating file size information with two processing strategies (slow and fast parsing methods)

### Foo3 (Foo3.java)
**Purpose**: File size and access time analysis utility with performance-optimized processing of large file datasets.

- **Core Function**: 
  - Processes large files containing file metadata, extracting file sizes and access times
  - Supports both slow (regex-based) and fast (manual parsing) processing modes
  - Generates statistical summaries of file sizes across different percentile ranges

### Foo4 (Foo4.java)
**Purpose**: File size and access time analysis utility that processes and reports storage access patterns across files.

- **Core Function**: 
  - Reads a delimited text file containing file metadata
  - Analyzes file sizes and last access times
  - Generates statistical summary of file access patterns

### Foo5 (Foo5.java)
**Purpose**: Advanced file metadata processing utility with comprehensive statistical analysis and time-based filtering.

- **Key Features**:
  - Multi-threaded file processing capabilities
  - Advanced statistical calculations including percentiles and distributions
  - Time-based file access pattern analysis

### Foo6 (Foo6.java)
**Purpose**: Reads a text file, filters lines, and generates statistical insights about file access and sizes.

- **Core Function**: 
  - Processes input text files with delimited data
  - Extracts and analyzes file size and access time information
  - Provides detailed statistical breakdowns of file metadata

### Foo7 (Foo7.java)
**Purpose**: A versatile text file processing utility that reads, filters, and analyzes file information with advanced statistical capabilities.

- **Core Function**: 
  - Reads a delimited text file containing file metadata
  - Filters and processes lines based on specific criteria
  - Generates statistical insights about file sizes and access times

---

## Package Usage
The fun package serves as a collection of:
- Mathematical and statistical calculation tools
- System performance analysis utilities
- Experimental algorithms and data processing tools
- Educational programming examples
- Specialized file analysis and metadata processing utilities
- Probabilistic simulation and modeling tools

## Performance Considerations
- Many utilities include performance benchmarking capabilities
- Monte Carlo simulations support configurable precision through iteration counts
- File processing tools offer both fast and thorough parsing strategies
- Memory-efficient implementations for large dataset processing

## Dependencies
- Relies on BBTools shared utilities and data structures
- Uses standard Java libraries for mathematical operations
- Integrates with file I/O utilities from the fileIO package
- Some utilities require specific input file formats

---
*Documentation generated using evidence-based analysis of source code*