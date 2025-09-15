# BBTools template Package

Template classes and utilities providing foundational frameworks for developing BBTools applications, including single-threaded and multi-threaded processing patterns, thread synchronization utilities, and standardized workflows for read processing and data manipulation.

---

## A_Sample (A_Sample.java)
**Purpose**: Template class designed to be easily modified into a program that processes reads in a single thread.
**Core Function**: Provides a structured framework for reading, processing, and writing read data with comprehensive input/output handling.
**Key Features**:
- Supports paired-end and single-end read processing
- Handles file extension and interleaving detection
- Configurable input/output stream management
- Built-in error checking and file validation
**Usage**: Used as a starting template for creating new single-threaded read processing tools in BBTools

## A_Sample2 (A_Sample2.java)
**Purpose**: Sample template class demonstrating basic BBTools read processing workflow
**Core Function**: Processes input reads using concurrent input/output streams with configurable parameters
**Key Features**:
- Supports single and paired-end read processing
- Configurable via command-line arguments
- Tracks read and base processing metrics
- Handles file input/output with format detection
**Usage**: Base template for developing new BBTools processing applications

## A_SampleBasic (A_SampleBasic.java)
**Purpose**: Reads a text file and prints it to another text file.
**Core Function**: Processes text files line by line, extracting content before the first tab and writing to output.
**Key Features**:
- Command-line argument parsing
- File input and output handling
- Line-by-line text processing
- Error tracking and reporting
**Usage**: Template for basic file reading and processing tasks in BBTools

## A_SampleByteFile (A_SampleByteFile.java)
**Purpose**: Reads a text file, prints it to another text file, and filters out invalid lines.
**Core Function**: Processes text files with line-by-line validation and optional invalid line routing
**Key Features**:
- Reads input files with multi-threaded ByteFile processing
- Filters lines based on simple validation rules
- Supports optional output for invalid lines
- Configurable line processing limits
**Usage**: Used as a template/sample for text file processing utilities in BBTools

## A_SampleByteFileMT (A_SampleByteFileMT.java)
**Purpose**: Multithreaded text file processing template with file input/output handling
**Core Function**: Loads text files, spawns worker threads to process input files, and manages concurrent file reading
**Key Features**:
- Supports multiple input/output file formats
- Configurable thread management for file processing
- Line-by-line file reading with validation
- Error tracking and performance timing
**Usage**: Base class for implementing custom multithreaded file processing workflows in BBTools

## A_SampleD (A_SampleD.java)
**Purpose**: Sample template class demonstrating BBTools application processing framework
**Core Function**: Processes sequence reads using distributed stream processing with configurable parameters
**Key Features**:
- Supports parallel read processing with concurrent input/output streams
- Handles both single-end and paired-end sequencing data
- Implements performance timing and MPI-compatible architecture
- Flexible command-line argument parsing with verbose mode
**Usage**: Base template for developing read processing applications in BBTools framework

## A_SampleMT (A_SampleMT.java)
**Purpose**: Template class for multi-threaded read processing programs
**Core Function**: Provides scaffolding for creating multi-threaded read processing applications
**Key Features**:
- Supports concurrent read processing with multiple threads
- Handles input and output stream management
- Includes placeholder for custom read pair processing method
- Tracks reads and bases processed
**Usage**: Base template for developing multi-threaded bioinformatics processing tools

## A_SampleMultipleInput (A_SampleMultipleInput.java)
**Purpose**: Accepts multiple input files and processes them sequentially, outputting to a single file.
**Core Function**: Reads input files, processes reads, and writes them to a consolidated output stream while handling file format detection and concurrent I/O.
**Key Features**:
- Supports multiple input file processing
- Handles FASTQ file format detection
- Uses concurrent read and write streams
- Tracks processing time and read statistics
**Usage**: Used for concatenating files while recompressing and avoiding standard I/O limitations.

## A_SampleSamStreamer (A_SampleSamStreamer.java)
**Purpose**: Provides a template for processing reads in multiple threads by filling in the processRead method
**Core Function**: Multithread SAM/BAM file processing with configurable read filtering and stream management
**Key Features**:
- Supports parallel file reading with configurable threads
- Implements flexible read processing and filtering
- Handles input/output streams for SAM/BAM files
- Provides extensible thread-safe read processing framework
**Usage**: Used as a starting template for developing new read processing programs in BBTools

## A_SampleSummary (A_SampleSummary.java)
**Purpose**: A template class demonstrating basic read processing workflow in BBTools
**Core Function**: Processes input reads, tracks processing metrics, and generates an output summary
**Key Features**:
- Supports concurrent read input stream processing
- Tracks reads and bases processed
- Provides performance timing metrics
- Flexible input/output file handling
**Usage**: Used as a template for creating new BBTools processing scripts with standard read processing pattern

## A_Sample_Generator (A_Sample_Generator.java)
**Purpose**: Template class for generating reads, potentially based on input sequences
**Core Function**: Provides a framework for multi-threaded read generation with configurable input/output processing
**Key Features**:
- Supports parallel input/output processing using concurrent streams
- Configurable read generation with thread-safe atomic read ID tracking
- Flexible file handling with extension and overwrite management
- Implements a multi-threaded ProcessThread inner class for read generation
**Usage**: Serves as a template for developing custom read generation tools in BBTools ecosystem

## A_Sample_Unpaired (A_Sample_Unpaired.java)
**Purpose**: Sample template for processing single-end sequencing reads using BBTools framework
**Core Function**: Demonstrates input/output stream processing for single-end read data with configurable parameters
**Key Features**:
- Concurrent read input and output stream management
- Configurable verbose and error handling modes
- Flexible file format detection (FASTQ default)
- Performance timing and read/base statistics tracking
**Usage**: Used as a template for developing new single-end read processing tools in BBTools

## Accumulator (Accumulator.java)
**Purpose**: Interface for accumulating statistics captured by threads
**Core Function**: Provides methods for thread-safe statistic accumulation and thread completion tracking
**Key Features**:
- Generically typed accumulation method
- Thread synchronization via ReadWriteLock
- Success status tracking
- Allows collecting results from parallel thread processing
**Usage**: Used in multi-threaded processing to gather and aggregate results from concurrent operations

## BBTool_ST (BBTool_ST.java)
**Purpose**: Abstract base class for single-threaded BBTools applications with standardized framework.
**Core Function**: Provides foundational infrastructure for command-line processing, I/O handling, and execution flow for BBTools.
**Key Features**:
- Standardized command-line argument parsing
- Consistent initialization and setup patterns
- Unified error handling and performance monitoring
- Integrated I/O and stream management
**Usage**: Serves as a template for creating single-threaded BBTools applications, enforcing consistent behavior across tools

## DoWorker (DoWorker.java)
**Purpose**: Interface defining a work execution contract for generic worker implementations
**Core Function**: Provides a standard method signature for performing work across different worker types
**Key Features**:
- Single method `doWork()` for executing worker-specific tasks
- Allows custom implementation of work logic
- Supports error handling within worker implementation
**Usage**: Used as a template for creating workers that perform specific computational or processing tasks in BBTools

## ThreadPoolJob (ThreadPoolJob.java)
**Purpose**: Generic job processing class for concurrent thread pool operations
**Core Function**: Provides a template for processing jobs with input/output queues and automatic cleanup
**Key Features**:
- Generic type support for flexible job processing
- Automatic input data return to destination queue
- Implements basic job lifecycle with doWork() method
- Built-in interrupt handling for queue operations
**Usage**: Extend this class to create specific concurrent job processing tasks in multi-threaded environments

## ThreadWaiter (ThreadWaiter.java)
**Purpose**: Utility class for managing and synchronizing thread operations
**Core Function**: Provides static methods to start, wait for, and synchronize thread execution
**Key Features**:
- Wait for threads to start using `waitForThreadsToStart()`
- Wait for threads to complete using `waitForThreadsToFinish()`
- Safely start and synchronize multiple threads
- Supports thread result accumulation with optional locking
**Usage**: Coordinates parallel thread execution in multi-threaded applications, ensuring proper thread lifecycle management