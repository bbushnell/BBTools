# fileIO Package - File Input/Output Operations
*Comprehensive utilities for efficient file reading, writing, and stream management across multiple formats*

## Package Overview
The `fileIO` package provides essential file I/O infrastructure for BBTools, offering high-performance, thread-safe file operations with support for multiple compression formats, streaming protocols, and specialized file handling capabilities.

---

## Core File Reading Classes

### ByteFile (fileIO/ByteFile.java)
**Purpose**: Abstract base class for efficient byte-level file reading and processing across different file formats and input strategies.

- **Core Function**: Provides a flexible, extensible framework for reading files as byte arrays with multiple implementation strategies and performance optimizations.

- **Key Features**:
  - Factory method-based file creation (`makeByteFile`) with dynamic selection of reading strategy
  - Support for different file input formats via `FileFormat` class
  - Ability to read files as byte lines or count lines efficiently
  - Automatic subprocess and memory-aware file reading strategies
  - Thread-aware file reading mode selection
  - Abstract methods for core file operations (reset, nextLine, close)

- **Key Methods**:
  - `makeByteFile()`: Static factory method for creating appropriate ByteFile instances
  - `toByteLines()`: Convert entire file to a list of byte arrays
  - `countLines()`: Count total number of lines in a file
  - `nextList()`: Read batches of lines (up to 200 lines)
  - `nextLine()`: Read next line as byte array

### ByteFile1 (ByteFile1.java)
**Purpose**: Efficient byte-based file reader with buffered input stream handling for line-by-line text file processing.

- **Core Function**: 
  - Provides low-level, memory-efficient line reading from input streams
  - Supports reading from files or standard input
  - Handles various line ending formats (Windows/Unix/Mac)
  - Implements dynamic buffer management for large file processing

- **Key Features**:
  - **Buffered Line Reading**: 
    - Dynamic buffer expansion (default 16KB)
    - Efficient line extraction with minimal memory overhead
    - Handles files with variable line lengths
  
  - **Stream Flexibility**:
    - Supports file and stdin input
    - Allows subprocess execution for compressed files
    - Provides `reset()` and `close()` methods for stream management

### ByteFile2 (ByteFile2.java)
**Purpose**: Multi-threaded byte file reader designed to enhance disk reading performance, particularly for compressed files.

- **Core Function**: Runs a ByteFile1 instance in a separate thread to optimize file reading operations, providing a drop-in compatible alternative to standard file reading with improved performance characteristics.

- **Key Features**:
  - **Threaded File Reading**: Utilizes a separate thread (`BF1Thread`) to read file contents concurrently
  - **Buffered I/O**: Implements a double-buffering mechanism using `ArrayBlockingQueue` for efficient file reading
  - **Configurable Parameters**:
    - `bufflen`: Maximum number of lines per buffer (default: 1000)
    - `buffs`: Number of buffers (default: 4)
    - `buffcapacity`: Maximum byte capacity per buffer (default: 256,000 bytes)

---

## Stream Writing Classes

### ByteStreamWriter (ByteStreamWriter.java)
**Purpose**: High-performance, thread-safe byte stream writing utility for handling various file output formats with buffering and ordered/unordered writing strategies.

- **Core Function**: Efficiently writes byte data to output streams using a threaded buffering mechanism, supporting multiple file formats and writing strategies.

- **Key Features**:
  - Supports multiple output formats (FASTQ, FASTA, SAM, BAM, text)
  - Threaded writing with configurable buffering
  - Ordered and unordered stream writing modes
  - Dynamic buffer management
  - Flexible print methods for various data types
  - Compressed file output support

- **Performance Considerations**:
  - Initial buffer size configurable (default 36000 bytes)
  - Maximum buffer length before flushing (default 32768 bytes)
  - Supports forced buffer flushing
  - Minimal synchronization overhead in writing process

### TextStreamWriter (TextStreamWriter.java)
**Purpose**: Threaded text output stream writer with advanced buffering and multi-format support for writing character sequences to files or standard output.

- **Core Function**: Provides a thread-safe, buffered writing mechanism for text-based output streams, supporting multiple file formats including FASTQ, FASTA, SAM, and custom text outputs.

- **Key Features**:
  - **Flexible File Output**: Supports writing to files or standard output with configurable overwrite and append modes
  - **Multi-Format Support**: Handles different output formats like FASTQ, FASTA, SAM, and custom text
  - **Buffered Writing**: Uses an `ArrayBlockingQueue` and buffer mechanism to optimize write operations
  - **Thread-Safe Design**: Extends `Thread` class with synchronized methods for safe concurrent writing
  - **Ordered Write Capability**: Provides `writeOrdered()` method to write character sequences in a specific sequence

---

## File Management and Utilities

### FileFormat (FileFormat.java)
**Purpose**: Provides comprehensive metadata and detection capabilities for various file formats and types in bioinformatics file processing.

- **Core Function**: Detect, validate, and manage metadata for input and output files across multiple bioinformatics file formats, supporting extensive format and compression type recognition.

- **Key Features**:
  - **Format Detection**: Identifies file types through extension and content analysis (FASTA, FASTQ, SAM, BAM, VCF, etc.)
  - **Compression Support**: Recognizes multiple compression formats (RAW, GZIP, ZIP, BZ2, XZ, 7z, DSRC)
  - **Interleaving Detection**: Determines whether sequence files are single-ended or interleaved
  - **Quality Offset Analysis**: Detects ASCII quality score offsets for FASTQ files
  - **File Metadata Extraction**: Provides detailed file metadata including barcode information, read length, and file characteristics

- **Notable Methods**:
  - `testInput()`: Comprehensive file format detection method
  - `testInterleavedAndQuality()`: Analyzes file content for interleaving and quality score characteristics
  - `getFirstOctet()`: Reads first 8 lines of a file for initial format detection

### ReadWrite (ReadWrite.java)
**Purpose**: Comprehensive utility class for reading, writing, and managing file input/output operations with advanced compression and multi-threading support.

- **Core Function**: Provides flexible, high-performance file reading and writing methods with support for multiple compression formats (gzip, zip, bzip2, etc.) and concurrent file operations.

- **Key Features**:
  - **Multi-Compression Format Support**: 
    - Handles reading/writing compressed files (.gz, .zip, .bz2, .xz, .dsrc)
    - Dynamic compression method selection based on file extension and available system tools
    - Supports gzip, pigz, bgzip, bzip2, and other compression utilities

  - **Flexible I/O Streams**:
    - `getInputStream()` and `getOutputStream()` methods for dynamic stream creation
    - Supports buffered and raw input/output streams
    - Handles various file sources including local files, JAR resources, and standard input/output

  - **Concurrent File Operations**:
    - Thread-safe writing methods with `writeStringInThread()` and `writeObjectInThread()`
    - Thread management for controlling active and running write operations
    - Maximum thread control with `maxWriteThreads` configuration

---

## Text File Processing

### TextFile (TextFile.java)
**Purpose**: Robust utility class for reading and processing text files with flexible input handling and advanced line processing capabilities.

- **Core Function**: Provides comprehensive text file reading operations supporting various input sources including files, standard input, and compressed file formats with advanced line-by-line processing.

- **Key Features**:
  - Supports multiple input sources: local files, standard input, JAR resources
  - Optional subprocess decompression for compressed files
  - Line-by-line reading with configurable blank line handling
  - Efficient memory management using BufferedReader
  - Static utility methods for file-to-string array conversions
  - Line counting and file existence verification

- **Advanced Processing Methods**:
  ```java
  // Convert entire file to string array
  public final String[] toStringLines()
  
  // Count total lines in file
  public final long countLines()
  
  // Static splitting utilities
  public static String[][] doublesplitTab(String[] lines, boolean trim)
  public static String[][] doublesplitWhitespace(String[] lines, boolean trim)
  ```

### GenericTextFile (GenericTextFile.java)
**Purpose**: Flexible text file handling utility with configurable parsing and processing capabilities for various text-based file formats.

- **Key Features**:
  - Generic text file operations with flexible parsing
  - Configurable processing capabilities
  - Support for various text-based file formats

---

## Specialized File Operations

### ArrayFile (ArrayFile.java)
**Purpose**: Specialized file reader for parsing array data files with metadata extraction and float array processing capabilities.

- **Key Features**:
  - Extends TextFile for base file reading capabilities
  - Specialized for reading float arrays with metadata
  - Uses assertions for strict parsing validation
  - Provides main method demonstrating basic usage

### MatrixFile (MatrixFile.java)
**Purpose**: Specialized file reader for parsing matrix data files with structured parsing and metadata extraction.

- **Core Function**: Reads matrix files with specific formatting, extracting matrix metadata and converting text representations into float[][] matrix objects.

- **Key Features**:
  - Extends TextFile for line-by-line reading of matrix files
  - Parses matrix metadata including name, dimensions, prefix, and count
  - Converts text-based matrix representations into Matrix objects
  - Skips comment lines starting with "//"
  - Handles matrix files with specific structural requirements

### SummaryFile (SummaryFile.java)
**Purpose**: Parses and validates genome summary files, extracting critical metadata about genome builds and file characteristics.

- **Core Function**: 
  - Reads and parses tab-delimited summary text files containing genome build information
  - Validates file metadata like chromosomes, bases, versions, and source characteristics
  - Provides comparison methods to check file integrity and matching

- **Key Methods**:
  - `compare(String refName)`: Validates summary file against a reference file
  - `getName()`: Generates default summary file path
  - `SummaryFile(String path)`: Constructor parsing summary file

---

## File Operation Utilities

### File Copying and Management

#### CopyFile (CopyFile.java)
**Purpose**: Utility class for file copying operations with explicit compression handling and performance tracking.

- **Core Function**: 
  - Copies files between source and destination paths with configurable compression and overwrite settings
  - Forces compression/decompression even for files with identical extensions, primarily for benchmarking purposes

- **Key Features**:
  - Synchronized file copying method with error handling
  - Configurable path creation for destination file
  - Supports overwrite and append modes
  - Performance timing and speed reporting

#### CopyFiles (CopyFiles.java)
**Purpose**: Utility class for processing and renaming files within directories, with a focus on batch file operations.

- **Core Function**: Recursively traverses directories and renames/copies files based on specific criteria, particularly targeting chromosome text files for conversion.

#### CopyFiles2 (CopyFiles2.java)
**Purpose**: Advanced file copying utility for selective batch file transfer and filtering across source directories.

- **Core Function**: 
  - Recursively copies files and directories from multiple source roots to a destination directory
  - Implements sophisticated filtering mechanisms to selectively copy files based on predefined rules
  - Supports batch file operations with configurable inclusion and exclusion criteria

### File Compression and Concatenation

#### CompressFiles (CompressFiles.java)
**Purpose**: Utility for batch file compression operations with support for multiple compression formats and recursive directory processing.

- **Key Features**:
  - Batch file compression with zip and gzip support
  - Recursive directory compression capabilities
  - Safety checks and filtering mechanisms
  - External compression utility integration

#### Concatenate (Concatenate.java)
**Purpose**: Accepts multiple input files and concatenates them into a single output file, with support for recompression and stream handling.

- **Core Function**: 
  - Sequentially reads multiple input files and writes their contents to a single output destination
  - Supports flexible file input through command-line arguments
  - Provides ability to concatenate files while avoiding standard input/output streams

---

## Threading and Chain Processing

### LoadThread (LoadThread.java)
**Purpose**: Implements a generic, thread-safe mechanism for asynchronous file loading with controlled concurrency.

- **Core Function**: Enables non-blocking file reading operations by loading files in separate threads while managing thread count and synchronization

- **Key Features**:
  - Generic type support for flexible object deserialization
  - Static thread management with active/waiting/running state tracking
  - Configurable thread limit based on system resources
  - Memory-aware thread allocation (adjusts for low-memory environments)
  - Synchronization mechanisms to prevent thread exhaustion

### PipeThread (PipeThread.java)
**Purpose**: A specialized thread for efficiently copying data between input and output streams, facilitating stream redirection and pipeline management.

- **Core Function**: Manages background stream copying, allowing non-blocking data transfer between input and output streams

- **Key Features**:
  - Efficient byte buffer-based data transfer (8196-byte buffer)
  - Thread-safe stream management with synchronization
  - Automatic stream closure for non-standard streams
  - Terminatable data pipeline
  - Error handling for stream operations

### Chain File Processing

#### ChainBlock (ChainBlock.java)
**Purpose**: Parses and manages genome build conversion chain file data from UCSC chain file format.

- **Core Function**: 
  - Represents a single chain block from genome alignment conversion files
  - Processes complex genome alignment information between reference and query sequences
  - Converts chain file entries into structured alignment chunks for cross-genome mapping

#### ChainLine (ChainLine.java)
**Purpose**: Represents individual alignment segments within genome build conversion chain files with coordinate mapping capabilities.

- **Core Function**: Manages line-level chain file data for genome coordinate conversion and alignment processing

---

## Utility Classes

### FindFiles (FindFiles.java)
**Purpose**: File searching utility with flexible pattern matching and recursive directory traversal capabilities.

- **Key Features**:
  - Support for prefix, suffix, and middle filename pattern matching
  - Case-insensitive search capabilities
  - Supports null, #, and * as wildcard patterns
  - Recursive directory traversal
  - Cross-platform path handling

### OpenFile (OpenFile.java)
**Purpose**: Utility class for opening and reading input streams with basic file handling capabilities.

- **Core Function**: Provides a demonstration of opening and reading input streams from files, with error handling for input and stream closure.

### QuickFile (QuickFile.java)
**Purpose**: Experimental file reading utility designed for performance testing and efficient line-by-line text file processing.

- **Core Function**: Low-overhead file reading mechanism with optimized buffering and line extraction capabilities

- **Key Features**:
  - Buffered line reading with dynamic buffer resizing
  - Support for reading from stdin, files, and compressed files
  - Multiple line extraction methods
  - Robust error handling and file close management
  - Configurable subprocess decompression support

---

## Package Usage
The fileIO package serves as the foundational I/O infrastructure for BBTools, providing:
- High-performance file reading and writing operations
- Multi-format and multi-compression support
- Thread-safe concurrent file operations
- Specialized bioinformatics file format handling
- Robust stream management and pipeline operations
- File utility operations for batch processing

## Performance Considerations
- Optimized buffering strategies across all file operations
- Thread-safe concurrent processing capabilities
- Memory-efficient handling of large files
- Support for compressed file formats without performance penalties
- Configurable threading limits and buffer sizes

## Dependencies
- Relies on Java NIO and standard I/O libraries
- Integrates with BBTools shared utilities and data structures
- Supports external compression utilities (pigz, bgzip, etc.)
- Compatible with multiple operating systems and file systems

---
*Documentation generated using evidence-based analysis of source code*