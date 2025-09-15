# hiseq Package - High-Throughput Sequencing Data Processing
*Specialized tools for analyzing and processing Illumina sequencing flow cell data*

## Package Overview
The `hiseq` package provides a complete suite of utilities for processing high-throughput Illumina sequencing data at the flow cell level. These tools enable quality control, format conversion, statistical analysis, and visualization of sequencing runs across different platforms including HiSeq, MiSeq, NextSeq, and NovaSeq systems.

---

## Flow Cell Analysis and Quality Control

### AnalyzeFlowCell (AnalyzeFlowCell.java)
**Purpose**: Analyzes a flow cell by identifying and filtering low-quality genomic sequencing reads

**Core Function**: Processes Illumina sequencing reads from FASTQ files, evaluating read quality across different tiles and lanes of a flow cell, and filtering out low-quality reads based on multiple statistical and quality criteria

**Key Features**:
- Supports paired-end and single-end read processing.
- Analyzes read quality using probabilistic error estimation
- Performs k-mer depth analysis using Bloom filter technique
- Supports optional read merging for paired-end reads
- Handles barcode detection and analysis
- Provides flexible quality filtering options (trimming, discarding)

**Technical Implementation**:
- Uses concurrent read input/output streams for efficient processing
- Employs multi-threaded processing with `ProcessThread` inner class
- Implements dynamic tile-based quality assessment using `FlowCell` and `MicroTile` classes
- Uses Bloom filter for k-mer counting and depth estimation
- Supports quality recalibration

**Usage/Command-Line**:
```bash
java AnalyzeFlowCell in1=input_R1.fastq in2=input_R2.fastq out1=output_R1.fastq out2=output_R2.fastq
```

**Dependencies**:
- Relies on multiple BBTools package classes: `Read`, `BloomFilter`, `SideChannel3`, `BarcodeStats`
- Utilizes Java concurrent utilities for thread-safe processing
- Depends on custom parsing and quality assessment utilities

### PlotFlowCell (PlotFlowCell.java)
**Purpose**: Analyzes and processes flow cell data in sequencing experiments, extracting quality metrics and filtering reads.

**Core Function**: Processes sequencing reads from FASTQ files, performing detailed analysis of flow cell tiles by computing k-mer statistics, barcode information, and generating quality reports. Supports multithreaded processing with optional Bloom filter-based k-mer counting.

**Key Features**:
- Parallel processing of sequencing reads using concurrent input streams
- Configurable k-mer analysis with multiple counting strategies
- Flow cell tile-level statistics generation
- Optional reference mapping for quality control
- Bloom filter-based k-mer frequency estimation
- Barcode error distance calculation

**Technical Implementation**:
- Uses concurrent read input streams for efficient processing
- Supports multiple k-mer counting modes:
  - All k-mers per read
  - Single random k-mer per read
- Bloom filter integration for memory-efficient k-mer tracking
- Tile-level statistics tracked via `MicroTile` objects
- Thread-safe accumulation of processing results

**Usage/Command-Line**:
```bash
java PlotFlowCell in=lane1.fq.gz size=999999 dump=dump.txt expected=expected.txt bloom multithreaded -Xmx40g
```

**Key Configuration Options**:
- `in`: Input FASTQ file(s)
- `size`: Tile grid size
- `dump`: Output file for flow cell statistics
- `expected`: Expected barcodes file
- `bloom`: Enable Bloom filter k-mer counting
- `multithreaded`: Enable parallel processing

**Dependencies**:
- `stream.Read`
- `fileIO.FileFormat`
- `bloom.BloomFilter`
- `barcode.BarcodeStats`
- `aligner.MicroAligner3`

---

## Format Conversion and Header Processing

### BGI2Illumina (BGI2Illumina.java)
**Purpose**: Converts BGI sequencing data file headers to Illumina format

**Core Function**: 
Processes FASTQ input files, parsing and transforming read headers from BGI sequencing format to standard Illumina header format. Supports single and paired-end read processing with flexible input/output file handling.

**Key Features**:
- Converts BGI sequencing headers to Illumina-compatible headers
- Supports single and paired-end read processing
- Handles input file path substitution using '#' placeholder
- Flexible file extension handling
- Validates input and output file accessibility
- Dynamically adjusts interleaving settings for input/output streams

**Technical Implementation**:
- Uses `BGIHeaderParser2` for header parsing and transformation
- Supports optional barcode insertion during header conversion
- Utilizes concurrent read input/output streams for efficient processing
- Tracks reads and bases processed/output with detailed metrics
- Thread-safe processing of read lists with `ConcurrentReadInputStream` and `ConcurrentReadOutputStream`

**Usage/Command-Line**:
```bash
java hiseq.BGI2Illumina in1=input_R1.fastq in2=input_R2.fastq out1=output_R1.fastq out2=output_R2.fastq
```

**Dependencies**:
- `fileIO.ByteFile`
- `fileIO.FileFormat`
- `fileIO.ReadWrite`
- `shared.Parse`
- `shared.Parser`
- `stream.ConcurrentReadInputStream`
- `stream.ConcurrentReadOutputStream`
- `stream.FASTQ`
- `stream.Read`

---

## Header Parsers

### ReadHeaderParser (ReadHeaderParser.java)
**Purpose**: Abstract base class for parsing Illumina sequencing read headers across different sequencing platforms.

**Core Function**: Provides a standardized framework for extracting metadata from sequencing read headers, including machine information, run details, lane, tile, coordinates, and barcode information. This is an abstract class that defines a contract for specific header parser implementations.

**Key Features**:
- **Header Parsing**: Extracts detailed metadata from sequencing read identifiers
- **Platform Flexibility**: Supports multiple Illumina platforms (MiSeq, HiSeq, NextSeq, NovaSeq)
- **Coordinate Extraction**: Parses tile, x/y positions, surface, and swath information
- **Barcode Handling**: Supports complex barcode parsing with delimiter and length detection
- **Pair Information**: Extracts pair number and chastity filtering details

**Technical Implementation**:
- Abstract methods enforce consistent parsing across different sequencing platforms
- Uses static utility methods for coordinate and barcode parsing
- Supports configurable parsing via static flags `PARSE_COORDINATES` and `PARSE_COMMENT`
- Provides default implementations for common parsing tasks
- Example header parsing patterns supported:
  - `/VP2-06:112:H7LNDMCVY:2:2437:14181:20134`
  - `MISEQ08:172:000000000-ABYD0:1:1101:18147:1925`
  - `A00178:38:H5NYYDSXX:2:1101:3007:1000`

**Usage/Command-Line**:
```java
// Typical usage in a subclass
public class SpecificHeaderParser extends ReadHeaderParser {
    @Override
    public ReadHeaderParser parse(String id) {
        // Implement platform-specific parsing logic
    }
    // Implement other abstract methods
}

// Test method for header parsing
ReadHeaderParser parser = new SpecificHeaderParser();
parser.test(readHeader); // Prints detailed metadata
```

**Dependencies**:
- `shared.Tools`: Utility methods for character classification
- `stream.Read`: Read object for parsing
- `structures.ByteBuilder`: Used for toString() representation

### IlluminaHeaderParser1 (IlluminaHeaderParser1.java)
**Purpose**: Parses Illumina sequencing machine read headers to extract positional and metadata information.

**Core Function**: A specialized parser that processes Illumina sequencing machine headers, extracting detailed information about each sequenced read, including lane, tile, coordinates, pair number, and barcode details. This class is designed to handle various Illumina sequencing platform header formats (MiSeq, HiSeq, NovaSeq).

**Key Features**:
- Extracts precise sequencing read metadata from header strings
- Supports parsing of complex header formats from multiple Illumina platforms
- Provides methods to retrieve specific header components like lane, tile, x/y coordinates
- Handles different header separators (' ' and '/')
- Robust error handling with try-catch mechanism

**Technical Implementation**:
- Uses recursive parsing methods to extract header components
- Key parsing methods:
  - `parseCoordinates()`: Extracts lane, tile, x, and y coordinates
  - `parseComment()`: Extracts pair code, chastity code, control bits, barcode
  - `findCommentSeparator()`: Locates header comment separator
- Supports parsing integer and character tokens with `parseInt()` and `parseChar()` methods
- Flexible parsing strategy accommodating variations in header formats

**Usage/Command-Line**:
```java
// Example usage in main method
IlluminaHeaderParser1 ihp = new IlluminaHeaderParser1();
ihp.parse("@A00178:38:H5NYYDSXX:2:1101:3007:1000 1:N:0:CAACCTA+CTAGGTT");
int lane = ihp.lane();
String barcode = ihp.barcode();
```

**Dependencies**:
- `shared.KillSwitch`
- `shared.Tools`
- `structures.ByteBuilder`
- Extends `ReadHeaderParser`

### IlluminaHeaderParser2 (IlluminaHeaderParser2.java)
**Purpose**: Faster parsing of Illumina sequencing read headers, extracting detailed machine-specific metadata.

**Core Function**: Provides a high-performance implementation for parsing and extracting structured information from complex Illumina sequencing read identifiers, supporting multiple machine types and header formats through efficient tokenization.

**Key Features**:
- Extends `ReadHeaderParser` for standardized header processing
- Supports parsing headers from various Illumina platforms (NovaSeq, HiSeq, MiSeq, NextSeq)
- Extracts precise metadata including:
  - Machine identifier
  - Run number
  - Flowcell details
  - Lane, tile, and coordinate information
  - Pair and chastity codes
  - Barcode and index information
- Provides coordinate encoding method
- Uses `LineParserS3` for efficient parsing with configurable delimiter

**Technical Implementation**:
- Uses delimiter-based parsing with `:` separator
- Implements coordinate encoding via bitwise operations
- Validates header structure through methods like `looksValid()` and `looksShrunk()`
- Supports flexible header parsing with configurable whitespace index tracking

**Usage/Command-Line**:
```java
IlluminaHeaderParser2 parser = new IlluminaHeaderParser2();
parser.parse("@LH00223:28:22GLGMLT3:1:1101:5928:1016 1:N:0:CTGCTTGGTT+CTAACGACAG");
String machine = parser.machine(); // Returns machine identifier
int lane = parser.lane(); // Returns lane number
```

**Dependencies**:
- `shared.LineParser`
- `shared.LineParserS3`
- `structures.ByteBuilder`
- Extends `ReadHeaderParser`

### BGIHeaderParser (BGIHeaderParser.java)
**Purpose**: Parses header metadata for BGI (Beijing Genomics Institute) sequencing data files

**Core Function**: Extends ReadHeaderParser to extract specific metadata components from BGI sequencing read headers, including parsing machine, run, flowcell, lane, tile, and other sequencing parameters

**Key Features**:
- Supports parsing complex BGI header formats like `v300056266_run28L3C001R0010057888/1`
- Converts BGI headers to Illumina-style format via `toIllumina()` method
- Uses `LineParserS4` to tokenize header strings with `_LCR/\t` delimiters
- Extracts specific header components: sample, flowcell, lane, tile, x/y positions, pair code

**Technical Implementation**:
- Inherits from `ReadHeaderParser`
- Uses `LineParserS4` for flexible string parsing
- `ByteBuilder` used for efficient string construction in `toIllumina()` method
- Multiple method overrides to extract different header components
- Supports handling headers with variable number of terms

**Dependencies**:
- `shared.LineParserS4`
- `structures.ByteBuilder`
- Extends `ReadHeaderParser`

**Usage Example**:
```java
BGIHeaderParser parser = new BGIHeaderParser();
parser.parse("v300056266_run28L3C001R0010057888/1");
String illuminaHeader = parser.toIllumina(null);
```

### BGIHeaderParser2 (BGIHeaderParser2.java)
**Purpose**: Parses BGI sequencing run headers using a reverse parser to extract metadata from complex file naming conventions.

**Core Function**: Transforms BGI sequencing header identifiers into standardized header formats, supporting conversion to Illumina-style headers with metadata extraction from complex naming patterns.

**Key Features**:
- Reverse parsing of BGI header identifiers using custom delimiter "_LCR/"
- Supports extracting machine, flowcell, lane, tile, x/y positions, and pair information
- Converts BGI headers to Illumina-style headers with optional barcode insertion
- Optional extra information parsing with whitespace detection
- Handles headers with format: `v300056266_run28L3C001R0010057888/1`

**Technical Implementation**:
- Uses `LineParserS4Reverse` for parsing complex header strings
- Implements header parsing with position-based token extraction
- Uses `ByteBuilder` for efficient string construction in `toIllumina()` method
- Static `PARSE_EXTRA` flag controls optional extra information parsing

**Usage/Command-Line**:
```java
BGIHeaderParser2 parser = new BGIHeaderParser2();
parser.parse("v300056266_run28L3C001R0010057888/1");
String illuminaHeader = parser.toIllumina("ACGTACGT");
```

**Dependencies**:
- `shared.LineParserS4Reverse`: Custom line parsing utility
- `structures.ByteBuilder`: Efficient string builder
- Extends `ReadHeaderParser` abstract class

---

## Flow Cell Data Structures

### FlowCell (FlowCell.java)
**Purpose**: A data structure representing an Illumina sequencing flow cell, managing lane and tile-level statistics for genomic read processing.

**Core Function**: Manages complex data aggregation and statistical computation for sequencing reads across multiple lanes and tiles, including read alignment, error rates, and quality metrics. Provides methods for processing, widening, and analyzing microtile-level genomic data.

**Key Features**:
- Dynamic lane and microtile management
- Comprehensive read statistics calculation
  - Tracks reads processed, aligned, and error rates
  - Calculates average and standard deviation for multiple quality metrics
- Adaptive data widening method for read targeting
  - Can expand tile dimensions to reach target read counts
- Statistical error rate computation
  - Calculates alignment and base error rates
- Supports loading from tile dump files
- Provides microtile retrieval methods with optional creation

**Technical Implementation**:
- Uses `ArrayList<Lane>` for storing sequencing lane data
- Implements methods for generating statistical points using `Tools.linearRegression()`
- Supports blurring/smoothing of microtile data
- Utilizes synchronized blocks for thread-safe lane manipulation

**Usage/Command-Line**: 
```java
// Create a FlowCell from a dump file
FlowCell flowCell = new FlowCell("sequencing_data.dump");

// Calculate statistics
flowCell.calcStats();

// Widen to target reads
FlowCell processedCell = flowCell.widenToTargetReads(10000);

// Get alignment rate
double alignmentRate = flowCell.alignmentRate();
```

**Dependencies**:
- `java.util.ArrayList`
- `java.util.Collections`
- `shared.Tools`
- `structures.FloatList`
- `structures.Point`
- Custom classes: `MicroTile`, `Lane`, `TileDump`, `IlluminaHeaderParser2`

### Lane (Lane.java)
**Purpose**: Represents a sequencing lane in a high-throughput sequencing dataset, managing tiles, microtiles, and associated genomic metrics.

**Core Function**: The Lane class serves as a container for sequencing data, organizing information at multiple granularities (lane, tile, microtile) and providing methods for data aggregation, iteration, and high-depth genomic kmer analysis.

**Key Features**:
- **Tile Management**: 
  - Dynamic tile storage using `ArrayList<Tile>`
  - Lazy tile initialization with `getTile()` method
  - Thread-safe tile addition with synchronized blocks

- **Genomic Metrics Calculation**: 
  - Calculates high-depth genomic kmer rate via `calcHighDepthGenomic()`
  - Tracks hits, misses, errors, and aligned bases across microtiles
  - Computes error rates and kmer correctness probabilities

- **Parallel Data Structures**:
  - Uses `AtomicLongArray` for thread-safe depth, match, and count tracking
  - Supports multiple depth-related metrics across two dimensions

- **Iterator Support**:
  - Implements `Iterable<Tile>` interface
  - Custom `TileIterator` for efficient tile traversal

**Technical Implementation**:
- **Genomic Kmer Rate Calculation**:
  - Error rate (E) = base errors / aligned bases
  - Per-base correctness (C) = 1 - E
  - Kmer correctness probability (P) = C^k
  - Unique kmer fraction (U) = misses / (hits + misses)
  - High-depth genomic kmer fraction constrained between 0.0001 and 0.9999

**Usage**:
```java
Lane lane = new Lane(1);  // Create lane with index 1
MicroTile microTile = lane.getMicroTile(tileIndex, x, y);  // Access specific microtile
double highDepthRate = lane.calcHighDepthGenomic(k);  // Compute high-depth genomic kmer rate
```

**Dependencies**:
- `java.util.ArrayList`
- `java.util.Iterator`
- `java.util.concurrent.atomic.AtomicLongArray`
- Local classes: `Tile`, `MicroTile`
- Utility classes: `fileIO.ByteStreamWriter`, `shared.Tools`

### Tile (Tile.java)
**Purpose**: Two-dimensional grid-based container for MicroTile objects representing spatial data in a sequencing lane and tile.

**Core Function**: Manages a hierarchical 2D spatial data structure that allows dynamic creation and retrieval of MicroTile instances organized in a grid-like arrangement, with support for iteration and text serialization.

**Key Features**:
- **Dynamic MicroTile Management**: 
  - Create MicroTiles on-demand using `get(int x, int y, boolean create)` method
  - Supports lazy initialization of grid elements
- **Spatial Organization**: 
  - Fixed grid size with static `xSize` and `ySize` of 500
  - Organizes MicroTiles into nested ArrayList structure
- **Iteration Support**: 
  - Implements `Iterable<MicroTile>` interface
  - Provides `iterator()` method for traversing all MicroTiles
- **Text Serialization**: 
  - `toText()` method converts MicroTile data to ByteBuilder
  - Supports optional parameters for text representation

**Technical Implementation**:
- Uses nested `ArrayList<ArrayList<MicroTile>>` for 2D grid representation
- Synchronized methods for thread-safe MicroTile addition
- Grid indices calculated by integer division of coordinates
- Lazy grid expansion with `getIndex()` method

**Usage**:
```java
Tile tile = new Tile(laneNumber, tileNumber);
MicroTile microTile = tile.get(x, y, true); // Create or retrieve MicroTile
Iterator<MicroTile> iterator = tile.iterator(); // Iterate through MicroTiles
```

**Dependencies**:
- `java.util.ArrayList`
- `java.util.Iterator`
- `structures.ByteBuilder`
- Depends on `MicroTile` class from the same package

### MicroTile (MicroTile.java)
**Purpose**: Represents a microchip tile in high-throughput sequencing, tracking detailed read and base-level metrics.

**Core Function**: Collects and manages statistical information about sequencing reads from a specific tile, including quality metrics, error rates, base composition, and alignment characteristics.

**Key Features**:
- Tracks per-tile sequencing metrics like read count, base count, error rates
- Calculates quality statistics including:
  - Average read quality
  - Expected base error rates
  - Alignment rates
  - K-mer error rates
- Supports adding read-level metrics through `add()` methods
- Generates detailed text output of tile statistics
- Provides comparison and sorting of tiles based on quality metrics

**Technical Implementation**:
- Uses specialized data tracking for:
  - Base composition tracking (`acgtn` array)
  - Homopolymer G-track
  - Kmer quality calculations
- Supports flexible metric accumulation via `add()` and `multiplyBy()` methods
- Implements `Comparable` for tile quality sorting
- Optional cycle tracking via `CycleTracker`

**Usage/Command-Line**: 
- Typically used within sequencing pipeline to aggregate and analyze tile-level read statistics
- Can generate detailed statistical headers and text representations of tile data

**Dependencies**:
- `QualityTools`
- `AminoAcid`
- `Tools`
- `Read` and `SamLine` from `stream` package
- `ByteBuilder` from `structures` package

---

## Specialized Tools

### FlowcellCoordinate (FlowcellCoordinate.java)
**Purpose**: Parses and represents coordinate information for Illumina sequencing flowcell data points.

**Core Function**: Extracts and manages lane, tile, x, and y coordinates from Illumina sequencing read identifiers using an IlluminaHeaderParser2. Provides thread-safe coordinate parsing and comparison capabilities for sequencing data processing.

**Key Features**:
- Parses complex Illumina sequencing identifiers into discrete coordinate components
- Implements `Comparable` interface for sorting and comparing flowcell coordinates
- Supports optional UMI (Unique Molecular Identifier) parsing via `parseUMI` flag
- Thread-safe coordinate retrieval using `ThreadLocal` mechanism in `getFC()` method
- Validates coordinate state through `isSet()` method

**Technical Implementation**:
- Stores coordinates as public integer fields: `lane`, `tile`, `x`, `y`
- Uses `IlluminaHeaderParser2` for robust identifier parsing
- Comparison method `compareTo()` prioritizes coordinates in order: lane, tile, y, x
- Default constructor and constructor accepting an identifier string

**Usage**:
```java
// Parse a flowcell identifier
FlowcellCoordinate fc = new FlowcellCoordinate("HISEQ07:419:HBFNEADXX:1:1101:1238:2072");
if(fc.isSet()) {
    int lane = fc.lane;    // Access parsed lane
    int tile = fc.tile;    // Access parsed tile
}

// Thread-safe retrieval
FlowcellCoordinate threadSafeFC = FlowcellCoordinate.getFC();
```

**Dependencies**:
- `hiseq.IlluminaHeaderParser2`: Primary parsing utility for sequencing identifiers
- `java.lang.ThreadLocal`: Enables thread-local storage of coordinate instances

### CycleTracker (CycleTracker.java)
**Purpose**: Tracks base and quality composition across sequencing cycles for DNA reads.

**Core Function**: Accumulates statistics about nucleotide bases (A, C, G, T, N) and their quality scores across different read positions, enabling comprehensive per-cycle base and quality distribution analysis.

**Key Features**:
- Tracks base frequencies for A, C, G, T, N at each read position
- Aggregates quality scores per read cycle
- Dynamically resizes internal storage to accommodate reads of varying lengths
- Computes cycle-wise base frequencies and maximum base occurrence
- Supports adding individual reads and merging data from multiple CycleTracker instances

**Technical Implementation**:
- Uses a 2D long array `acgtnq` to store base and quality counts
- Converts bases to numeric indices using `AminoAcid.baseToNumberACGTN`
- Calculates cycle averages by dividing base counts by total cycle counts
- Computes overall base and quality distribution statistics

**Usage**:
```java
CycleTracker tracker = new CycleTracker();
tracker.add(read);  // Add individual reads
tracker.process();  // Compute statistics
float maxA = tracker.max('A');  // Get max frequency for A
float avgC = tracker.avg('C');  // Get average frequency for C
```

**Dependencies**:
- `dna.AminoAcid`: Base-to-number conversion
- `shared.Tools`: Maximum calculation utility
- `shared.Vector`: Sum calculation
- `stream.Read`: Input read representation
- `java.util.Arrays`: Array manipulation

---

## Data Analysis and Visualization

### TileDump (TileDump.java)
**Purpose**: Processes and analyzes sequencing tile data from flow cell experiments, providing detailed statistical reporting and filtering of microtiles.

**Core Function**: Parses tile dump files from high-throughput sequencing runs, extracting detailed statistics about read quality, error rates, and tile characteristics across different lanes and tiles.

**Key Features**:
- Loads tile dump files with multiple version support
- Calculates detailed microtile statistics including:
  - Read counts and alignment rates
  - Error rates for reads and bases
  - Quality metrics like error-free percentage
- Implements advanced tile filtering mechanisms
  - Discards microtiles based on multiple quality thresholds
  - Supports configurable discard parameters like quality deviations, error rates
- Writes comprehensive output with tile and flow cell metadata
- Supports command-line configuration of analysis parameters

**Technical Implementation**:
- Uses custom parsing utilities (LineParser1, LineParser2) for efficient file reading
- Supports multiple dump file versions with version-specific loading logic
- Calculates statistical metrics using deviation and fraction-based thresholds
- Uses AtomicLongArray for thread-safe statistical tracking

**Usage/Command-Line**:
```bash
java TileDump in=input.dump out=output.txt x=100 y=100 \
  targetreads=1000 targetalignedreads=250 \
  qdeviations=2.4 udeviations=1.5
```

**Dependencies**:
- fileIO package: ByteFile, ByteStreamWriter
- shared package: LineParser, Timer, Tools
- align2 package: QualityTools
- java.util: ArrayList, Collections
- java.util.concurrent: AtomicLongArray

### PlotHist (PlotHist.java)
**Purpose**: Generates histograms from data matrix files, specifically designed for plotting TileDump data.

**Core Function**: 
PlotHist processes input text files containing numerical data, creating histogram distributions by binning values across multiple columns. It dynamically calculates bin ranges based on the maximum values in each column and generates output TSV files representing the histogram counts.

**Key Features**:
- Dynamic binning of numerical data with configurable bin count (default 1000 bins)
- Supports multi-column input files with tab-separated values
- Automatically calculates maximum values for each column to normalize histogram ranges
- Generates separate histogram output files for each column
- Handles compressed input files (.gz, .bz2) automatically
- Configurable verbosity and file overwrite options

**Technical Implementation**:
- Uses `LineParser1` for parsing tab-separated input lines
- Stores data in `long[][] countMatrix` for histogram counts
- Calculates bin indices using `(d/maxArray[term])*bins)` formula
- Writes histogram data using `ByteStreamWriter` to individual TSV files

**Usage/Command-Line**:
```bash
java hiseq.PlotHist in=input.txt bins=1000
```
- `in`: Input text file (required)
- `bins`: Number of histogram bins (optional, default 1000)
- `verbose`: Enable verbose output (optional)

**Dependencies**:
- `fileIO`: ByteFile, ByteStreamWriter, FileFormat, ReadWrite
- `shared`: LineParser1, Parse, Tools
- `structures`: ByteBuilder

### PlotReadPosition (PlotReadPosition.java)
**Purpose**: Generates a plot of read positions from Illumina sequencing data, analyzing barcode and header information.

**Core Function**: Processes input FASTQ files to extract and map read positions, barcode information, and Hamming distance between expected and actual barcodes. Outputs a tab-delimited file with x-coordinate, y-coordinate, and barcode distance.

**Key Features**:
- Parses Illumina sequencing read headers to extract x and y coordinates
- Finds closest matching barcode using Hamming distance algorithm
- Supports concurrent input stream processing for efficient read handling
- Generates output with read position and barcode distance metrics
- Tracks reads processed and processing time

**Technical Implementation**:
- Uses `PCRMatrixHDist` for barcode distance calculations
- Employs `ConcurrentReadInputStream` for parallel read processing
- Utilizes `ByteBuilder` for efficient string construction
- Supports both single-end and paired-end sequencing data

**Usage/Command-Line**:
```bash
java hiseq.PlotReadPosition in=input.fastq out=output.header expected=barcodes.txt
```

**Dependencies**:
- `barcode.Barcode`
- `barcode.BarcodeCounter`
- `barcode.PCRMatrixHDist`
- `fileIO.FileFormat`
- `stream.ConcurrentReadInputStream`
- `stream.Read`
- `shared.Parser`

---

## Package Usage

The hiseq package provides essential functionality for:

**Illumina Sequencing Analysis**:
- Quality control and filtering of flow cell data
- Per-tile and per-lane statistical analysis
- Read quality assessment across sequencing cycles

**Format Conversion**:
- BGI to Illumina header format conversion
- Cross-platform sequencing data compatibility
- Standardized header parsing for multiple platforms

**Flow Cell Management**:
- Hierarchical data organization (FlowCell → Lane → Tile → MicroTile)
- Spatial coordinate tracking and analysis
- Thread-safe tile-level statistics aggregation

**Quality Control Pipelines**:
- K-mer based quality assessment
- Bloom filter integration for memory efficiency  
- Barcode validation and error distance calculation

**Data Visualization**:
- Flow cell heatmap generation
- Histogram plotting for statistical distributions
- Read position mapping and clustering analysis

## Performance Considerations

- Multi-threaded processing support for large datasets
- Memory-efficient Bloom filter implementations
- Concurrent stream processing for I/O optimization
- Thread-safe data structures for parallel analysis
- Configurable processing parameters for resource management

## Dependencies

The hiseq package integrates with several BBTools components:
- `stream` package for read processing and I/O
- `fileIO` package for efficient file handling
- `shared` package for parsing and utility functions
- `bloom` package for memory-efficient k-mer counting
- `barcode` package for barcode analysis
- `structures` package for specialized data structures

---
*Documentation generated using evidence-based analysis of source code*