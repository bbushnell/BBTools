# gff Package - Genomic Feature Format Processing Tools
*Comprehensive suite of utilities for parsing, converting, and analyzing genomic annotation files*

## Package Overview
The `gff` package provides essential tools for working with genomic feature formats including GFF3, GTF, VCF, and GenBank flat files. These utilities enable format conversion, annotation comparison, feature extraction, and specialized processing for genomic research and bioinformatics pipelines.

---

## File Format Comparison and Analysis

### CompareGff (CompareGff.java)
**Purpose**: Compares GFF (Genomic Feature Format) files to evaluate gene-calling accuracy.

**Core Function**: 
The CompareGff class performs precision and recall analysis on genomic annotation files by comparing a reference GFF file against a query GFF file. It calculates true positives, false positives, and false negatives for gene feature starts and stops across different genomic feature types like CDS, rRNA, and tRNA.

**Key Features**:
- Supports comparing GFF files with feature types CDS, rRNA, and tRNA.
- Tracks multiple performance metrics including:
  - True Positive Start/Stop counts (ref-relative and query-relative)
  - False Positive Start/Stop counts (ref-relative and query-relative)
  - False Negative Start/Stop counts
- Generates a Signal-to-Noise Ratio (SNR) for gene prediction performance.
- Handles compressed input files with automatic .gz/.bz2 extension handling.
- Supports configurable line processing limits.

**Technical Implementation**:
- Uses `ByteFile` for efficient file reading.
- Builds hash maps to track reference and query line matches:
  - `lineMap`: Maps sequence and stop position to reference GFF lines
  - `startCountMap`: Tracks start position match counts
  - `stopCountMap`: Tracks stop position match counts
- Performs feature type filtering using `ProkObject.processType()`
- Compares query lines against reference lines, checking:
  - Stop position match
  - Strand consistency
  - Feature type consistency

**Usage/Command-Line**:
```bash
java gff.CompareGff ref=reference.gff query.gff [options]
```
Options include:
- `ref=`: Reference GFF file path
- `lines=`: Maximum number of lines to process
- `verbose=`: Enable verbose output mode

**Dependencies**:
- `fileIO.ByteFile`
- `fileIO.FileFormat`
- `prok.ProkObject`
- `shared.Tools`
- `shared.Parser`
- `structures.StringNum`

### CompareGff_old (CompareGff_old.java)
**Purpose**: Legacy tool for comparing GFF (Gene Feature Format) files to evaluate gene-calling accuracy.

**Core Function**: Performs detailed comparison between reference and query GFF files, focusing on CDS (Coding Sequence) features, and calculates various statistical metrics about gene start and stop positions.

**Key Features**:
- Processes GFF files line-by-line using ByteFile reading
- Supports compressed input files via automatic extension handling
- Generates comparative metrics for gene start/stop positions including:
  - True Positive Start/Stop counts (reference and query-relative)
  - False Positive Start/Stop counts (reference and query-relative)
  - Signal-to-Noise Ratio (SNR) calculation
- Configurable via command-line arguments with support for:
  - Reference file specification
  - Maximum lines processing
  - Verbose mode

**Technical Implementation**:
- Uses HashSet and HashMap for efficient position tracking
- Supports strand-aware position comparison
- Calculates metrics by comparing reference and query CDS feature positions
- Implements a two-pass comparison algorithm with mutable sets

**Usage/Command-Line**:
```bash
java gff.CompareGff_old input.gff ref.gff
# Optional parameters:
# lines=1000 (limit processing)
# verbose=true (enable detailed logging)
```

**Dependencies**:
- fileIO package: ByteFile, FileFormat, ReadWrite
- shared package: Parse, Parser, PreParser, Shared, Timer, Tools
- Custom GffLine class for parsing GFF entries

---

## Feature Extraction and Processing

### CutGff (CutGff.java)
**Purpose**: GFF (Gene Feature Format) file processing utility for extracting and filtering genomic feature annotations.

**Core Function**: 
Processes GFF and FASTA files to extract, filter, and transform genomic feature annotations based on multiple configurable criteria. Supports single-threaded and multi-threaded processing of multiple input files.

**Key Features**:
- Extract genomic features by type (e.g., CDS - Coding Sequences)
- Filter features based on:
  - Attribute presence/absence
  - Length constraints
  - Feature type selection
- Optional taxonomic renaming of features
- Support for inverting feature selection
- Handles compressed (.gz) input files
- Multi-threaded processing for improved performance

**Technical Implementation**:
- Uses `HashMap<String, Read>` for rapid read lookup
- Implements custom attribute filtering with `hasAttributes()` method
- Supports optional rRNA alignment and strand verification
- Flexible taxonomic identification modes:
  - Accession mode
  - GI (GenBank Identifier) mode
  - Header mode
  - Taxid mode

**Usage/Command-Line**:
```bash
java CutGff in=input.fna gff=input.gff types=CDS attributes=gene_name maxlen=1000
```

**Dependencies**:
- BBTools packages: fileIO, shared, stream, tax
- Java standard libraries: java.io, java.util
- Requires GffLine and Read classes from BBTools framework

### CutGff_ST (CutGff_ST.java)
**Purpose**: Command-line tool for processing genomic feature files (GFF) and associated FASTA sequences by extracting or modifying regions based on specific attributes.

**Core Function**: Filters and extracts genomic regions from FASTA files using corresponding GFF annotations, supporting strand-specific processing, attribute-based filtering, and optional region inversion.

**Key Features**:
- Extract genomic regions matching specific feature types (default: "CDS")
- Filter regions by attribute presence or absence
- Support for processing both forward and reverse strand annotations
- Configurable minimum and maximum region lengths
- Optional region inversion (replace region with 'N' bases)
- Handles gzipped input files automatically

**Technical Implementation**:
- Uses `GffLine` class for parsing genomic feature annotations
- Processes files pair-wise (matching FASTA and GFF files)
- Strand-specific processing with coordinate adjustment
- Uses `HashMap<String, Read>` for efficient sequence lookup
- Supports multi-threaded file reading via `ByteFile` configurations

**Usage/Command-Line**:
```bash
java CutGff_ST in=sequences.fna gff=annotations.gff types=CDS minlen=10 maxlen=1000 attributes=gene_name banattributes=partial=true
```

**Dependencies**:
- `fileIO` package: File and stream handling
- `stream` package: Read and sequence processing
- `shared` package: Utility classes and parsing
- `prok` package: Static parsing tools

---

## GenBank Format Processing

### GbffFile (GbffFile.java)
**Purpose**: Parses GenBank flat file (.gbff) format and converts it to GFF (General Feature Format) annotation files.

**Core Function**: Provides a sequential file parser for GenBank files that extracts genomic feature annotations from LOCUS records, transforming them into a standardized GFF3 output format. The parser reads GenBank files line-by-line, processing each locus record and generating corresponding GFF annotations.

**Key Features**:
- Converts GenBank (.gbff) files to GFF3 format
- Supports command-line usage for direct file conversion
- Sequentially processes genomic locus records
- Generates GFF3 headers with BBTools version information
- Thread-safe file processing with synchronized reset method

**Technical Implementation**:
- Uses `ByteFile` for low-level file reading
- Implements sequential parsing with byte-level processing
- Handles multi-line GenBank record parsing
- Detects LOCUS record boundaries using byte-level comparisons
- Skips sequence data sections during parsing

**Usage/Command-Line**:
```bash
java gff.GbffFile input.gbff [output.gff]
# Examples:
# java gff.GbffFile genome.gbff
# java gff.GbffFile genome.gbff annotations.gff
```

**Dependencies**:
- `fileIO.ByteFile`
- `fileIO.ByteStreamWriter`
- `fileIO.FileFormat`
- `shared.Shared`
- `shared.Tools`
- Custom `GbffLocus` class for record processing

### GbffLocus (GbffLocus.java)
**Purpose**: Parser and processor for GenBank locus records, extracting genomic annotation details from text-based file formats.

**Core Function**: Systematically parses different sections of a GenBank locus record, including metadata, sequence information, and genomic features. Processes records line-by-line, identifying and extracting specific annotations such as accession, organism, species, and feature details.

**Key Features**:
- Parses multiple GenBank record sections including LOCUS, DEFINITION, ACCESSION, VERSION, SOURCE, and FEATURES
- Extracts genomic metadata like accession number, organism name, and species
- Supports feature parsing for genomic elements like CDS, tRNA, and rRNA
- Converts parsed locus information to GFF (General Feature Format) output
- Handles multi-line record parsing with advanced line advancement techniques

**Technical Implementation**:
- Uses `ArrayList<byte[]>` for efficient text line processing
- Implements block-based parsing with state tracking via `num` variable
- Utilizes utility methods like `nextLine()`, `advanceBlock()`, and `toFeatureType()` for flexible parsing
- Creates `GbffFeature` objects to represent individual genomic features
- Filters features based on type and error status before GFF conversion

**Usage**:
```java
ArrayList<byte[]> gbffLines = // load GenBank file lines
GbffLocus locus = new GbffLocus(gbffLines);
locus.toGff(byteStreamWriter); // convert to GFF format
```

**Dependencies**:
- `fileIO.ByteStreamWriter` for output streaming
- `shared.Tools` for string and array utility methods
- Internal `GbffFeature` class for feature representation

### GbffFeature (GbffFeature.java)
**Purpose**: Parses and represents GenBank feature annotations, converting complex feature data into a standardized GFF (Generic Feature Format) representation.

**Core Function**: Processes raw GenBank feature annotation lines, extracting biological feature details like coordinates, strand information, and qualifiers. Transforms these annotations into a consistent, machine-readable format suitable for genomic analysis.

**Key Features**:
- Parses complex coordinate specifications for genomic features
- Supports multiple feature types: gene, CDS, rRNA, tRNA, ncRNA, etc.
- Handles strand-specific information (complement and join annotations)
- Converts GenBank annotations to GFF format with precise attribute extraction
- Supports pseudo-gene identification
- Extracts feature details like product, locus tag, and subtypes

**Technical Implementation**:
- Uses byte array processing for memory efficiency
- Implements custom parsing methods for coordinate extraction:
  - `parseStartStop()`: Parses complex coordinate strings
  - `fixLines()`: Handles multi-line feature annotations
  - `append()`: Processes feature line continuations
- Supports genomic feature types via static type arrays
- Uses `ByteBuilder` for efficient string manipulation

**Usage/Command-Line**: 
```java
// Example constructor usage
ArrayList<byte[]> featureLines = // ... get feature lines
String type = "gene";
String accession = "NC_000913.3";
GbffFeature feature = new GbffFeature(featureLines, type, accession);

// Convert to GFF format
ByteStreamWriter writer = // ... initialize writer
feature.toGff(writer);
```

**Dependencies**:
- `java.util.ArrayList`
- `java.util.Arrays`
- `fileIO.ByteStreamWriter`
- `shared.Shared`
- `shared.Tools`
- `structures.ByteBuilder`

---

## Core Format Parsers

### GffLine (GffLine.java)
**Purpose**: Parses and represents GFF3 format annotation lines for genomic features, supporting multiple input formats and providing flexible feature processing.

**Core Function**: 
Implements a comprehensive parser for GFF (General Feature Format) files, capable of converting between different genomic annotation formats like GTF, VCF, and ORF predictions. Provides methods to load, filter, and manipulate genomic feature annotations with high precision.

**Key Features**:
- Supports parsing GFF3, GTF, VCF, and ORF prediction formats
- Flexible constructors for multiple input types:
  - Byte array input
  - GTF line conversion
  - VCF variant conversion
  - Variant object conversion
  - Prokaryotic ORF prediction conversion
- Static methods for file loading and filtering
- Range mapping for genomic feature overlap queries
- Prokaryotic feature type classification

**Technical Implementation**:
- Implements `Comparable<GffLine>`, `Feature`, and `Cloneable` interfaces
- Byte-based parsing with precise coordinate extraction
- Supports 1-based coordinate system
- Strand orientation handling (0=+, 1=-, 2=?, 3=.)
- Handles missing fields with '.' placeholder
- Interning of string fields for memory efficiency

**Usage/Command-Line**:
```java
// Load GFF file, filtering by feature types
ArrayList<GffLine> features = GffLine.loadGffFile("annotations.gff", "CDS,tRNA", false);

// Generate GFF header with variant calling statistics
String header = GffLine.toHeader(properPairRate, qualityAvg, mapqAvg, rarity, 
                                  minAlleleFraction, ploidy, reads, pairs, 
                                  properPairs, bases, refGenome);
```

**Dependencies**:
- `dna.AminoAcid`
- `dna.Data`
- `fileIO.ByteFile`
- `fileIO.FileFormat`
- `ml.CellNet`
- `prok.Orf`
- `prok.ProkObject`
- `shared.Parse`
- `shared.Shared`
- `shared.Tools`
- `structures.ByteBuilder`
- `structures.Feature`
- `structures.Range`
- `var2.ScafMap`
- `var2.VCFLine`
- `var2.Var`

### GtfLine (GtfLine.java)
**Purpose**: Parses GTF (Gene Transfer Format) file lines and converts them to GFF (Gene Feature Format) representation

**Core Function**: Converts GTF file entries into a structured format by parsing tab-separated values from each line, extracting genomic feature information such as sequence name, source, start and end positions, strand, and attributes.

**Key Features**:
- Line parsing using `LineParser1` for tab-separated values
- Extracts key genomic feature components:
  - Sequence name
  - Source
  - Feature type
  - Start and end positions
  - Score handling with `.` as default
  - Strand information
  - Frame/phase parsing
  - Attribute extraction

**Technical Implementation**:
- Uses `LineParser1` for robust tab-separated value parsing
- Handles missing or default values by converting to -1 (score, frame)
- Supports byte array input for efficient parsing
- Main method demonstrates conversion from GTF to GFF format

**Usage/Command-Line**:
```bash
java GtfLine input.gtf [output.gff]
```
Converts GTF file to GFF format, writing to stdout.gff by default

**Dependencies**:
- `fileIO.ByteFile`
- `fileIO.ByteStreamWriter`
- `shared.LineParser1`
- `shared.Tools`
- `structures.ByteBuilder`
- `gff.GffLine`

---

## Format Conversion Utilities

### VcfToGff (VcfToGff.java)
**Purpose**: Converts Variant Call Format (VCF) files to Gene Feature Format (GFF) files.

**Core Function**: A command-line utility that parses input VCF files and transforms variant information into GFF format, preserving metadata and generating a standard GFF3 output with genomic feature annotations.

**Key Features**:
- Translates VCF variant records to GFF3 format
- Preserves VCF metadata by selectively copying header comments
- Supports input/output file specification via command-line arguments
- Generates a standard GFF3 header with version 3
- Implements robust file handling with overwrite and append options
- Uses byte-level file processing for memory efficiency

**Technical Implementation**:
- Uses `ByteFile` for input reading and `ByteStreamWriter` for output
- Converts `VCFLine` objects to `GffLine` objects for transformation
- Implements buffered writing with a `ByteBuilder` of 17000 bytes
- Filters VCF headers, selectively preserving non-standard comments
- Dynamic command-line argument parsing with flexible input handling

**Usage/Command-Line**:
```bash
java gff.VcfToGff in=input.vcf out=output.gff
```

**Dependencies**:
- `fileIO.ByteFile`: File reading
- `fileIO.ByteStreamWriter`: File writing
- `var2.VCFLine`: VCF line parsing
- `gff.GffLine`: GFF line generation
- `shared.Tools`: Utility methods
- `shared.Parser`: Argument parsing

---

## Specialized Processing Tools

### ParseCrispr (ParseCrispr.java)
**Purpose**: Command-line utility for parsing and extracting CRISPR annotation sequences from GFF files.

**Core Function**: Processes input GFF files to extract direct repeat CRISPR sequences, converting them into a FASTA-like format with sequential identifiers. Supports filtering CRISPR annotations, handling file compression, and separating valid and invalid entries.

**Key Features**:
- Extracts CRISPR direct repeat sequences from GFF files
- Converts sequences to uppercase FASTA format with sequential numbering
- Supports optional maximum line processing limit
- Handles compressed input files via ByteFile and ReadWrite utilities
- Generates separate output streams for valid and invalid annotations
- Provides runtime statistics including processed lines and bytes

**Technical Implementation**:
- Uses ByteFile for efficient file reading
- Implements line-by-line processing with strict CRISPR annotation validation
- Validates entries based on specific GFF annotation criteria
- Extracts sequence data using precise string indexing and character checking
- Utilizes ByteBuilder for efficient string manipulation

**Usage/Command-Line**:
```bash
java gff.ParseCrispr in=input.gff out=output.fasta invalid=invalid_entries.txt lines=1000
```
Supported parameters:
- `in`: Input GFF file
- `out`: Output FASTA file
- `invalid`: Optional file for invalid entries
- `lines`: Optional maximum number of lines to process
- `verbose`: Enable detailed logging

**Dependencies**:
- fileIO: ByteFile, ByteStreamWriter, FileFormat, ReadWrite
- shared: Parse, Parser, PreParser, Tools
- structures: ByteBuilder
- java.io: PrintStream

---

## Package Usage

The gff package provides essential functionality for:

**File Format Processing**:
- Converting between GFF3, GTF, VCF, and GenBank formats
- Parsing and validating genomic feature annotations
- Extracting specific feature types and attributes

**Quality Assessment**:
- Comparing annotation accuracy between reference and query files
- Computing precision, recall, and signal-to-noise ratios
- Evaluating gene prediction performance

**Feature Extraction**:
- Extracting genomic regions based on annotation criteria
- Filtering features by type, length, and attributes
- Processing strand-specific information

**Specialized Applications**:
- CRISPR sequence extraction and formatting
- GenBank to GFF conversion for genomic databases
- Legacy format support for older annotation systems

## Performance Considerations

- Most utilities support compressed input files (.gz, .bz2)
- Multi-threaded processing available for large-scale operations
- Memory-efficient byte-level parsing for large files
- Configurable processing limits for testing and debugging

## Dependencies

The gff package integrates with several BBTools components:
- `fileIO` package for efficient file handling
- `shared` package for parsing and utility functions
- `stream` package for sequence processing
- `var2` package for variant handling
- `prok` package for prokaryotic-specific features

---
*Documentation generated using evidence-based analysis of source code*