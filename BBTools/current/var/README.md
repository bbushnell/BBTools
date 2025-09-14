# Var Package

This package provides comprehensive genomic variation analysis tools for variant calling, filtering, and processing. It includes utilities for detecting SNPs, insertions, deletions, and complex variants from sequencing data.

## ApplyVarsToReference (ApplyVarsToReference.java)
**Purpose**: Apply genetic variations to a reference genome
**Core Function**: Process and modify chromosome sequences by applying a set of genetic variations
**Key Features**:
- Supports multiple variation types: SNP, insertion, deletion, replacement
- Can regenerate N-blocks in chromosome sequences
- Validates reference sequences during variation application
- Handles genome build transformations
**Usage**: Modifies reference genome by applying a list of genetic variations across specified chromosomes

## GenerateConsensusVariations (GenerateConsensusVariations.java)
**Purpose**: Processes and filters genomic variations across chromosomes
**Core Function**: Filters input variant files using coverage depth and consensus thresholds
**Key Features**:
- Filters variants based on minimum coverage depth
- Removes overlapping variations by selecting higher quality variant
- Supports multiple chromosome processing
- Tracks variation input and output statistics
**Usage**: Genomic variation filtering and quality control in sequence analysis

## GenerateVarlets (GenerateVarlets.java)
**Purpose**: Generates varlets (small variants) from aligned sequencing reads by identifying genomic variations
**Core Function**: Multi-threaded variant detection process that processes mapped reads to detect SNPs, insertions, deletions, and complex substitutions
**Key Features**:
- Multi-threaded variant processing with configurable thread count
- Support for single-end and paired-end read processing
- Configurable variant condensation and merging strategies
- Chromosome-specific output file generation
- Quality-based variant filtering
**Usage**: Primary variant calling tool in BBTools suite for detecting genomic variations from sequencing data

## GenerateVarlets2 (GenerateVarlets2.java)
**Purpose**: Splits output files across genomic blocks for low-memory variant detection with parallel processing
**Core Function**: Multi-threaded genomic variant caller that processes sequencing reads, extracts variants, and writes block-based output files
**Key Features**:
- Parallel processing of sequencing reads using configurable thread count
- Block-based memory management for large-scale variant detection
- Supports single-end and paired-end read processing
- Configurable variant condensation and merging options
**Usage**: Used in genomic variant calling pipelines to generate variant lists from sequencing data with minimal memory overhead

## GenerateVarlets3 (GenerateVarlets3.java)
**Purpose**: Splits output files across blocks for low memory usage with id-sorted site lists
**Core Function**: Generates genomic variants by processing read alignments and filtering variants using multiple quality control techniques
**Key Features**:
- Supports multi-threaded variant processing with configurable threads
- Performs site filtering and variant condensation
- Handles single-end and paired-end read processing
- Supports block-based memory partitioning for variant writing
**Usage**: Generates variant files from sequencing read alignment data with memory-efficient processing

## StackVariations (StackVariations.java)
**Purpose**: Main entry point for variant stacking and filtering across chromosomes
**Core Function**: Processes variant files by applying multi-level chromosome-specific quality filtering, merging, and stacking of genetic variants
**Key Features**:
- Supports parallel processing of chromosome variant files
- Performs variant merging and de-duplication
- Implements configurable filtering criteria for variant quality
- Tracks detailed variant statistics (SNPs, deletions, insertions)
**Usage**: Used in genomic variant analysis to filter and consolidate variant calls from sequencing data across multiple chromosomes

## StackVariations2 (StackVariations2.java)
**Purpose**: Processes genomic variants with configurable thread-based filtering and merging
**Core Function**: Filters and merges genomic variants across chromosomes using multi-threaded processing
**Key Features**:
- Thread-based chromosome variant processing
- Configurable quality filtering for genomic variants
- Merges variants with similar genomic positions
- Tracks detailed statistics about variant processing
**Usage**: Command-line tool for processing and filtering genomic variant files across entire chromosomes

## VarLine (VarLine.java)
**Purpose**: Represents a single genomic variant with detailed parsing and conversion capabilities
**Core Function**: Parses and manages variant data from tab-delimited files, supporting multiple version formats
**Key Features**:
- Converts variant data into standardized internal representation
- Supports heterozygous variant splitting
- Generates multiple string representations (toString, toSourceString)
- Handles complex variant type mapping and coordinate transformations
**Usage**: Used for parsing, representing, and manipulating genomic variant information in BBTools variant analysis pipeline

## Varlet (Varlet.java)
**Purpose**: Detailed representation of a genetic variation with comprehensive read and quality metadata
**Core Function**: Stores and processes variant information from sequencing reads, capturing strand, quality, and mapping details
**Key Features**:
- Packs quality metrics into single integer using bitwise operations
- Tracks variant location, read alignment, and strand information
- Supports parsing variant data from text files
- Implements complex comparison and scoring methods for variants
**Usage**: Used in variant calling and analysis to represent genetic variations with high-precision metadata

## Variation (Variation.java)
**Purpose**: Represents genomic variations with precise coordinate and sequence tracking
**Core Function**: Tracks genomic variations by chromosome, location, type, reference sequence, and variant sequence
**Key Features**:
- Supports multiple variation types: SNP, INS, DEL, DELINS
- Precise location tracking with 0-based coordinates
- Advanced intersection and overlap detection methods
- Comprehensive variant type classification
**Usage**: Used for tracking and analyzing genomic variations across chromosomes in genomic analysis pipelines