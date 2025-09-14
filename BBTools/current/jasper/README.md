# jasper Package - K-mer Position Analysis Tools
*Specialized utilities for analyzing k-mer positional distribution in sequencing data*

## Package Overview
The `jasper` package provides tools for analyzing k-mer occurrence patterns and positional distribution within DNA sequencing reads.

---

## Core Components

### KmerPosition (KmerPosition.java)
**Purpose**: Analyzes positional k-mer distribution in sequencing reads
**Core Function**: Processes DNA sequencing reads, identifying k-mer occurrences at specific positions using HashSet storage for efficient lookup
**Key Features**: 
- Supports paired-end sequencing reads
- Calculates k-mer position percentages
- Generates detailed output report
**Usage**: Command-line tool for genomic k-mer positional analysis
**Dependencies**: java.util, fileIO, stream, structures

### KmerPosition3 (KmerPosition3.java)
**Purpose**: Identifies kmer sequence matches between reads and reference sequences
**Core Function**: Processes high-throughput read files to find and report kmer positions that match a reference sequence, tracking occurrence frequencies
**Key Features**: 
  - Converts nucleotide sequences to 2-bit notation for fast comparison
  - Supports paired-end read analysis
  - Handles non-degenerate base kmer construction
**Usage**: Command-line tool for genomic sequence analysis
**Dependencies**: fileIO, stream, structures packages

---

## Package Usage
The jasper package enables k-mer positional analysis for understanding sequence composition patterns and reference matching in genomic data.

---
*Documentation generated using evidence-based analysis of source code*