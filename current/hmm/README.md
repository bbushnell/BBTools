# hmm Package - Hidden Markov Model Processing
*Tools for analyzing and processing HMM search results*

## Package Overview
The `hmm` package provides utilities for processing Hidden Markov Model search results, enabling protein sequence analysis and hit tracking functionality.

---

## Core Components

### HMMSearchReport (HMMSearchReport.java)
**Purpose**: Processes and summarizes Hidden Markov Model (HMM) search results
**Core Function**: Loads HMM search output, parses input files, and generates a protein summary map
**Key Features**: 
- Parses command-line HMM search results (line 211)
- Builds protein summary from search lines (lines 184-199)
- Supports compressed input files
**Usage**: CLI tool for analyzing HMM search results
**Dependencies**: fileIO, shared packages

### ProteinSummary (ProteinSummary.java)
**Purpose**: Tracks protein sequence hit lengths across HMM search results
**Core Function**: Maintains a mapping of reference model names to their maximum hit lengths, allowing tracking of best protein sequence matches
**Key Features**:
  - Stores best hit length for each reference model name
  - Allows adding new hit lengths, updating only if longer
  - Preserves query sequence name
**Usage**: Create with query name, add HMMSearchLine instances to track best matches
**Dependencies**: java.util.HashMap

### HMMSearchLine (HMMSearchLine.java)
**Purpose**: Precise parser for multi-field delimited Hidden Markov Model (HMM) search result lines
**Core Function**: Extracts 23 distinct fields from space-delimited byte arrays using robust parsing with error assertions and type conversion
**Key Features**: 
- Parses complex multi-space delimited lines with strict field extraction
- Supports multiple numeric types (int, float, double)
- Generates compact text representation via `toText()` method
**Usage**: Internal parsing for HMM search result processing
**Dependencies**: shared.Parse, shared.Tools, structures.ByteBuilder

---

## Package Usage
The hmm package enables processing of HMM search output files, tracking best protein hits, and analyzing sequence search results in bioinformatics pipelines.

---
*Documentation generated using evidence-based analysis of source code*