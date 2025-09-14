# covid Package - COVID-19 Genomic Coverage Analysis
*Specialized utilities for processing and analyzing COVID-19 sequencing coverage data*

## Package Overview
The `covid` package provides tools specifically designed for analyzing genomic coverage data from COVID-19 sequencing projects, enabling researchers to assess sequencing depth and genome coverage quality.

---

## Core Classes

#### SummarizeCoverage (SummarizeCoverage.java)
**Purpose**: A utility for processing and summarizing genomic coverage data, specifically designed for analyzing COVID-19 related sequencing coverage.

- **Core Function**: Reads coverage files and generates statistical summaries of genome coverage across different depth thresholds.

- **Key Features**:
  - Processes multiple input coverage files in text format
  - Calculates average coverage and percentage of genome covered at various depth levels
  - Supports compressed input files (via .gz or .bz2 extensions)
  - Customizable processing with command-line parameters
  - Outputs tab-separated summary statistics

- **Usage**: Analyze genomic sequencing depth for COVID-19 samples, helping researchers understand genome coverage.

- **Technical Implementation**:
  - Generates coverage histogram with depth levels up to 20x
  - Calculates percentage of genome covered at thresholds: 1x, 2x, 3x, 4x, 5x, 10x, and 20x
  - Handles multiple input files with flexible file naming
  - Uses ByteFile and ByteStreamWriter for efficient file processing

- **Input Processing**:
  ```
  Input Format: Tab-separated coverage file
  Example Header: Chromosome   Position   Coverage
  ```

- **Output Format**:
  ```
  #Sample  AvgCov  %>=1x  %>=2x  %>=3x  %>=4x  %>=5x  %>=10x  %>=20x
  SampleName  15.23  99.5   98.7   97.2   95.6   94.1   85.3    70.2
  ```

- **Key Parameters**:
  - `lines=N`: Limit processing to N lines
  - `refbases=X`: Set reference genome length
  - `verbose`: Enable detailed logging

- **Error Handling**:
  - Validates input and output file accessibility
  - Throws runtime exceptions for configuration errors
  - Tracks and reports processing statistics

- **Performance Considerations**:
  - Multi-threaded file reading
  - Efficient byte-level processing
  - Minimal memory footprint

- **Limitations**:
  - Requires pre-generated coverage files
  - Assumes tab-separated input format
  - Designed for text-based coverage reports

---

## Package Usage
The covid package is specifically designed for COVID-19 genomic analysis workflows, particularly useful for:
- Quality assessment of COVID-19 sequencing data
- Coverage analysis for viral genome sequencing
- Statistical summarization of sequencing depth
- Batch processing of multiple COVID-19 samples
- Research applications requiring coverage metrics

## Dependencies
- Relies on BBTools core file I/O utilities (ByteFile, ByteStreamWriter)
- Integrates with standard BBTools data processing frameworks
- Supports standard compressed file formats

---
*Documentation generated using evidence-based analysis of source code*