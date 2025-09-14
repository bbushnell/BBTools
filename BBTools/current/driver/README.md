# Driver Package

BBTools driver package contains **62 command-line utilities** for bioinformatics data processing, file manipulation, and sequence analysis. These tools provide essential functionality for genomic workflows including file concatenation, filtering, format conversion, and specialized genomic sequence operations.

---

## Classes

#### BBVersion (BBVersion.java)
**Purpose**: Utility class for printing BBTools version information
**Core Function**: Outputs version strings for BBTools and BBMap
**Key Features**:
  - Prints Shared.BBTOOLS_VERSION_STRING by default
  - Optionally prints Shared.BBMAP_VERSION_NAME with arguments
  - Static main method for standalone version checking
  - Simple, single-purpose version reporting utility
**Usage**: Used to quickly retrieve and display version information for BBTools and BBMap packages

#### ClearRam (ClearRam.java)
**Purpose**: Memory allocation testing utility that determines maximum allocatable memory
**Core Function**: Repeatedly allocates long arrays to trigger OutOfMemoryError and report memory limits
**Key Features**:
  - Dynamically creates 1MB long arrays until memory exhaustion
  - Runs two test cycles with forced garbage collection
  - Tracks and reports total memory allocation in MB
  - Includes writeJunk method for generating test data files
**Usage**: Diagnostic tool for testing Java memory allocation and heap size limits

#### CollateSpikeIn (CollateSpikeIn.java)
**Purpose**: Parses BBMap log files to extract spike-in control statistics for sequencing jobs.
**Core Function**: Reads log files, extracts job IDs, mapping percentages, and accuracy rates from BBMap alignment logs.
**Key Features**:
  - Processes log files with complex line parsing
  - Extracts normalized job ID from file paths
  - Calculates mapping and accuracy percentages
  - Outputs tab-delimited results for spike-in analysis
  - Handles variable log file formats
**Usage**: Post-processing tool for analyzing BBMap alignment logs, typically used in bioinformatics pipeline quality control

#### CompareReferenceGenomes (CompareReferenceGenomes.java)
**Purpose**: Base-by-base comparison tool for genome reference files across chromosomes 1-25
**Core Function**: Reads two genome file patterns, compares ChromosomeArray files for each chromosome
**Key Features**:
  - Supports flexible filename patterns with # placeholder
  - Compares chromosomes from 1 to 25 automatically
  - Detects index range and base-level mismatches
  - Prints detailed mismatch information (chromosome, position, bases)
  - Returns boolean indicating file similarity
**Usage**: Validate genome reference files by comparing two sets of chromosome files

#### CompareSequences (CompareSequences.java)
**Purpose**: Compare two chromosome sequences and analyze their base-level differences
**Core Function**: Compares two ChromosomeArray inputs, tracking base-by-base similarity, case differences, and N-base transformations
**Key Features**:
  - Tracks total matching and different bases
  - Detects base-to-N and N-to-base transformations
  - Calculates percentage of different bases
  - Identifies case-sensitive changes (toUpper/toLower)
  - Handles sequences of different lengths
**Usage**: Command-line tool for comparative genomic sequence analysis, useful for detecting mutations, sequencing errors, or genome variations

#### ConcatenateFiles (ConcatenateFiles.java)
**Purpose**: Command-line utility for concatenating multiple files from a directory or using a pattern
**Core Function**: Reads input files and merges their contents into a single output file or stream
**Key Features**:
  - Supports directory-based file concatenation
  - Handles file input streams with buffered reading
  - Supports writing to stdout or specified output file
  - Sorts input files before concatenation
  - Provides processing progress via console output
**Usage**: Used for combining multiple files into a single file, typically for data consolidation or preprocessing tasks

#### ConcatenateTextFiles (ConcatenateTextFiles.java)
**Purpose**: Concatenate multiple text files or directories into a single output file.
**Core Function**: Reads multiple input files/directories, buffers their contents, and writes them sequentially to an output file using a multi-threaded approach.
**Key Features**:
- Supports single files and entire directory concatenation
- Multi-threaded writing with buffered input processing
- Configurable append and overwrite modes
- Handles large files with chunked buffer processing
**Usage**: Command-line utility for merging text files, useful for log consolidation, data preprocessing, and file aggregation.

#### Concatenator (Concatenator.java)
**Purpose**: Command-line utility for merging multiple text files into a single output
**Core Function**: Reads multiple input files and writes their contents sequentially to a target file or standard output
**Key Features**:
  - Accepts comma-separated list of input file names
  - Supports writing to file or printing to console
  - Provides static methods for both command-line and programmatic file merging
  - Handles file reading and writing with robust TextFile and TextStreamWriter
**Usage**: Useful for combining text files, log consolidation, and batch file processing

#### ConvertSamToAln (ConvertSamToAln.java)
**Purpose**: Converts SAM alignment files to compact three-column ALN format
**Core Function**: Batch converts SAM files by extracting essential mapping information and creating compressed ALN representations
**Key Features**:
  - Handles compressed input files (zip, gz, bz2)
  - Extracts chromosome, location, and strand information
  - Supports batch conversion of multiple SAM files
  - Writes compressed ALN output files
  - Filters out unmapped reads
  - Preserves mapping orientation (forward/reverse)
**Usage**: Command-line tool for transforming detailed SAM alignment files into lightweight, compressed text format for further processing

#### CorrelateIdentity (CorrelateIdentity.java)
**Purpose**: Correlate identity matrices by sampling and shuffling corresponding matrix elements
**Core Function**: Reads two input matrices, extracts and shuffles matching lower triangular elements, and writes correlated pairs to output
**Key Features**:
- Supports command-line input/output file specification
- Implements sampling with configurable sample rate
- Performs cross-matrix element correlation
- Supports file overwrite and append modes
- Provides flexible matrix element extraction
**Usage**: Used for comparing and sampling corresponding elements from two identity matrices in bioinformatics or data analysis workflows

#### CountRNAs (CountRNAs.java)
**Purpose**: Count and categorize genes across human chromosomes by type
**Core Function**: Iterates through chromosomes 1-24, classifying genes as pseudogenes, untranslated, or coding
**Key Features**:
  - Processes gene data for 24 chromosomes
  - Distinguishes between pseudogenes, coding, and untranslated genes
  - Uses predefined gene mapping configuration
  - Outputs detailed gene type statistics
**Usage**: Provides genetic classification summary for genome analysis, useful in bioinformatics research

#### CountSharedLines (CountSharedLines.java)
**Purpose**: Filters and counts shared lines between text files
**Core Function**: Compares text files to determine number of common lines
**Key Features**:
  - Supports case-sensitive and case-insensitive line matching
  - Allows substring comparison between lines
  - Configurable input and output file handling
  - Supports multiple input file sets
  - Optional line replacement and prefix mode
**Usage**: Utility for comparing line contents across multiple text files, useful for identifying shared text elements

#### EstherFilter (EstherFilter.java)
**Purpose**: BLAST-based sequence filtering tool that processes and filters genomic sequences by BLAST alignment score
**Core Function**: Runs BLAST command, parses tabular output, filters sequences based on alignment score cutoff
**Key Features**:
  - Supports both name-only and FASTA sequence output modes
  - Configurable score cutoff for sequence filtering
  - Binary search-based efficient sequence matching
  - Handles BLAST (-m 8) tabular output format
**Usage**: Command-line utility for genomic sequence filtering using BLAST results, useful in bioinformatics data preprocessing

#### FilterAssemblySummary (FilterAssemblySummary.java)
**Purpose**: Filters assembly summary files based on taxonomic criteria
**Core Function**: Processes text files, applying taxonomic filters to retain or remove lines
**Key Features**:
  - Supports command-line input for filter configuration
  - Processes files with customizable read limits
  - Filters lines based on taxonomic node numbers
  - Tracks processing statistics like lines processed and retained
  - Handles compressed and uncompressed input files
**Usage**: Used in bioinformatics workflows to filter assembly summary files by taxonomic classification

#### FilterLines (FilterLines.java)
**Purpose**: Filters text lines by exact match or substring using flexible matching criteria.
**Core Function**: Processes input text files, filtering lines based on configurable name matching rules.
**Key Features**:
  - Supports exact name matching
  - Substring matching (line to name, name to line)
  - Case-sensitive and case-insensitive modes
  - Optional prefix matching
  - Text replacement during filtering
  - Maximum line processing limit
**Usage**: Command-line text file filtering tool for selective line extraction or removal

#### FilterReadsByName (FilterReadsByName.java)
**Purpose**: Filters sequencing reads based on read names using flexible matching criteria.
**Core Function**: Processes input reads and selectively retains or excludes reads matching specified name patterns.
**Key Features**:
  - Supports exact name matching with case sensitivity controls
  - Allows substring matching between read names and headers
  - Configurable inclusion/exclusion filtering
  - Supports prefix-based matching
  - Handles paired-end read filtering
  - Configurable minimum read length filtering
**Usage**: Used in bioinformatics pipelines to filter sequencing datasets based on read identifiers, enabling targeted read selection or removal

#### FindMotifs (FindMotifs.java)
**Purpose**: Identifies and analyzes gene start motifs across chromosomes
**Core Function**: Detects gene start locations using probabilistic motif matching across different chromosomes
**Key Features**:
  - Supports analysis of multiple chromosomes (1-22)
  - Uses probabilistic motif models for gene start detection
  - Provides detailed histogram of gene start occurrences
  - Supports strand-specific gene start analysis
  - Implements multiple analysis strategies (in-frame, stronger matches)
**Usage**: Computational genomics tool for identifying potential gene start sites using probabilistic sequence matching

#### FixChr (FixChr.java)
**Purpose**: Converts genomic coordinate reference files from grch38 to hg19 format
**Core Function**: Adds "chr" prefix to chromosome names in SAM file headers and contig definitions
**Key Features**:
  - Reads input SAM file line by line
  - Adds "chr" prefix to non-header lines
  - Modifies contig ID lines to include "chr" prefix
  - Writes modified lines to output file
  - Preserves original file structure
**Usage**: Command-line utility for genomic file coordinate system conversion

#### FixDumbFile (FixDumbFile.java)
**Purpose**: Transforms and consolidates tabular data files by reorganizing entries based on unique keys
**Core Function**: Processes input files, groups entries by first column, and generates a restructured output file
**Key Features**:
  - Uses LinkedHashMap to group entries by first column
  - Skips header lines starting with "library_name"
  - Preserves order of first encountered entries
  - Generates a new file with reorganized data columns
  - Supports tab-separated input files
**Usage**: Command-line utility for reformatting and consolidating tabular data files with specific structural requirements

#### GenerateNoCallsFromCoverage (GenerateNoCallsFromCoverage.java)
**Purpose**: Generates no-call variation lines based on genome coverage levels
**Core Function**: Creates VarLine objects for genomic regions with insufficient coverage
**Key Features**:
  - Handles different ploidy scenarios (haploid/diploid chromosomes)
  - Tracks coverage levels across chromosome arrays
  - Generates no-call variations for low-coverage genomic regions
  - Supports multiple chromosome types (autosomes, sex chromosomes)
  - Removes duplicate no-call variations
**Usage**: Part of variant calling pipeline to mark low-coverage genomic regions as no-call

#### GetSequence (GetSequence.java)
**Purpose**: Command-line utility for extracting and manipulating genomic sequences from chromosomes
**Core Function**: Retrieves genomic sequence segments by chromosome, location, and optional strand
**Key Features**:
  - Supports flexible chromosome and location input parsing
  - Handles genomic build selection via command-line arguments
  - Outputs original and reverse complement sequences
  - Converts retrieved sequences to amino acid representations
  - Supports base-0 and base-1 coordinate systems
**Usage**: Genomic sequence extraction and analysis in bioinformatics research and sequence processing workflows

#### GetUniquePrefixes (GetUniquePrefixes.java)
**Purpose**: Extracts unique sequence prefixes from a FASTA-formatted file
**Core Function**: Reads a text file, generates unique sequence prefixes of specified length
**Key Features**:
  - Trims sequences to specified prefix length
  - Uses HashSet to ensure uniqueness of prefixes
  - Preserves original sequence headers
  - Handles multi-line sequence entries
**Usage**: Filters sequence files to get unique initial subsequences, useful for reducing redundancy in genomic data processing

#### Grep (Grep.java)
**Purpose**: Simple command-line text file line matching utility
**Core Function**: Searches text files for lines containing a specified substring
**Key Features**:
  - Reads input file from first command-line argument
  - Searches for substring specified in second argument
  - Prints matching lines to standard output
  - Uses TextFile for efficient file reading
  - Automatically closes file after processing
**Usage**: Command-line tool for filtering text files by substring, similar to Unix grep command

#### Life (Life.java)
**Purpose**: Command-line Game of Life simulation with configurable grid parameters
**Core Function**: Simulates cellular automaton using randomized initial state and Conway's Game of Life rules
**Key Features**:
  - Configurable grid dimensions (rows/columns)
  - Random initial cell state generation
  - Configurable maximum simulation cycles
  - Optional display mode with console visualization
  - Supports delay between simulation cycles
  - Implements periodic boundary conditions
**Usage**: Scientific visualization and simulation of cellular automaton evolution

#### LineCount (LineCount.java)
**Purpose**: Command-line utility for counting non-blank lines in text files
**Core Function**: Reads a text file and returns the count of non-blank lines
**Key Features**:
  - Uses TextFile class for efficient file line counting
  - Handles single file input from command line arguments
  - Prints total non-blank line count to console
  - Simple, lightweight line counting mechanism
**Usage**: Quickly determine the number of non-blank lines in a text file

#### LoadReads (LoadReads.java)
**Purpose**: Analyzes memory usage during sequence data loading operations
**Core Function**: Profiles memory consumption, reads, and I/O performance for sequence files
**Key Features**:
  - Calculates memory used during read processing
  - Estimates file memory requirements
  - Tracks reads, bases, and quality score processing
  - Supports single and paired-end input files
  - Generates detailed memory and performance statistics
**Usage**: Diagnostic tool for understanding memory overhead in sequence data processing workflows

#### LookAtID (LookAtID.java)
**Purpose**: Analyzes text files containing SiteScoreR data to inspect numeric IDs and detect integer overflow
**Core Function**: Reads input text files, parses SiteScoreR entries, and checks for numeric ID integer overflow
**Key Features**:
  - Detects numeric IDs exceeding Integer.MAX_VALUE
  - Tracks maximum observed numeric ID
  - Prints detailed overflow information including raw data
  - Uses SiteScoreR parsing for structured data analysis
  - Handles text file processing with TextFile utility
**Usage**: Diagnostic tool for detecting potential integer overflow in large numeric ID datasets

#### MakeTestScript (MakeTestScript.java)
**Purpose**: Generates shell scripts for benchmarking sequence alignment tools
**Core Function**: Creates standardized test scripts for various genomic mapping algorithms like BWA, Bowtie, GSNAP, and BBMap
**Key Features**:
  - Supports multiple alignment tools (BWA, Bowtie2, GSNAP, BBMap, SNAP, etc.)
  - Configurable read count and read length parameters
  - Generates time tracking and SAM file grading scripts
  - Supports multiple benchmarking modes and configurations
**Usage**: Generates reproducible alignment tool performance test scripts for bioinformatics research

#### MakeTestScriptScoreOnly (MakeTestScriptScoreOnly.java)
**Purpose**: Generates shell scripts for scoring and evaluating alignment test results
**Core Function**: Dynamically creates test scripts with configurable read counts and alignment score parameters
**Key Features**:
  - Generates scripts for multiple alignment tools (SSAHA2, BWA, BBMap, Bowtie, BFAST)

#### MeasureGene (MeasureGene.java)
**Purpose**: Analyzes gene structures and characteristics across chromosomes
**Core Function**: Calculates frequency and matching probabilities of gene and exon start/stop sequences
**Key Features**:
  - Analyzes gene structures for specific chromosomes
  - Measures exon frequency and match strengths
  - Supports gene filtering for normal genes
  - Computes average motif match probabilities
**Usage**: Used for detailed genomic sequence analysis, particularly for chromosome-specific gene investigations

#### MergeBigelow (MergeBigelow.java)
**Purpose**: Merge and process Bigelow Laboratory data files with advanced text manipulation
**Core Function**: Merge two input text files by matching first column entries and combining data rows
**Key Features**:
  - Supports delimiter-based file parsing
  - Performs case-sensitive and case-insensitive data matching
  - Handles multiple input and output file formats
  - Configurable read limit and verbosity settings
  - Error handling for file and input validation
**Usage**: Used for merging scientific or research data files with complex text processing requirements

#### MergeCoverageOTU (MergeCoverageOTU.java)
**Purpose**: Merge and aggregate Operational Taxonomic Unit (OTU) coverage statistics from input files
**Core Function**: Reads input coverage statistics, aggregates multiple entries for the same OTU, and writes merged results to output file
**Key Features**:
- Supports flexible input file parsing with command-line input/output specification
- Aggregates coverage statistics for duplicate OTU entries
- Preserves header information from input files
- Uses LinkedHashMap for efficient OTU tracking and merging
- Handles multiple coverage statistic files
**Usage**: Consolidate OTU coverage data from multiple sources into a single comprehensive output file, useful in microbial ecology and genomic research

#### MergeTextFiles (MergeTextFiles.java)
**Purpose**: Merge two tab-delimited text files by matching key columns
**Core Function**: Combines two text files using a specified key column, merging rows with matching keys
**Key Features**:
  - Supports merging files with different column structures
  - Preserves original header from first input file
  - Sorts merged output alphabetically by key
  - Handles missing values in either input file
  - Creates a tab-delimited output combining both files
**Usage**: Typically used for combining related data files with a common identifier column

#### MergeTextFiles2 (MergeTextFiles2.java)
**Purpose**: Merges two text files by matching specific columns and generating a combined tabular output
**Core Function**: Combines text files by matching keys in specified columns, handling missing entries and maintaining original order
**Key Features**:
  - Supports merging files with variable column widths
  - Preserves header from first input file
  - Sorts output by merged keys
  - Handles missing entries by padding with null values
  - Uses hashtables for efficient key matching
**Usage**: Merging tabular data files with a common key column, such as joining datasets from different sources

#### MoveFiles (MoveFiles.java)
**Purpose**: Organize files into chromosome-specific directories
**Core Function**: Moves files matching specific chromosome patterns to dedicated subdirectories
**Key Features**:
  - Processes files in a root directory
  - Creates subdirectories for chromosomes 1-22
  - Filters files based on chromosome-specific naming patterns
  - Copies files to appropriate chromosome subdirectories
  - Handles file path and name normalization
**Usage**: Used for organizing genomic or chromosomal data files into structured directory layout

#### ParseCrossblockResults (ParseCrossblockResults.java)
**Purpose**: Parse and analyze crossblock results to track contigs and bases kept or discarded
**Core Function**: Process text file with crossblock filtering results, tracking length and removal status of contigs
**Key Features**:
  - Reads tabbed input files with crossblock filtering information
  - Tracks total bases and contigs kept and discarded
  - Supports verbose logging and error state tracking
  - Flexible input/output file handling with overwrite and append options
**Usage**: Used in bioinformatics workflows to summarize filtering results from crossblock operations, providing statistics on sequence retention

#### PlotGC (PlotGC.java)
**Purpose**: Calculates and outputs GC content distribution for genomic sequences
**Core Function**: Processes input sequence files, dividing reads into intervals and computing GC content percentages
**Key Features**:
  - Supports interval-based GC content calculation
  - Processes FASTA/FASTQ input files with concurrent reading
  - Configurable interval size and offset
  - Optional short bin handling
  - Outputs detailed GC content statistics
**Usage**: Bioinformatics tool for analyzing sequence composition and GC distribution across genomic reads

#### PrintEnv (PrintEnv.java)
**Purpose**: Utility for displaying system environment information and current machine details
**Core Function**: Retrieves and prints sorted system environment variables and local machine hostname
**Key Features**:
  - Prints current system timestamp
  - Retrieves and alphabetically sorts all environment variables
  - Displays each environment variable and its value
  - Attempts to retrieve and print local machine hostname
**Usage**: Diagnostics and system configuration verification tool, useful for understanding runtime environment

#### ProcessFragMerging (ProcessFragMerging.java)
**Purpose**: Processes and formats performance data for BBMerge comparisons
**Core Function**: Parses input log files and extracts key performance metrics like execution time, read usage, mapping rates, and error statistics
**Key Features**:
  - Converts time formats from minutes:seconds to total seconds
  - Extracts performance metrics from log files
  - Prints formatted performance data using tab-separated values
  - Handles multiple log entry types (real time, reads used, mapping, error rates)
**Usage**: Command-line utility for parsing and summarizing BBMerge performance benchmark logs

#### ProcessSpeed (ProcessSpeed.java)
**Purpose**: Parses and processes performance timing logs for BBMerge benchmarks
**Core Function**: Extracts and formats time metrics from command line timing outputs
**Key Features**:
  - Converts time strings from minutes:seconds to decimal seconds
  - Parses real, user, and system time measurements
  - Prints formatted performance metrics
  - Handles log files with specific timing and statistic entries
**Usage**: Performance log processing utility for BBMerge benchmark comparisons

#### ProcessSpeed2 (ProcessSpeed2.java)
**Purpose**: Parse and convert time measurements from command output
**Core Function**: Converts time log lines with minutes and seconds into decimal seconds
**Key Features**:
  - Reads time log input from text files
  - Converts time formats (m:ss) to total seconds
  - Prints formatted time measurements
  - Handles real, user, and system time types
**Usage**: Used for processing and standardizing time measurement logs from command executions

#### ProcessWebcheck (ProcessWebcheck.java)
**Purpose**: Web-based quality control data processing utility for parsing and analyzing web request metrics
**Core Function**: Processes input files containing web request data, extracts latency and response code statistics
**Key Features**:
  - Parses web request lines with timestamp, URL, response code, and latency
  - Generates extended statistics on processed lines
  - Separates valid and invalid request entries
  - Tracks pass/fail latency and response codes
  - Supports configurable input/output file handling
**Usage**: Used for analyzing web request performance, generating metrics on response times and error rates

#### ReduceSilva (ReduceSilva.java)
**Purpose**: Command-line tool for filtering FASTQ reads based on taxonomic information
**Core Function**: Filters reads by unique taxonomic identifiers in a specified column of read headers
**Key Features**:
  - Processes reads from FASTQ input files
  - Allows configurable column for taxonomic identifier extraction
  - Tracks and filters out duplicate taxa
  - Supports both paired and unpaired read formats
  - Provides detailed processing statistics
**Usage**: Used in bioinformatics workflows to reduce redundancy in Silva database reads by filtering out duplicate taxa

#### RenameAndMux (RenameAndMux.java)
**Purpose**: Rename and multiplex FASTQ read files with flexible input/output handling
**Core Function**: Processes FASTQ read files, renaming reads and merging multiple input files into single or paired output streams
**Key Features**:
  - Supports single and paired-end read file processing
  - Dynamically handles input file naming with '#' replacement
  - Multithreaded file processing
  - Configurable interleaved output settings
  - Atomic read and base counting
**Usage**: Command-line utility for renaming and consolidating sequencing read files during bioinformatics data preprocessing

#### RenameByHeader (RenameByHeader.java)
**Purpose**: Renames files based on their header information, specifically for biological sequence files
**Core Function**: Extracts metadata from file headers to generate new filenames with taxonomic or sample information
**Key Features**:
  - Supports processing individual files or entire directories
  - Extracts species and sample names from file headers
  - Handles FASTQ/FASTA files with specific header formats
  - Preserves original file extension in new filename
  - Supports verbose logging for troubleshooting
**Usage**: Used in bioinformatics workflows to standardize file naming conventions based on sequence metadata

#### RenameNcbiToTid (RenameNcbiToTid.java)
**Purpose**: Converts NCBI sequence headers to TaxID-based headers in text files
**Core Function**: Transforms FASTA headers from ">ncbi..." format to ">tid..." format, adding a taxonomic identifier
**Key Features**:
  - Replaces ">ncbi" prefix with ">tid"
  - Inserts a "|" separator after the TID
  - Processes text files line by line
  - Supports input/output file specification
  - Handles compressed input/output files
  - Configurable verbosity and error handling
**Usage**: Standardizes sequence headers in bioinformatics sequence files, particularly for taxonomic identification and parsing

#### RenameRefseqFiles (RenameRefseqFiles.java)
**Purpose**: Utility for standardizing RefSeq genome file naming conventions
**Core Function**: Renames RefSeq genome files by prefixing with "refseq_" and taxonomy ID
**Key Features**:
  - Loads taxonomy tree from default taxonomy file
  - Iterates through all taxonomy nodes
  - Checks for existing genome files
  - Renames files with "refseq_" prefix
  - Uses taxonomy ID in filename
**Usage**: Standardizes genome file naming for bioinformatics research and data organization

#### Sample (Sample.java)
**Purpose**: A simple file processing utility for reading and writing text files.
**Core Function**: Reads input file line by line and writes contents to output file with optional parsing capability.
**Key Features**:
  - Line-by-line file reading using BufferedReader
  - Line-by-line file writing using PrintWriter
  - Separate methods for input and output stream creation
  - Basic error handling with System.exit() on file access failures
**Usage**: Command-line utility for basic file copying or transformation, with placeholder for custom parsing logic

#### Search (Search.java)
**Purpose**: Provides utility methods for genomic feature searching and range point location.
**Core Function**: Implements binary and linear search algorithms for finding genes and ranges that intersect with specific genomic points.
**Key Features**:
  - findGenes(): Searches for genes intersecting a specific point
  - findGenesBinary(): Performs binary search to find nearby genes
  - findPointBinary(): Efficient binary search for locating ranges
  - findPointLinear(): Linear search method for range point location
  - containsPointBinary(): Checks if a point is within a range threshold
**Usage**: Used for efficient genomic feature localization in genomic data processing workflows

#### SelectReads (SelectReads.java)
**Purpose**: Filters SAM/BAM files to select reads with specific CIGAR string characteristics
**Core Function**: Processes input SAM/BAM files, selectively writing reads matching specified deletion or alignment criteria
**Key Features**:
  - Supports filtering by CIGAR string symbols (M, S, D, I, C)
  - Configurable minimum length for selected reads
  - Limits total number of reads processed
  - Handles compressed input files transparently
  - Preserves header information from input file
**Usage**: Command-line tool for extracting reads with specific alignment properties from sequencing data files

#### SniffSplices (SniffSplices.java)
**Purpose**: Detects and analyzes splice sites in genomic sequences
**Core Function**: Calculates match strengths and percentiles for different genomic motifs at sequence positions
**Key Features**:
  - Supports multiple motif types (exon starts/stops, gene starts/stops)
  - Calculates motif match strength and normalized percentile
  - Handles both forward and reverse complement sequences
  - Flexible command-line motif selection
**Usage**: Genomic sequence analysis for identifying potential splice sites and regulatory regions

#### SummarizeContamReport (SummarizeContamReport.java)
**Purpose**: Summarize contamination report files by processing taxonomic classification data
**Core Function**: Reads contamination reports, aggregates sequence unit and read counts for different taxonomic entities
**Key Features**:
  - Processes multiple input contamination report files
  - Aggregates sequence counts by taxonomic name
  - Supports custom minimum read and sequence unit filtering
  - Integrates with TaxTree for detailed taxonomic information
  - Outputs comprehensive tab-delimited summary
**Usage**: Bioinformatics tool for consolidating and analyzing contamination report data from sequencing experiments

#### SummarizeCoverage (SummarizeCoverage.java)
**Purpose**: Analyzes and summarizes coverage statistics from input files
**Core Function**: Processes multiple input files to identify primary and secondary coverage information
**Key Features**:
  - Parses command-line arguments for input files
  - Identifies primary and secondary coverage metrics
  - Calculates count and megabyte statistics for each file
  - Outputs tab-separated summary with file details
  - Supports multiple input file formats
  - Handles file existence and path resolution
**Usage**: Used to aggregate and summarize coverage information from bioinformatics or file analysis workflows

#### SummarizeCrossblock (SummarizeCrossblock.java)
**Purpose**: Summarizes processing statistics for multiple crossblock-processed input files
**Core Function**: Iterates through input files, runs ParseCrossblockResults for each, and generates summary report of file processing metrics
**Key Features**:
  - Processes multiple input files in sequence
  - Generates tabular output with per-file metrics
  - Tracks contigs and bases kept/discarded
  - Handles error cases for individual file processing
  - Supports input file lists via comma-separated arguments
**Usage**: Utility for aggregating and reporting statistics from crossblock processing across multiple input files

#### SummarizeMSDIN (SummarizeMSDIN.java)
**Purpose**: Summarizes match/substitution/insertion/deletion/N rates for BBMap alignment runs
**Core Function**: Parses BBMap log files to extract and aggregate alignment quality metrics
**Key Features**:
  - Extracts match, substitution, deletion, insertion, and N rates
  - Supports configurable rate tracking via boolean flags
  - Uses TextFile for file reading
  - Aggregates metrics across multiple BBMap runs
  - Outputs tab-separated summary statistics
**Usage**: Command-line utility to analyze BBMap alignment performance metrics, typically used for batch processing alignment logs

#### SummarizeQuast (SummarizeQuast.java)
**Purpose**: Processes and summarizes QUAST (Quality Assessment Tool for Genome Assemblies) output files
**Core Function**: Aggregates and optionally normalizes metrics from multiple QUAST assembly evaluation reports
**Key Features**:
  - Reads QUAST TSV files with assembly metrics
  - Supports optional metric normalization
  - Can output summary statistics in boxplot or raw formats
  - Configurable metric filtering via command-line parameters
  - Handles multiple input files simultaneously
**Usage**: Command-line tool for consolidating and analyzing genome assembly quality metrics across multiple assemblies or runs

#### SummarizeSealCrosstalk (SummarizeSealCrosstalk.java)
**Purpose**: Calculates and reports contamination statistics across sequencing seal mapping results.
**Core Function**: Process multiple sequencing stat files to quantify correct, contaminated, and ambiguous read mappings.
**Key Features**:
  - Reads multiple SealStats files from different sequencing runs
  - Calculates reads per million (PPM) for various mapping categories
  - Supports primary mapping or custom mapping strategies
  - Generates sorted contamination statistics output
  - Tracks total, matched, primary, and contaminated reads
**Usage**: Used in bioinformatics to analyze cross-sample contamination in sequencing experiments, providing quantitative insights into mapping accuracy and potential sample interference.

#### SummarizeSealStats (SummarizeSealStats.java)
**Purpose**: Command-line tool for summarizing statistical data from Seal files
**Core Function**: Processes input files to summarize primary and secondary sequence counts and base statistics
**Key Features**:
  - Supports multiple input file processing
  - Calculates total and per-file sequence metrics
  - Optional filtering for taxa, barcodes, and location
  - Customizable output denominator calculation
**Usage**: Generates summary statistics for bioinformatics sequencing data analysis

#### TestCompressionSpeed (TestCompressionSpeed.java)
**Purpose**: Benchmark compression performance across different ZIP compression levels
**Core Function**: Systematically tests file compression and decompression efficiency by measuring time and output size for compression levels 0-9
**Key Features**:
  - Measures compression time for each ZIP level
  - Calculates resulting compressed file size
  - Supports file compression and decompression timing
  - Uses ReadWrite utility for flexible output stream handling
  - Prints detailed performance metrics for each compression level
**Usage**: Performance benchmarking tool for evaluating text file compression efficiency at various compression settings

#### TestLockSpeed (TestLockSpeed.java)
**Purpose**: Benchmark and compare performance of different concurrent counting mechanisms in Java
**Core Function**: Measure execution speed for various thread synchronization techniques like locks, atomic operations, and volatile fields
**Key Features**:
  - Supports multiple thread synchronization modes (unlocked, locked, atomic, volatile)
  - Configurable thread count and maximum iteration count
  - Precise timer to measure operation speed
  - Calculates operations per second (giga operations)
**Usage**: Performance testing and comparing thread synchronization overhead in multi-threaded Java applications

#### Translator (Translator.java)
**Purpose**: Transforms genomic variations and variant lines between different genome builds
**Core Function**: Translates genomic location references using chain line mapping for coordinate conversion
**Key Features**:
  - Supports translation of variant lines and variations between genome builds
  - Handles coordinate mapping using binary search on chain lines
  - Processes single and multi-chromosome variant sets
  - Supports both plus and minus strand coordinate conversions
  - Handles point and range-based variant translations
  - Performs base complement conversions for reverse strand mappings
**Usage**: Used in genomic coordinate system conversion, enabling cross-build variant analysis and mapping

#### Translator2 (Translator2.java)
**Purpose**: Genome coordinate translation utility for converting genomic locations between different genome builds
**Core Function**: Translates genomic coordinates between genome builds 36 and 37, handling chromosome and strand mapping
**Key Features**:
  - Supports translation between genome builds 36 and 37
  - Uses ChainLine mapping for precise coordinate conversion
  - Handles chromosome-specific coordinate translation
  - Supports both command-line and programmatic coordinate mapping
  - Manages strand orientation during translation
**Usage**: Used in bioinformatics research for converting genomic coordinates between different genome assembly versions

#### TransposeTextFile (TransposeTextFile.java)
**Purpose**: Transposes text files by converting rows to columns and vice versa
**Core Function**: Reads a text file, restructures its contents by swapping rows and columns
**Key Features**:
  - Supports optional line skipping during transposition
  - Handles chromosome-based filename substitution
  - Uses StringBuilder for efficient string manipulation
  - Writes transposed output to a new file with ".transposed" suffix
**Usage**: Command-line tool for converting tabular data between row and column orientations, useful in bioinformatics and data processing workflows

#### TrimSamFile (TrimSamFile.java)
**Purpose**: Filters and trims SAM files by removing reads mapped to a specific genomic region
**Core Function**: Identifies and removes SAM file entries based on read mapping positions within a specified genomic interval
**Key Features**:
  - Filters reads mapped to specific scaffold/chromosome
  - Supports genomic region trimming via command-line coordinates
  - Handles paired-end read filtering
  - Preserves original SAM file headers
  - Removes reads mapped across specified genomic boundaries
**Usage**: Command-line utility for selective SAM file read filtering, useful for extracting reads from specific genomic regions

---

*This documentation covers all 62 command-line utilities in the BBTools driver package, providing comprehensive coverage of bioinformatics tools for file processing, sequence analysis, and genomic data manipulation.*