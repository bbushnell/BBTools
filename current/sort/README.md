# BBTools sort Package

Read sorting, shuffling, and comparison utilities providing comprehensive infrastructure for organizing sequencing data by various criteria including genomic position, quality, taxonomy, and custom comparators.

---

## ContigLengthComparator (ContigLengthComparator.java)
**Purpose**: Sorts contigs based on multiple criteria with customizable order
**Core Function**: Compares Contig objects using length, coverage, base sequence, and ID
**Key Features**:
- Sorts primarily by contig length
- Secondary sort by coverage
- Tertiary sort by base sequence
- Final tiebreaker by contig ID
**Usage**: Used for ordering and sorting Contig objects in assembly and analysis workflows

## MergeSorted (MergeSorted.java)
**Purpose**: Sorts reads by name, potentially from multiple input files.
**Core Function**: Recursively merges and sorts read files using configurable comparators.
**Key Features**:
- Supports multiple sorting methods (length, name, quality, position)
- Can merge files from multiple input sources
- Handles file compression and temporary file management
- Configurable sorting order (ascending/descending)
**Usage**: Used for preprocessing sequencing reads, enabling efficient sorting and merging of large read datasets

## ReadComparator (ReadComparator.java)
**Purpose**: Abstract base class for implementing comparators to sort Read objects
**Core Function**: Provides a template method pattern for configurable Read sorting
**Key Features**:
- Abstract method to set sorting order (ascending/descending)
- Implements Java's Comparator interface for Read objects
- Allows custom sorting implementations for Read streams
**Usage**: Used as a base class for creating specialized Read sorting strategies in BBTools stream processing

## ReadComparatorClump (ReadComparatorClump.java)
**Purpose**: Sorts reads similarly to Clumpify algorithm using k-mer based comparison
**Core Function**: Compares reads based on their numeric ID, strand, and start position with k-mer hashing
**Key Features**:
- Compares both read and mate reads
- Uses 31-base k-mer for generating comparison keys
- Handles forward and reverse strand comparisons
- Supports short read comparison with alternative method
**Usage**: Used for sorting reads in BBTools sorting operations, particularly in read clustering and organization

## ReadComparatorCrispr (ReadComparatorCrispr.java)
**Purpose**: Comparator for sorting and ranking DNA reads with optional neural network scoring
**Core Function**: Compares DNA reads based on embedded count and score metadata
**Key Features**:
- Extracts read count from ID using parse methods
- Optional neural network-based sequence scoring
- Thread-safe neural network model loading
- Supports ascending and descending sort orders
**Usage**: Used for sorting and ranking DNA reads in bioinformatics processing pipelines

## ReadComparatorFlowcell (ReadComparatorFlowcell.java)
**Purpose**: Comparator for sorting reads based on flowcell coordinates and pair numbers
**Core Function**: Sorts Read objects by their flowcell lane, tile, y, and x coordinates
**Key Features**:
- Handles null read IDs with specific sorting rules
- Uses thread-local storage for coordinate parsing
- Supports ascending and descending sort orders
- Provides fallback to pair number comparison for equal coordinates
**Usage**: Used in sorting and organizing sequencing reads from high-throughput sequencing data

## ReadComparatorID (ReadComparatorID.java)
**Purpose**: Hierarchical comparator for sorting sequencing reads by their numeric and lexicographic identifiers
**Core Function**: Implements a three-level sorting strategy for Read objects using numeric ID, pair number, and string ID
**Key Features**:
- Sorts reads primarily by numeric ID
- Secondary sorting by read pair number
- Tertiary sorting by lexicographic string ID
- Supports ascending and descending sort order
**Usage**: Used in read sorting and organizing sequencing data, particularly for maintaining read pair relationships during processing

## ReadComparatorList (ReadComparatorList.java)
**Purpose**: Custom comparator for sorting reads based on predefined list order
**Core Function**: Maps read identifiers to ordered positions and compares reads accordingly
**Key Features**:
- Supports sorting from file or comma-separated string
- Handles reads not found in the original list
- Tie-breaks by pair number when positions are equal
- Configurable ascending/descending sort direction
**Usage**: Specialized sorting of reads in bioinformatics processing workflows

## ReadComparatorMapping (ReadComparatorMapping.java)
**Purpose**: Implements a sophisticated comparator for sorting sequencing reads based on genomic mapping characteristics.
**Core Function**: Custom comparison of Read objects for sorting, handling single-end and paired-end read scenarios with multiple prioritization criteria.
**Key Features**:
- Handles single-end and paired-end read comparisons
- Prioritizes mapped reads over unmapped reads
- Compares reads by chromosome, strand, position, and mapping quality
- Supports complex tie-breaking strategies for read ordering
**Usage**: Used in genomic read sorting and alignment processing workflows, particularly for organizing sequencing data by genomic coordinates and read characteristics

## ReadComparatorName (ReadComparatorName.java)
**Purpose**: Compares Read objects lexicographically by their identifier and pair number
**Core Function**: Implements a name-based comparator for Read objects with specific sorting rules
**Key Features**:
- Handles null identifier cases by sorting null IDs first
- Uses lexicographic comparison of Read identifiers
- Tie-breaks by pair number when identifiers are identical
- Supports ascending and descending sort directions
**Usage**: Used for sorting Read objects in sequence processing and alignment workflows

## ReadComparatorPosition (ReadComparatorPosition.java)
**Purpose**: Comparator for sorting genomic reads based on their position and metadata
**Core Function**: Provides multi-level comparison of genomic reads using scaffold, position, strand, and pair information
**Key Features**:
- Compares reads by scaffold number
- Secondary comparison by genomic position
- Handles strand orientation differences
- Supports ascending and descending sort orders
**Usage**: Used in read sorting and alignment processing workflows to organize genomic reads systematically

## ReadComparatorRandom (ReadComparatorRandom.java)
**Purpose**: Comparator for randomly ordering Read objects
**Core Function**: Compares two Read objects based on their random value attribute
**Key Features**:
- Enables random sorting of Read objects
- Supports ascending and descending sort directions
- Uses a random value comparison method
**Usage**: Used in sorting Read streams or collections with random order preference

## ReadComparatorTaxa (ReadComparatorTaxa.java)
**Purpose**: Comparator for sorting genomic reads based on their taxonomic classification
**Core Function**: Compares reads by hierarchically traversing taxonomic nodes from species to family levels
**Key Features**:
- Handles unclassified reads with name-based fallback
- Supports multi-level taxonomic node comparison
- Configurable ascending/descending sort order
- Uses TaxTree for precise taxonomic hierarchy navigation
**Usage**: Used for sorting sequencing reads based on their taxonomic assignments in bioinformatics analysis

## ReadComparatorTopological (ReadComparatorTopological.java)
**Purpose**: Performs topological comparison of genomic reads with multi-level sorting strategy.
**Core Function**: Compares Read objects using sequence, length, quality, and ID as hierarchical comparison criteria.
**Key Features**:
- Compares primary and mate sequence bases
- Handles variable-length reads with length-based comparison
- Provides quality-based sorting with inverted quality preference
- Supports ascending/descending sort direction
**Usage**: Used for sorting genomic reads in bioinformatics processing pipelines, ensuring consistent ordering of sequencing data

## ReadComparatorTopological5Bit (ReadComparatorTopological5Bit.java)
**Purpose**: Specialized comparator for sorting reads using a 5-bit topological encoding
**Core Function**: Compares reads based on k-mer, sequence content, length, and quality
**Key Features**:
- Generates 12-mer using 5-bit nucleotide encoding
- Compares reads by numeric ID, sequence content, and mate information
- Handles null and short sequences with specialized comparison methods
- Supports ascending and descending sort orders
**Usage**: Used for efficient read sorting in bioinformatics sequence processing

## ReadErrorComparator (ReadErrorComparator.java)
**Purpose**: Comparator for sorting Read objects based on error metrics
**Core Function**: Implements a multi-stage comparison of DNA reads using error count, read length, and expected error rates
**Key Features**:
- Compares total actual errors for read and its mate
- Prioritizes reads with fewer errors
- Ranks longer reads higher when error counts are equal
- Uses expected error rates as secondary sorting criterion
**Usage**: Used in sorting and filtering DNA reads during bioinformatics processing, particularly in read quality assessment

## ReadLengthComparator (ReadLengthComparator.java)
**Purpose**: Sorts reads primarily by their sequence length
**Core Function**: Custom Comparator that orders Read objects based on multiple criteria
**Key Features**:
- Sorts reads by descending length by default
- Provides secondary sorting by mate length
- Tertiary sorting by string ID
- Quaternary sorting by numeric ID
**Usage**: Used in read sorting operations where length is the primary sorting criterion

## ReadQualityComparator (ReadQualityComparator.java)
**Purpose**: Sorts reads based on quality and length, prioritizing high-quality and longer reads.
**Core Function**: Compares Read objects using a multi-stage sorting algorithm that evaluates error rates, total length, and read IDs.
**Key Features**:
- Calculates expected errors for reads and their mates
- Sorts reads with lower error rates first
- Prioritizes longer reads when error rates are equal
- Provides lexicographic and numeric ID tiebreakers
**Usage**: Used for sorting read data in sequencing analysis, ensuring highest quality reads are processed first

## Shuffle (Shuffle.java)
**Purpose**: Randomizes the order of reads.
**Core Function**: Processes input read files and modifies their order using various sorting/shuffling modes.
**Key Features**:
- Supports multiple modes: shuffle, sort by name, sequence, coordinate, or ID
- Handles both single-end and paired-end read files
- Multi-threaded input stream processing
- Configurable output formats (FASTQ, SAM/BAM)
**Usage**: Used in bioinformatics to randomize read order for statistical analysis or preprocessing

## Shuffle2 (Shuffle2.java)
**Purpose**: Sorts reads by name, potentially from multiple input files.
**Core Function**: Performs shuffling and merging of read data from FASTQ files with advanced memory management and multi-file handling
**Key Features**:
- Supports sorting reads from single or paired input files
- Handles external memory shuffling with temporary file management
- Supports deterministic shuffling using random seed
- Enables recursive merging of large file sets
**Usage**: Command-line tool for randomizing and sorting read sequence files, particularly useful in bioinformatics data preprocessing

## SortByName (SortByName.java)
**Purpose**: Sorts reads by name, potentially from multiple input files.
**Core Function**: Implements a memory-efficient read sorting mechanism that can handle multiple input files, with flexible sorting comparators and memory management.
**Key Features**:
- Supports sorting by multiple criteria: name, length, quality, position, taxonomic classification
- Uses dynamic memory allocation and temporary file management for large datasets
- Allows configurable sorting order (ascending/descending)
- Handles both single and paired-end read files
**Usage**: Used in bioinformatics workflows to organize read data before further processing, such as sequence alignment or taxonomic analysis

## SortReadsByID (SortReadsByID.java)
**Purpose**: Sorts sequencing reads by their numeric ID across multiple files and output formats
**Core Function**: Reads input files, bins reads by ID ranges, then sorts and writes output files
**Key Features**:
- Supports paired-end and single-end read sorting
- Bins reads into blocks for efficient memory management
- Handles multiple input and output file formats (FASTQ, FASTA, text)
- Configurable block size for read binning
**Usage**: Command-line tool for organizing sequencing reads by their numeric identifiers, useful for read preprocessing and data organization