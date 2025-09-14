# Assemble Package - Genome Assembly and Scaffolding Tools

## Overview
The `assemble` package contains tools for genome assembly, scaffolding, and assembly improvement.

## Core Components

### AbstractBuildThread.java - Multithreaded Assembly Worker Foundation

**Purpose**: Abstract base class providing essential infrastructure for parallel De Bruijn graph construction in Tadpole assembler, implementing thread-safe k-mer processing and contig building operations.

**Core Function**: Establishes common framework for multithreaded assembly workers with shared data structures, performance monitoring, and resource management. Coordinates concurrent read processing through synchronized input streams while maintaining thread-local storage for nucleotide frequency analysis and sequence construction.

**Threading Architecture**:

- **Parallel Execution Framework**: Extends Thread class enabling direct integration into multithreaded assembly pipelines
- **Thread Identification System**: Unique integer identifier (id) assigned at construction for coordination, debugging, and performance tracking across parallel workers
- **Concurrent Input Processing**: Array of ConcurrentReadInputStream objects supporting thread-safe access to input data across multiple workers without synchronization overhead
- **Assembly Mode Configuration**: Integer mode parameter encoding processing flags for coverage thresholds, error correction settings, and assembly behavior configuration

**Key Data Structures**:

- **Input Stream Management**: 
  - `ConcurrentReadInputStream[] crisa`: Thread-safe input stream array enabling parallel read access without blocking
  - Supports multiple input files with automatic load balancing across worker threads

- **Nucleotide Frequency Analysis**: 
  - `leftCounts[4]`: Memory-tracked array for left-side base frequency counting (A,C,G,T)
  - `rightCounts[4]`: Memory-tracked array for right-side base frequency counting (A,C,G,T)
  - KillSwitch.allocInt1D(4) allocation ensures proper memory accounting in multithreaded environment

- **Sequence Construction Infrastructure**:
  - `ByteBuilder builderT`: Thread-local sequence builder for dynamic contig construction with automatic capacity management
  - Optimized for repeated append operations during k-mer extension algorithms

- **Assembly Analytics**:
  - `LongList insertSizes`: Dynamic collection for tracking insert size statistics during paired-read processing
  - `ArrayList<Contig> contigs`: Thread-local contig storage for assembled sequences before aggregation

**Performance Monitoring System**:

- **Read Processing Metrics**:
  - `readsInT`: Thread-local counter tracking total reads processed by this worker
  - `basesInT`: Thread-local counter for total nucleotides processed

- **Quality Control Statistics**:
  - `lowqReadsT`: Counter for reads failing quality thresholds
  - `lowqBasesT`: Counter for individual bases below quality cutoffs

- **Thread-Specific Tracking**: All counters maintain thread-local state enabling parallel performance monitoring without synchronization bottlenecks during high-throughput assembly operations

**Memory Management**:

- **Resource Tracking**: KillSwitch allocation system for count arrays ensures memory usage monitoring and proper cleanup in multithreaded assembly environment
- **Thread-Local Storage**: Each worker maintains private data structures preventing synchronization overhead during intensive k-mer processing
- **Automatic Initialization**: Constructor establishes complete thread state with input streams, mode configuration, and performance counters

**Usage**: Abstract base class requiring concrete implementation of run() method for specific assembly algorithms. Provides standardized infrastructure for:

- **Parallel Assembly Workers**: Multiple threads processing different k-mer table segments simultaneously
- **Graph Construction Pipelines**: Coordinated k-mer processing with shared input streams and aggregated output
- **Performance-Critical Operations**: Thread-local data structures eliminate synchronization bottlenecks during high-throughput processing

**Thread Safety Design**: All shared resources (input streams) implement thread-safe access patterns, while thread-local resources (counters, builders, arrays) eliminate contention. The architecture supports arbitrary thread counts with linear performance scaling for k-mer processing workloads.

**Implementation Notes**: The minCountSeedCurrent field suggests dynamic threshold adjustment during assembly progression. Memory-tracked allocation via KillSwitch ensures proper resource management in long-running assembly processes. Thread identifier enables sophisticated debugging and performance analysis across parallel workers.

### AbstractRemoveThread.java
**Purpose**: Abstract worker thread framework for parallel removal of k-mers based on frequency thresholds. Implements multithreaded k-mer filtering to eliminate low-coverage k-mers (sequencing errors) and high-coverage k-mers (repetitive regions) from hash tables.

**Core Function**: Dual-threshold frequency filtering where k-mers with count < min or count > max are removed. Coordinates parallel processing across multiple hash tables using atomic work distribution.

**Key Features**:

- **Parallel Processing Architecture**: Static process() method creates multiple worker threads based on k-mer length - RemoveThread1 for short k-mers (≤31bp) and RemoveThread2 for long k-mers (>31bp)

- **Atomic Work Distribution**: Uses AtomicInteger nextTable for thread-safe table assignment. Each thread claims next available table via getAndAdd(1) operation

- **Dual Filtering Strategy**: 
  - Primary hash array filtering: Direct iteration through values array, zeroing counts outside thresholds
  - Collision resolution filtering: Recursive traversal of binary trees in victims array to handle hash collisions

- **Memory Management**: Calls clearOwnership() to reset thread flags and regenerate() to compact tables and count removed k-mers

- **Performance Monitoring**: Aggregates removal counts from all threads and provides timing information

**Usage**: Called by assembly tools to clean k-mer frequency tables before graph construction. Takes minimum and maximum frequency thresholds to define acceptable k-mer coverage range for downstream assembly algorithms.

**Implementation Details**:
- Thread safety via atomic counter starting at 0
- Polymorphic design supporting both KmerTableSet and KmerTableSetU hash table types
- Binary tree traversal for collision resolution maintains O(log n) removal complexity
- Zero-based removal marking followed by table regeneration for memory efficiency

### ErrorTracker.java - Assembly Error Detection and Correction Statistics

**Purpose**: Tracks and manages sequencing errors detected during assembly processing operations.

**Core Function**: Maintains comprehensive statistics on error detection, correction attempts, and rollback operations for quality assessment and debugging purposes. Provides centralized counter management for multiple error detection algorithms used in assembly workflows.

**Key Features**:

- **Multi-Algorithm Error Detection Tracking**: Monitors four distinct error detection methods:
  - `detectedPincer`: Pincer-based error detection count
  - `detectedTail`: Tail-based error detection count  
  - `detectedBrute`: Brute force detection algorithm count
  - `detectedReassemble`: Reassembly-based detection count

- **Correction Method Statistics**: Tracks success rates for five correction approaches:
  - `correctedPincer`: Pincer-based corrections applied
  - `correctedTail`: Tail-based corrections applied
  - `correctedBrute`: Brute force corrections applied
  - `correctedReassembleInner`: Inner reassembly corrections
  - `correctedReassembleOuter`: Outer reassembly corrections

- **Quality Assessment Metrics** (lines 77, 90, 92): Provides error position tracking with:
  - `suspected`: Count of suspected error positions
  - `marked`: Number of positions marked for correction
  - `rollback`: Boolean flag indicating rollback operations required

- **Aggregation Methods**: 
  - `corrected()`: Returns total corrections across all methods
  - `detected()`: Returns total detections (NOTE: line 53 uses `correctedTail` instead of `detectedTail` - potential bug)
  - `correctedReassemble()`: Sums inner and outer reassembly corrections

- **Counter Management**: Provides selective clearing operations:
  - `clear()`: Resets all counters and flags
  - `clearDetected()`: Resets only detection counters
  - `clearCorrected()`: Resets only correction counters

**Usage**: Integrated into assembly pipelines for real-time error tracking, quality metrics reporting via `toString()` method, and assembly algorithm performance evaluation. Essential for debugging assembly workflows and assessing correction algorithm effectiveness.

### KmerCompressor.java
**Purpose**: Compresses overlapping k-mers into contigs by identifying and extending linear chains in De Bruijn graph structures

**Core Function**: Implements k-mer graph traversal to build contigs from shared k-mer overlaps, performing atomic ownership claiming to prevent duplicate processing in multithreaded environments

**K-mer Processing Architecture**:
- **Hash Table Distribution**: Uses KmerTableSet with atomic access control via nextTable and nextVictims arrays (lines 211-216, 424-425)
- **Collision Handling**: Processes both main hash tables and overflow victim structures (HashForest) for complete k-mer coverage (lines 424-425, 444-451)
- **Atomic Ownership**: Thread-safe k-mer claiming system prevents race conditions during contig extension (lines 464-466, 691-692)
- **Bidirectional Extension**: Extends contigs in both directions using forward and reverse complement k-mers

**Contig Assembly Algorithm**:
- **Seed Selection**: Processes k-mers within count thresholds (minCount/maxCount) as potential contig seeds (lines 456, 487)
- **Linear Chain Detection**: Follows unique k-mer overlaps to extend contigs until branching points or dead ends
- **Extension Logic**: Implements fillRightCounts() to evaluate all 4 possible nucleotide extensions
- **Termination Conditions**: Stops extension at junctions (no valid next k-mer), ownership conflicts, or maximum length limits (lines 676, 711-713)

**Key Features**:
- **Coverage Filtering**: Configurable min/max k-mer count thresholds for quality control (lines 164-167, 456, 487)
- **Same-Count Extension**: Optional REQUIRE_SAME_COUNT mode for uniform coverage contigs (lines 168-169, 684, 800)
- **Reverse Complement Handling**: Dual-mode extension supporting both canonical and reverse-complement-only processing (lines 566-570, 729-837)
- **Contig Fusion**: Optional post-assembly fusion concatenating short contigs with N-base separators
- **Thread Safety**: Multithreaded processing with atomic operations and ownership tracking (lines 324-326, 428-451)

**Memory Management**:
- **ThreadLocal Storage**: Per-thread data structures (rightCounts, ByteBuilder, etc.) to avoid synchronization overhead
- **Bit Masking**: Efficient k-mer representation using bit shifts and masks for DNA encoding (lines 617, 630-631)
- **Dynamic Memory**: Configurable k-mer size with automatic memory calculation (lines 68-81, 225-232)

**Output Generation**:
- **FASTA Format**: Generates sorted contigs by length in descending order
- **Contig Naming**: Sequential numbering with optional coverage annotation
- **Statistics Reporting**: Comprehensive assembly metrics including total bases, contig count, and longest contig

**Usage**: Primary entry point for k-mer-based genome assembly. Accepts input reads, builds k-mer tables, and outputs assembled contigs:
- Input: FASTA/FASTQ files via command-line arguments
- Processing: Multithreaded k-mer loading followed by parallel contig building
- Output: Assembled contigs in FASTA format with optional assembly statistics

**Implementation Notes**: Uses atomic integer arrays for thread coordination, ensuring work distribution without locks. The extendToRight() method implements the core assembly algorithm using bit-level DNA manipulation and ownership claiming for thread safety.

### TadpoleWrapper.java - Multi-K Assembly Optimization Framework

**Purpose**: Assembles genomes with multiple k-mer lengths to systematically determine the optimal k-mer size for highest quality assembly results. Implements automated k-mer optimization through comprehensive assembly parameter exploration and quality metric evaluation.

**Core Function**: Orchestrates sequential Tadpole assembly runs across user-specified k-mer ranges, evaluating each assembly using L50, L90, maximum contig length, and total contigs metrics to identify the optimal k-mer size. Provides optional bisection and expansion algorithms for refined k-mer search within promising parameter ranges.

**Multi-K Assembly Architecture**:

- **K-mer Parameter Processing**: Parses comma-separated k-mer values and validates using Kmer.getKbig() to ensure valid odd k-mer sizes. Maintains HashSet for unique k-mer collection with automatic sorting for systematic evaluation (lines 42, 89-96).

- **Assembly Quality Evaluation**: Captures comprehensive assembly statistics from AssemblyStats2 after each Tadpole run:
  - L50: Minimum contig length containing 50% of total assembly bases
  - L90: Minimum contig length containing 90% of total assembly bases  
  - Total assembly size, contig count, and maximum contig length
  - Statistics stored in Record objects for comparative analysis

- **Systematic K-mer Testing**: Sequential assembly execution with parameter injection:
  - Dynamic argument modification inserting current k-mer size
  - Template-based output filename generation with k-mer substitution
  - Memory optimization via System.gc() before each assembly run
  - Automatic best assembly tracking with compareTo() evaluation

**Key Features**:

- **Assembly Quality Ranking System**: Multi-criteria Record comparison algorithm prioritizing:
  - **Primary Metrics**: L50 and L90 with 1% tolerance thresholds for significant differences
  - **Secondary Metrics**: Maximum contig length comparison with 1% tolerance
  - **Tie-Breaking**: Contig count minimization for equal primary metrics
  - **Fine Resolution**: 0.2% tolerance for precise ranking when coarse metrics are equal
  - **K-mer Preference**: Smaller k-mer sizes preferred for equivalent assemblies

- **Bisection Algorithm**: Recursive k-mer optimization for refined parameter selection:
  - **Midpoint Calculation**: Uses Kmer.getKbig((left.k+mid.k+1)/2) ensuring valid odd k-mers
  - **Recursive Improvement**: Continues bisection when intermediate k-mers show superior metrics (lines 265-267, 281-283)
  - **Boundary Validation**: Prevents infinite recursion by checking k-mer validity at boundaries

- **Expansion Search**: Extends k-mer search beyond initial range when edge values are optimal:
  - **Left Expansion**: Tests smaller k-mers using 0.7x scaling factor when leftmost k-mer is best
  - **Right Expansion**: Evaluates larger k-mers via min(k+40, 10+k*1.25) formula for rightmost optimum
  - **Recursive Extension**: Continues expansion while improvements are detected (lines 198-200, 223-225)

- **File Management System** (lines 44-45, 79-82, 143-175): Comprehensive output handling:
  - **Template Processing**: % symbol replacement with k-mer values in output filenames (lines 44, 62-63, 109)
  - **Final Assembly Selection**: Automatic copying of best assembly to user-specified output filename
  - **Intermediate Cleanup**: Optional deletion of suboptimal assemblies via delete parameter
  - **Cross-Platform Compatibility**: Windows-aware file operations with ReformatReads fallback

- **Early Termination Optimization** (lines 66-67, 126-129): Configurable quitEarly parameter stops k-mer testing when assembly quality begins declining, reducing computational overhead for clear optimization peaks.

**Technical Implementation**:

- **Argument Routing System**: Filters TadpoleWrapper-specific parameters while preserving all Tadpole arguments:
  - TadpoleWrapper parameters: k, out, outfinal, quitearly, delete, bisect, expand
  - All other arguments passed through to Tadpole unchanged
  - Dynamic argument list modification for k-mer and output injection (lines 100-101, 108-109)

- **Memory Management**: Strategic garbage collection before each assembly run to prevent memory accumulation across multiple Tadpole executions

- **Default Configuration**: Provides k=31 default when no k-mer values specified, ensuring functional operation without explicit parameters

**Assembly Quality Metrics**:

- **L50/L90 Analysis**: Contig length thresholds containing specified percentages of total assembly
- **Contiguity Assessment**: Maximum contig length indicates assembly contiguity quality  
- **Fragmentation Measurement**: Total contig count reflects assembly fragmentation level
- **Comprehensive Evaluation**: Multi-metric approach avoids single-criterion optimization artifacts

**Usage**: Command-line tool for automated k-mer optimization in genome assembly projects:
- **Input**: Standard Tadpole parameters plus comma-separated k-mer list via k=21,31,41,51
- **Processing**: Sequential assembly execution with quality evaluation and optional refinement
- **Output**: Best assembly in final output file with optimal k-mer size recommendation

**Performance Features**:
- **Parallel Tadpole Execution**: Each assembly run utilizes full Tadpole threading capabilities
- **Quality-Based Termination**: Early stopping prevents unnecessary computation when optimum is found
- **Incremental Improvement**: Bisection and expansion algorithms focus computation on promising parameter ranges
- **Memory Efficiency**: Garbage collection and selective file retention minimize memory footprint

**Implementation Notes**: Uses float arithmetic in Record comparison with percentage-based tolerance to handle assembly metric variability. The algorithm strongly favors L50/L90 improvements over contig count reductions, prioritizing assembly contiguity over fragmentation reduction. Template-based filename generation (% substitution) enables systematic organization of intermediate results while preserving final assembly in user-specified location.

### Contig.java - Assembly Contig Data Structure

**Purpose**: Core data structure representing assembled contigs with sequence data, coverage statistics, and graph connectivity metadata for Tadpole assembler.

**Core Function**: Stores DNA sequence as byte arrays with comprehensive assembly metadata including coverage depth, connectivity codes for graph traversal, branch ratios for decision making, and bidirectional edge lists for assembly graph navigation.

**Key Features**:
- **Sequence Storage**: DNA sequences stored as byte arrays with length calculation and GC content analysis using AminoAcid base conversion
- **Coverage Tracking**: Average, minimum, and maximum coverage values with single decimal precision formatting in headers
- **Graph Connectivity**: Left/right connectivity codes referencing Tadpole.codeStrings with branch detection and branch ratio tracking when applicable (lines 95-96, 99-100)
- **Edge Management**: Bidirectional edge lists with add/remove operations (lines 239-267, 277-293) supporting destination matching, orientation tracking, and depth merging
- **K-mer Operations**: Left/right k-mer extraction supporting both long encoding for k<32 and Kmer object population
- **Strand Operations**: Canonical strand determination via lexicographic comparison and reverse complement transformation with metadata swapping (lines 203-220, 346-371)
- **FASTA Output**: Formatted FASTA generation with configurable line wrapping and comprehensive headers including length, coverage, GC content, and connectivity information (lines 51-68, 85-111)
- **State Management**: Used/associate flags for assembly processing with edge cleanup and contig renumbering with edge reference updates

**Usage**: Primary data structure in Tadpole assembler for storing assembly results with full connectivity preservation, supporting graph traversal algorithms, bubble detection, and scaffold construction through comprehensive metadata tracking.

### Edge.java - De Bruijn Graph Edge Data Structure

**Purpose**: Represents directed edges in the De Bruijn graph connecting two contigs with orientation and overlap information.

**Core Function**: Stores connection data between contigs including orientation encoding, overlap length, coverage depth, and optional sequence data for graph traversal and assembly operations. Implements 2-bit orientation encoding for efficient storage of connection directionality.

**Graph Edge Properties**:
- **Connection Identifiers**: Source and destination contig IDs for graph connectivity (lines 26-27, 90-91)
- **Overlap Specification**: Overlap length in bases between connected contigs (lines 28, 92)
- **Orientation Encoding**: 2-bit encoding system where bit 0 represents source side and bit 1 represents destination side (lines 29, 93)
  - Value 0: Left source to left destination
  - Value 1: Right source to left destination  
  - Value 2: Left source to right destination
  - Value 3: Right source to right destination
- **Coverage Information**: Edge depth representing support evidence from read data (lines 30, 96)
- **Sequence Storage**: Optional byte array containing overlap sequence data (lines 31, 89)

**Key Features**:

- **Orientation Management**: Provides bit-level orientation testing and manipulation:
  - `destRight()`: Tests bit 1 for destination side connection
  - `sourceRight()`: Tests bit 0 for source side connection
  - `flipSource()`: Toggles source orientation bit and reverse complements sequence data
  - `flipDest()`: Toggles destination orientation bit

- **Edge Merging Operations**: Combines duplicate edges with identical topology:
  - Verifies matching origin, destination, and orientation before merging
  - Preserves data from higher-coverage edge for length and sequence
  - Sums coverage depths from both edges (lines 83, 85)

- **Graph Visualization Support**: Generates Graphviz DOT format output:
  - Creates directed edge notation with arrow syntax
  - Labels edges with orientation (LEFT/RIGHT), length, and orientation value
  - Enables graph visualization for debugging and analysis

- **String Representation**: Provides formatted output for debugging:
  - `toString()`: Delegates to ByteBuilder for efficient string construction
  - `appendTo()`: Formats edge as "(destination-orientation-length-depth-bases)"

**Sequence Consistency**: When flipping source orientation, automatically applies reverse complement transformation to stored sequence data maintaining biological accuracy.

**Usage**: Core data structure for De Bruijn graph representation in genome assembly algorithms. Enables:
- Graph construction from overlapping contigs
- Orientation-aware graph traversal for scaffold building
- Coverage-based edge filtering and merging
- Graph visualization and debugging through DOT format export

**Implementation Notes**: Uses bit manipulation for efficient orientation storage and testing. Sequence data storage is optional (may be null) allowing memory optimization for large-scale assemblies where only connectivity information is required.

### AbstractExploreThread.java - Parallel Graph Exploration and Dead-End Detection

**Purpose**: Abstract worker thread for parallel graph exploration and dead-end detection in De Bruijn graph assembly structures

**Core Function**: Implements multithreaded graph traversal algorithms to identify terminable paths and problematic graph structures during assembly processing. Provides infrastructure for concurrent detection and removal of assembly artifacts including dead ends and bubble structures.

**Graph Exploration Architecture**:
- **Dual K-mer Processing**: Maintains primary exploration k-mer and secondary comparison k-mer for graph traversal operations
- **Two-Phase Algorithm**: Sequential processing with k-mer table analysis followed by victim validation
- **Thread Coordination**: Unique thread identifier system for parallel execution coordination and statistics tracking (lines 24-25, 74)
- **Abstract Implementation**: Requires subclass implementation of core traversal methods processNextTable() and processNextVictims() (lines 57, 60)

**Dead-End Detection Protocol**:
- **Phase 1 Processing**: processNextTable() identifies candidate dead-end k-mers from hash table structures
- **Phase 2 Validation**: processNextVictims() validates and processes confirmed dead-end k-mers for removal
- **Bubble Classification**: Matrix-based analysis distinguishes genuine branches from bubble artifacts using removeMatrixT branching patterns
- **Branch Type Analysis**: F_BRANCH and B_BRANCH combination detection for bubble counting

**Key Features**:
- **Performance Tracking**: Comprehensive statistics with thread-local counters for kmersTestedT, deadEndsFoundT, and bubblesFoundT
- **Memory Management**: KillSwitch-allocated count arrays for nucleotide frequency tracking
- **Matrix Statistics**: 2D count and removal matrices for tracking branch type combinations and removal decisions
- **Dynamic Sequence Building**: ByteBuilder for constructing sequences during exploration operations
- **Thread Safety**: Thread-local data structures prevent race conditions during parallel exploration (lines 77-79, 81-82)

**Graph Structure Analysis**:
- **Branch Detection**: Maintains left and right nucleotide count arrays for identifying branch points
- **Artifact Classification**: Distinguishes dead ends from bubbles using branch pattern analysis in removal matrix
- **Termination Logic**: Abstract methods allow subclass-specific termination criteria for different graph exploration strategies (lines 57, 60)

**Usage**: Abstract base class for implementing specialized graph exploration algorithms in assembly pipelines:
- Subclasses implement processNextTable() for initial dead-end candidate identification
- Subclasses implement processNextVictims() for validated dead-end processing
- Thread management via standard start(), join(), and getState() operations
- Statistics aggregation across multiple worker threads for assembly quality metrics

**Implementation Notes**: Contains TODO indicating known issues with processNextVictims stability and nondeterministic behavior. The algorithm expects 40% dead-end rate but observes lower counts, suggesting potential issues with k-mer removal counting accuracy.

### AbstractShaveThread.java

**Purpose**: Abstract worker thread framework for parallel removal of dead-end k-mers from assembly graphs.

**Core Function**: Implements multithreaded graph cleaning operations through coordinated dead-end trimming to remove spurious branches and improve assembly quality. The thread continuously processes work units by polling a shared queue until empty.

**Key Features**:
- **Thread Pooling Architecture**: Each thread maintains unique identifier for logging and statistics tracking (lines 17-18, 32)
- **Work Queue Processing**: Implements polling-based work distribution with `processNextTable()` method called continuously until queue exhaustion
- **K-mer Removal Tracking**: Maintains per-thread counter `kmersRemovedT` for monitoring cleanup progress
- **Template Method Pattern**: Abstract `processNextTable()` allows concrete implementations to define specific graph trimming algorithms

**Usage**: Serves as base class for specialized shaving thread implementations. Concrete subclasses implement `processNextTable()` to define specific dead-end detection and removal strategies. Multiple threads coordinate through shared work queue to parallelize graph cleaning operations across large assembly datasets.

**Technical Implementation**: Thread runs continuously until work queue is depleted, with each iteration processing one table/work unit through the abstract method interface.

### BubblePopper.java - Assembly Graph Bubble Detection and Simplification

**Purpose**: Detects and collapses bubbles in assembly graphs to simplify complex branching structures caused by sequencing errors, heterozygosity, or repeats (lines 9-16, 18-46).

**Core Function**: Implements sophisticated bubble detection algorithms to identify alternative paths between shared endpoints in contig graphs, then merges these paths to create linear assemblies while preserving sequence accuracy and coverage information (lines 154-308, 311-393).

**Bubble Detection Architecture**:
- **Direct Linear Extension**: Simple path merging for single-edge connections using expandRightSimple() algorithm
- **Complex Bubble Detection**: Multi-step validation for branching structures via expandRight() with mutual destination analysis
- **Bidirectional Processing**: Processes both left and right directions by flipping contigs for comprehensive bubble detection
- **Representative Node Selection**: Chooses optimal intermediate nodes based on length and depth criteria using findRepresentativeMidEdge()

**Key Features**:

- **Minimum Length Filtering**: Enforces minLen = 2*kbig-1 to ensure reliable overlap detection for k-mer based assemblies
- **Multi-Algorithm Approach**: Combines direct merging (popDirect) and indirect bubble detection (popIndirect) for comprehensive graph simplification
- **Dead-End Branch Removal**: Implements debranching functionality to remove short dead-end paths based on:
  - Length thresholds (400bp cutoff at lines 65, 79)
  - Connection analysis for left/right sides
  - Quality metrics combining length and coverage

- **Topology Validation**: Ensures intermediate nodes have consistent connectivity patterns:
  - All mid-nodes must connect to same left destination
  - All mid-nodes must connect to same right destination
  - No self-loops allowed in bubble structure (lines 610-614, 621-625)

- **Mutual Destination Finding**: Core algorithm for identifying bubble convergence points:
  - Traces paths through intermediate nodes to find common endpoints
  - Validates orientation consistency across all paths
  - Returns -1 for invalid topologies, >=0 for valid convergence points

**Sequence Merging Algorithm**:
- **Path Concatenation** (lines 325-340, 446-453): Builds merged sequences by:
  - Appending left contig bases completely
  - Adding bridge sequences from edges (excluding final k-mer overlap)
  - Including intermediate node internal sequences (excluding k-mer overlaps at both ends)
  - Concatenating right contig bases
- **Coverage Calculation** (lines 369-370, 477-478): Weighted average based on unique base contributions from each contig
- **Edge Redirection**: Updates graph topology by redirecting incoming edges to merged contig

**Graph Maintenance**:
- **Loop Detection**: Identifies self-looping contigs by checking edge connectivity patterns
- **Dead Edge Cleanup**: Removes edges pointing to consumed/associated contigs
- **Validation Framework**: Comprehensive integrity checking for graph consistency (disabled for performance)

**Usage**: Integrated into assembly pipelines for graph simplification after initial contig construction:
- Input: Contig collection with edge connectivity map (destMap) and k-mer size parameter
- Processing: Iterative bubble detection and merging via expand() method
- Output: Simplified graph with merged contigs and updated connectivity

**Performance Metrics**:
- `expansions`: Count of successful bubble merges (lines 372, 486)
- `contigsAbsorbed`: Total contigs consumed during merging (lines 373, 487)
- `branchesRemoved`: Dead-end branches eliminated during debranching

**Implementation Notes**: Uses ByteBuilder for efficient sequence concatenation (lines 321-340, 442-453) and atomic contig state management via setUsed()/setAssociate() methods to prevent double-processing in multithreaded environments. The algorithm prioritizes longer, higher-coverage paths when multiple alternatives exist, preserving assembly quality while reducing complexity.

### AbstractProcessContigThread.java - Parallel Contig Connectivity Analysis Framework

**Purpose**: Abstract worker thread foundation for parallel processing of contig connectivity and extension operations in assembly graph traversal algorithms

**Core Function**: Implements multithreaded framework for exploring connections between assembled contigs to identify opportunities for further extension and merging through bidirectional connectivity analysis (lines 11-12, 35-41)

**Threading Architecture**:
- **Atomic Work Distribution**: Uses AtomicInteger counter for lock-free work allocation across threads (lines 25-27, 36)
- **Abstract Processing Model**: Defines template for left/right connectivity analysis requiring concrete implementation
- **Thread-Safe Operation**: Extends Thread class with shared contig list and thread-local processing arrays (lines 18, 25-28)

**Key Features**:

- **Bidirectional Connectivity Analysis** (lines 38-39, 43-45): 
  - `processContigLeft()`: Abstract method for analyzing left-side connections
  - `processContigRight()`: Abstract method for analyzing right-side connections
  - Each direction processed independently with dedicated count arrays

- **Connection Type Tracking**: Three count arrays with size 4 each for categorizing connectivity:
  - `leftCounts`: Tracks different left-side connection types
  - `rightCounts`: Tracks different right-side connection types 
  - `extraCounts`: Extended connectivity metrics for complex analysis

- **Performance Optimization Cache**: Thread-local variables for avoiding redundant calculations:
  - `lastLength`: Cached length value for performance optimization
  - `lastTarget`: Cached target identifier
  - `lastExitCondition`: Tracks processing termination reason
  - `lastOrientation`: Forward(0) or reverse-complement(1) orientation tracking

- **Memory Management**: 
  - **KillSwitch Allocation**: Uses KillSwitch.allocInt1D(4) for memory-tracked array allocation
  - **ByteBuilder Reuse**: Thread-local ByteBuilder instance avoids repeated memory allocation
  - **Statistics Tracking**: Thread-local edge creation counter for connectivity metrics

**Processing Algorithm**:
- **Work Distribution Loop**: Processes contigs using atomic counter increment for thread coordination
- **Sequential Processing**: Each thread processes assigned contigs in order, analyzing both left and right connectivity
- **Termination Handling**: Tracks exit conditions and orientations for debugging and optimization

**Usage**: Abstract base class for implementing specific contig connectivity algorithms:
- Subclasses must implement `processContigLeft()` and `processContigRight()` methods
- Provides shared infrastructure for parallel graph traversal operations
- Used in assembly pipelines requiring bidirectional contig analysis

**Implementation Notes**: The abstract design allows specialization for different connectivity algorithms while maintaining consistent multithreaded processing patterns. The size-4 arrays suggest tracking of DNA nucleotides (A,C,G,T) or connection states, with thread-local caching optimizing repeated graph operations. Edge statistics enable performance monitoring across parallel workers.

### Shaver.java - K-mer Graph Trimming and Dead-End Removal System

**Purpose**: Abstract factory and coordination system for multithreaded removal of dead-end structures (aka hairs) and bubble artifacts from k-mer-based assembly graphs (lines 12, 17). Implements sophisticated graph trimming algorithms to improve assembly contiguity by eliminating spurious branches caused by sequencing errors or low-coverage regions.

**Core Function**: Orchestrates two-phase parallel processing with exploration phase identifying removable k-mer structures followed by atomic removal phase. Uses dynamic dispatch based on k-mer table implementation to instantiate appropriate worker thread types. Coordinates thread-safe access to hash tables through atomic counters while maintaining comprehensive statistics for algorithm evaluation.

**Key Features**:

- **Factory Design Pattern**: Dynamic thread instantiation based on k-mer table type:
  - `KmerTableSet` → `Shaver1` implementation for standard k-mer processing
  - `KmerTableSetU` → `Shaver2` implementation for unsigned k-mer handling
  - Conservative default parameters: low thresholds with both hair and bubble removal enabled

- **Dual-Phase Shaving Algorithm**:
  - **Phase 1 - Exploration**: AbstractExploreThread workers identify dead-end candidates and bubble structures using graph traversal
  - **Phase 2 - Removal**: AbstractShaveThread workers perform atomic k-mer removal from hash tables
  - Independent timing and statistics tracking for each phase (lines 82, 130, 136, 166-167)

- **Multithreaded Coordination** (lines 100-125, 146-161):
  - Atomic work distribution via `nextTable` and `nextVictims` counters (lines 96-97, 142-143, 223-226)
  - Thread lifecycle management with complete termination verification (lines 106-113, 151-158)
  - Cross-thread statistics aggregation using matrix accumulation

- **Coverage-Based Filtering** (lines 47-62, 71-75):
  - `minCount`/`maxCount`: K-mer frequency bounds for inclusion in shaving operations
  - `minSeed`: Minimum coverage for seed k-mer selection during exploration
  - `minCountExtend`: Coverage threshold for extension during graph traversal
  - `branchMult2`: Coverage ratio threshold for branch removal decisions

- **Graph Structure Parameters**:
  - `maxLengthToDiscard`: Maximum path length eligible for removal to prevent over-trimming
  - `maxDistanceToExplore`: Search radius limit during dead-end detection to control complexity
  - `removeHair`/`removeBubbles`: Boolean flags enabling specific artifact removal types

- **Performance Optimization Flags**:
  - `startFromHighCounts=true`: Begin exploration from high-coverage k-mers for faster convergence
  - `shaveFast=true`: Enable optimized trimming algorithms
  - `shaveVFast=false`: Ultra-fast mode with potential contiguity trade-offs

**Threading Architecture**:

- **Abstract Thread Factories**: Subclasses must implement `makeExploreThread()` and `makeShaveThread()` for k-mer table-specific operations
- **Table Ownership Management**: Initializes thread ownership flags for atomic k-mer claiming during processing
- **Event Matrix Tracking**: 8x8 matrices recording exploration and removal event types for algorithm analysis
- **Thread-Safe Statistics**: Accumulates per-thread counters: kmersTestedT, deadEndsFoundT, bubblesFoundT, kmersRemovedT

**Usage**: Primary interface for graph cleaning in assembly pipelines:
- Basic usage: `Shaver.makeShaver(tables, threads)` with conservative defaults
- Advanced usage: Custom parameters for coverage thresholds, distance limits, and removal type selection
- Runtime parameter adjustment via `shave(minCount, maxCount)` method

**Algorithm Validation**: Enforces parameter consistency with assertion checking:
- Requires: `minCount ≤ minSeed ≤ maxCount` for logical coverage bounds
- Ensures at least one removal type (hair or bubbles) is enabled

**Statistical Output** (lines 132-167, 170-185):
- Exploration metrics: k-mers tested, dead ends found, bubbles detected
- Removal metrics: k-mers removed with timing information
- Optional event matrices: 8x8 tables showing detailed operation breakdown

**Implementation Notes**: Uses atomic integer coordination for lock-free thread synchronization. The 8x8 event matrices enable sophisticated analysis of shaving algorithm behavior, tracking specific combinations of graph topology and removal decisions. Conservative defaults ensure safe operation while advanced parameters allow tuning for specific datasets or assembly quality requirements.

### Tadpole1.java - Short K-mer De Bruijn Graph Assembler

**Purpose**: Short k-mer assembler optimized for k-mers ≤31 bases based on KmerCountExact algorithm. Implements high-performance De Bruijn graph assembly with exact k-mer counting and multithreaded contig construction for small to medium genome assembly projects.

**Core Function**: Builds De Bruijn graphs from exact k-mer frequency tables using KmerTableSet hash structures, then performs guided graph traversal to construct contigs through bidirectional extension with ownership claiming for thread safety (lines 70-77, 113-116, 272-285).

**Key Features**:

- **Memory-Optimized Storage**: Dynamic memory allocation based on assembly mode:
  - Base storage: 12 bytes per k-mer for value, count, and collision handling
  - Thread ownership: +4 bytes when useOwnership enabled for parallel processing
  - Assembly metadata: +1 byte for contig/extend modes to track graph state

- **Multithreaded Contig Building**: Parallel BuildThread workers with atomic work distribution:
  - Multiple assembly passes with decreasing coverage thresholds
  - Atomic table and victim processing via nextTable/nextVictims arrays (lines 176-177, 182-183)
  - Thread ownership claiming prevents duplicate processing of k-mers (lines 230-237, 261-268)

- **Progressive Assembly Strategy**: Multi-pass algorithm with adaptive thresholds:
  - Initial passes use elevated coverage requirements (minCountSeed + i)
  - Geometric decay factor (contigPassMult) reduces thresholds progressively
  - Final pass uses minimum seed threshold for comprehensive assembly

- **Bidirectional Extension Algorithm**: Sophisticated contig building from k-mer seeds:
  - Forward extension with extendToRight() exploring all possible nucleotide additions
  - Reverse complement processing for bidirectional growth
  - Junction detection and branch ratio calculation for assembly decision making (lines 424-431, 460-467)
  - Dynamic termination handling for loops, dead ends, and ownership conflicts

- **Graph Traversal and Connectivity**: Advanced graph exploration for scaffold construction:
  - Left/right k-mer connectivity analysis using fillLeftCounts/fillRightCounts (lines 637-641, 674-677)
  - Junction threshold analysis with isJunction() for branch point detection (lines 646, 683)
  - Path exploration up to 500 bases with exploreRight() algorithm
  - Orientation tracking for proper contig joining (0=left-to-left, 2=left-to-right)

- **Coverage-Based Quality Control** (lines 275-276, 480-481): Multiple filtering layers:
  - K-mer frequency thresholds (minCountSeed, minCountExtend) for error removal
  - Contig coverage bounds (minCoverage, maxCoverage) calculated by tables.calcCoverage()
  - Minimum contig length requirements (minContigLen) with extension validation
  - Branch ratio analysis for complex junction resolution

**Assembly Algorithm Details**:

- **K-mer Processing Pipeline**: Cell-by-cell hash table traversal:
  - Coverage threshold filtering eliminates low-quality k-mers
  - Thread ownership claiming via setOwner() prevents race conditions
  - Contig construction from high-confidence k-mer seeds
  - Atomic contig ID assignment and storage

- **Extension Termination Logic**: Comprehensive stopping criteria:
  - BAD_SEED: Insufficient k-mer coverage for reliable extension (lines 901, 924, 932-935)
  - DEAD_END: No valid next k-mers available for continuation (lines 971, 1078-1081)
  - LOOP: Self-intersection detected via ownership checking
  - Branch codes (F_BRANCH, B_BRANCH, D_BRANCH): Junction points requiring decision logic (lines 972-980, 1039-1051)

- **Read Extension Mode**: Alternative assembly from existing reads:
  - Insert size estimation using measureInsert() path finding
  - Read-based extension with extendRead() for gap filling
  - Quality score propagation for extended bases using FAKE_QUAL

**Error Correction Integration**: Comprehensive error detection and correction:
- **Multi-Algorithm Approach**: Combines pincer, tail, brute force, and reassembly methods
- **Quality-Based Validation**: Uses expected error rates and quality scores for correction decisions (lines 1450, 1496-1500)
- **Rollback Protection**: ECC_ROLLBACK system prevents over-correction damage
- **Base Quality Marking**: MARK_BAD_BASES flag system for problematic position tracking

**Junk Detection System**: Read quality assessment:
- **K-mer Presence Analysis**: Detects reads with insufficient k-mer support
- **Coverage Threshold Testing**: hasKmersAtOrBelow() identifies low-coverage regions
- **Paired-Read Awareness**: Adjusts quality thresholds for paired vs single reads (lines 1302, 1320-1323)

**Usage**: Primary assembler for projects requiring high accuracy with moderate computational resources:
- Input: FASTA/FASTQ files with reads for k-mer table construction
- Processing: Multithreaded k-mer loading followed by parallel contig building
- Output: High-quality contigs with coverage metadata and connectivity information
- Optimal for: Bacterial genomes, small eukaryotes, targeted assemblies with k ≤ 31

**Implementation Notes**: Extends abstract Tadpole class providing KmerTableSet-specific implementations. Uses bit manipulation for efficient k-mer encoding with shift operations and mask values. The architecture supports both assembly modes (contig construction) and correction modes (read error fixing) through unified k-mer frequency analysis. Thread safety achieved through atomic operations and ownership claiming rather than traditional locking mechanisms.

### Postfilter.java - Two-Phase Assembly Quality Filtering Pipeline

**Purpose**: Implements post-assembly quality filtering through sequential read alignment and coverage-based contig filtering to produce high-quality assembly outputs.

**Core Function**: Executes two-phase pipeline combining BBMap read alignment with FilterByCoverage contig filtering. Phase 1 aligns input reads against assembly contigs to generate coverage statistics. Phase 2 applies multiple quality thresholds to retain only well-supported contigs.

**Key Features**:

- **Two-Phase Processing Architecture**:
  - **Phase 1 - Coverage Analysis**: BBMap alignment with permissive settings for comprehensive coverage calculation
  - **Phase 2 - Quality Filtering**: FilterByCoverage application with configurable thresholds for final contig selection
  - **Optional Statistics Generation**: AssemblyStats2 integration for final assembly quality metrics

- **Comprehensive Quality Thresholds**:
  - `minCoverage=2.0`: Minimum average coverage depth across contig length (lines 96, 171, 207)
  - `minReads=6`: Minimum number of supporting reads for contig validation (lines 100, 172, 209)
  - `minCoveredPercent=95.0`: Minimum percentage of bases with coverage above threshold (lines 98, 173, 208)
  - `minLength=400`: Minimum contig length in bases for inclusion (lines 102, 174, 210)

- **Alignment Parameter Optimization**:
  - `ambig=all`: Processes all ambiguous alignments for maximum coverage capture
  - `maxindel=0`: Strict indel tolerance for accurate coverage measurement (lines 91, 156, 199)
  - `minhits=2`: Minimum seed hits requirement for alignment validation (lines 93, 155, 200)
  - `nodisk`: Memory-only operation for improved performance

- **Flexible Input/Output Configuration** (lines 77-86, 192-197):
  - Paired-end read support with automatic file detection using # placeholder
  - Primary output for high-quality contigs (lines 83-84, 195)
  - Optional `outdirty` parameter for capturing filtered-out contigs (lines 85-86, 169, 196)
  - Configurable coverage statistics file output (lines 89-90, 197)

- **Optional End Trimming** (lines 105-111, 175, 205):
  - Configurable base trimming from contig ends for quality improvement
  - Boolean mode (trim=true sets 100bp default) or explicit base count specification
  - Applied during filtering phase to remove potentially problematic terminal regions

**Usage**: Command-line tool for assembly post-processing requiring input reads, reference contigs, and output paths:
```bash
# Basic filtering with default thresholds
java Postfilter in=reads.fq ref=contigs.fa out=filtered.fa covstats=coverage.txt

# Custom quality thresholds
java Postfilter in=reads.fq ref=contigs.fa out=filtered.fa mincov=5 minr=10 minp=90 minl=500

# Paired-end input with rejected contig capture
java Postfilter in=reads#.fq ref=contigs.fa out=clean.fa outdirty=rejected.fa
```

**Implementation Details**:
- **File Validation Framework**: Comprehensive input/output file validation with overwrite protection and duplicate detection
- **Argument Processing**: Dynamic argument parsing with automatic forwarding of unrecognized parameters to underlying tools
- **Memory Management**: Explicit data unloading between phases via Data.unloadAll() to prevent memory accumulation
- **Statistics Integration**: Optional AssemblyStats2 invocation provides N50, total length, and contig distribution metrics

**Quality Control Pipeline**: The filtering cascade removes contigs failing any threshold, ensuring output contains only high-confidence assembly results. The permissive alignment phase (ambig=all) maximizes coverage capture, while strict filtering phase eliminates poorly-supported sequences. This approach balances sensitivity during coverage calculation with specificity during final selection.

### TadPipe.java - Integrated Assembly Pipeline Coordinator

**Purpose**: Orchestrates a complete assembly pipeline by sequentially executing read preprocessing, error correction, merging, and multi-k assembly operations through coordinated BBTools component invocation (lines 19-24, 113-396).

**Core Function**: Implements a standardized assembly workflow that processes raw sequencing reads through a series of quality improvement steps before final assembly, managing intermediate file creation, parameter routing, and temporary file cleanup for streamlined genome assembly operations.

**Pipeline Architecture**:

- **Sequential Tool Coordination**: Executes seven distinct processing stages in fixed order:
  1. **Adapter/Quality Trimming**: BBDuk with default adapters, k=23/mink=11, quality trimming at Q10
  2. **ECCO Error Correction**: BBMerge in strict ECCO mode with default adapters
  3. **Read Clumpification**: Clumpify with error correction, 8 passes, unpair/repair operations
  4. **Read Merging**: BBMerge with k=75, extend2=120, ordered output
  5. **K-mer Error Correction**: Tadpole ECC with k=50, deadzone=2, junk removal
  6. **Read Extension**: Tadpole extend mode with k=81, 100bp left/right extension
  7. **Optional Second Extension**: Tadpole k=124 with 60bp extension when extend2=true
  8. **Multi-K Assembly**: TadpoleWrapper with k=210,250,290 and graph optimization

- **Parameter Routing System**: Prefix-based argument distribution to appropriate tools:
  - `merge_*`: Routed to BBMerge operations
  - `ecco_*`: Directed to ECCO error correction
  - `ecc_*` / `correct_*`: Sent to Tadpole error correction
  - `extend_*` / `extend1_*`: First extension pass parameters
  - `extend2_*`: Second extension pass parameters
  - `clump_*` / `clumpify_*`: Clumpify configuration
  - `trim_*`: BBDuk trimming parameters
  - `assemble_*`: Final assembly parameters

**Key Features**:

- **Temporary File Management**: Creates intermediate files with proper extensions based on compression settings:
  - Trimmed reads, ECCO output, clumped reads
  - Merged/unmerged streams, error-corrected versions
  - Extended reads (first and optional second pass)
  - Multi-k assembly intermediates with template pattern

- **Compression Control**: Global gz parameter controls intermediate file compression with .fq.gz vs .fq extensions (lines 80-81, 140)

- **Resource Optimization**: Configures BBTools components for optimal performance:
  - Enables parallel gzip (PIGZ/UNPIGZ) for faster I/O
  - Allocates half available threads to compression operations
  - Disables automatic quality detection for processing speed

- **Cleanup Management**: Selective temporary file deletion controlled by deleteTemp parameter:
  - Preserves original input files (lines 251-252, 407)
  - Removes intermediate files after each processing stage (lines 274, 297, 314, 334, 354, 378)
  - Safe deletion with existence verification

- **Intelligent Pipeline Adaptation**: 
  - Skips adapter trimming if adapter database unavailable (lines 195, 232-234)
  - Optional second extension pass based on extend2 flag (lines 154-157, 337-355)
  - TadpoleWrapper selects optimal k-mer size from multi-k assembly results

**Technical Implementation**:

- **Argument Processing**: Standard BBTools key=value parsing with lowercase normalization
- **Tool Integration**: Direct main() method invocation of BBTools components with reversed argument arrays for proper precedence (lines 228, 249, 272, 295, 313, 333, 353, 377)
- **File Safety**: Comprehensive path validation preventing accidental deletion of input files
- **Error Handling**: IOException catching for temporary file creation with stack trace output

**Default Processing Parameters**:
- **BBDuk Trimming**: k=23, mink=11, hdist=1, right k-mer trimming, Q10 quality trimming, 62bp minimum length
- **ECCO Correction**: Strict mode, mixed reads, default adapter handling, ordered output
- **Clumpification**: 8 error correction passes, short naming, unpair/repair operations
- **BBMerge**: k=75 overlap detection, 120bp extension, adapter removal, ordered processing
- **Tadpole ECC**: k=50, deadzone=2, junk removal, ordered output
- **Tadpole Extension**: k=81, 100bp bidirectional extension, deadzone=0
- **TadpoleWrapper Assembly**: k=210,250,290 multi-k approach with expand, bisect, shave, rinse, pop optimization

**Usage**: Command-line assembly pipeline accepting standard BBTools parameters with prefix routing for component-specific configuration. Primary input via in=/in1= and out= parameters, with automatic intermediate file management and optional temporary directory specification via tmpdir= parameter.

**Implementation Notes**: Uses Collections.reverse() on argument lists before tool invocation to ensure proper BBTools parameter precedence (lines 228, 247, 270, 293, 311, 331, 351, 375). The pipeline design assumes adequate disk space for multiple intermediate files and benefits from SSD storage for intensive I/O operations. TadpoleWrapper integration enables automatic k-mer size optimization for best assembly results.

### ShaveObject.java - Assembly Graph Operation Base Class and Constants

**Purpose**: Abstract base class providing shared constants and common functionality for assembly graph processing throughout the Tadpole assembly pipeline

**Core Function**: Defines standardized operation modes, exploration result codes, and branch type classifications to ensure consistent graph traversal and manipulation across all assembly components (lines 8-9, 23-32)

**Key Features**:

- **Assembly Mode Constants**: Five distinct operation modes for different assembly contexts:
  - `contigMode=0`: Standard contig assembly from input reads
  - `extendMode=1`: Extension of existing contig sequences
  - `correctMode=2`: Error correction processing mode
  - `insertMode=3`: Gap filling between existing contigs
  - `discardMode=4`: Low-quality sequence removal operations

- **Graph Exploration Result Codes**: Seven standardized outcome codes for graph traversal operations:
  - `KEEP_GOING=0`: Continue graph exploration (normal progress)
  - `DEAD_END=1`: Path terminates without connections
  - `TOO_SHORT=2`: Path length below minimum threshold
  - `TOO_LONG=3`: Path length exceeds maximum threshold
  - `TOO_DEEP=4`: Recursion depth limit reached
  - `LOOP=7`: Circular path detected during traversal
  - `SUCCESS=8`: Target reached successfully

- **Branch Type Classification** (lines 32, 34): Bit-masked branch detection system for graph topology analysis:
  - `BRANCH_BIT=16`: Bit mask for branch code detection
  - `F_BRANCH=17`: Forward branch (BRANCH_BIT|1)
  - `B_BRANCH=18`: Backward branch (BRANCH_BIT|2)  
  - `D_BRANCH=19`: Dead-end branch (BRANCH_BIT|3)
  - `isBranchCode()`: Bit mask test function for branch identification

- **Extension Error Codes**: Specialized error conditions for contig extension operations:
  - `BAD_OWNER=11`: Ownership conflict during extension
  - `BAD_SEED=12`: Invalid seed sequence for extension

- **Traversal State Management**: Four-state progression system for graph exploration:
  - `STATUS_UNEXPLORED=0`: Initial unprocessed state
  - `STATUS_EXPLORED=1`: Successfully analyzed
  - `STATUS_REMOVE=2`: Marked for deletion
  - `STATUS_KEEP=3`: Confirmed for retention

- **Diagnostic Framework** (lines 41-46, 50, 53-55): Comprehensive debugging infrastructure:
  - `codeStrings[]`: Human-readable names for all numeric codes
  - `MAX_CODE`: Array bounds safety constant
  - `printEventCounts`: Performance monitoring toggle
  - `verbose`/`verbose2`: Two-level diagnostic output control

- **Logging Infrastructure**: Centralized error stream output for assembly operations using `System.err` for consistent debugging across all components

**Usage**: Inherited by all assembly graph processing classes to ensure consistent:
- Operation mode handling across different assembly contexts
- Standardized result code interpretation for graph algorithms
- Uniform branch type classification for topology analysis
- Shared diagnostic and logging functionality

**Implementation Notes**: Uses bit manipulation for efficient branch detection with `BRANCH_BIT` mask, enabling rapid topology classification during graph traversal. The comprehensive constant definitions eliminate magic numbers throughout the assembly codebase, improving maintainability and reducing errors in graph processing algorithms.

### Shaver2.java - Advanced K-mer Graph Trimming with Unlimited K-mer Support

**Purpose**: Advanced dead-end removal and bubble elimination tool for De Bruijn graph cleaning with support for unlimited k-mer lengths through KmerTableSetU data structures (lines 17-22, 39).

**Core Function**: Implements sophisticated graph shaving algorithms that detect and remove dead-end paths and bubble artifacts from k-mer graphs using bidirectional exploration and path classification. Provides thread-safe multithreaded processing with atomic work distribution for scalable performance on large assemblies (lines 23, 31-41).

**Advanced Shaving Architecture**:

- **Unlimited K-mer Support**: Extends base Shaver class with KmerTableSetU integration for arbitrary k-mer lengths beyond 31 bases, enabling complex genome assemblies requiring longer k-mers
- **Dual-Phase Processing**: Sequential exploration followed by targeted removal using ExploreThread and ShaveThread implementations (lines 48-50, 309-595)
- **Bidirectional Path Analysis**: Complete path exploration in both directions using reverse complement operations for comprehensive dead-end detection

**Key Features**:

- **Path Classification System**: Sophisticated termination code analysis distinguishing:
  - `DEAD_END`: Single-ended paths suitable for removal
  - `B_BRANCH`: Backward branching points indicating complex connectivity (lines 126, 136-142)
  - `F_BRANCH`: Forward branching points requiring preservation (lines 110, 115)
  - `TOO_LONG`, `TOO_DEEP`, `LOOP`: Problematic structures requiring exploration termination (lines 110-118, 120-123)

- **Configurable Removal Strategy**: Dual removal modes supporting:
  - **Dead-End Removal**: Eliminates simple dead-end paths (DEAD_END to DEAD_END or DEAD_END to B_BRANCH) when removeHair=true
  - **Bubble Removal**: Collapses B_BRANCH to B_BRANCH connections when removeBubbles=true

- **Graph Exploration Algorithm**: Advanced path traversal implementing:
  - **Seed K-mer Validation**: Ensures k-mer population and count threshold compliance (lines 82-85, 173-174)
  - **Dynamic Extension**: Continuous path building using nucleotide frequency analysis
  - **Loop Detection**: First k-mer comparison to prevent infinite cycles
  - **Branch Point Recognition**: Multi-count threshold analysis for identifying topology changes (lines 186-196, 247-275)

**Multithreaded Processing Framework**:

- **ExploreThread Implementation**: Specialized worker for dead-end identification:
  - **Hash Table Processing**: Atomic table distribution with count-based filtering
  - **Collision Resolution**: Binary tree traversal through victim structures
  - **High-Count Optimization**: Dedicated processing path for high-coverage k-mers with branch analysis acceleration
  - **Left-Side Extension Testing**: Four-nucleotide extension testing for comprehensive branch detection

- **ShaveThread Implementation**: Targeted k-mer removal worker:
  - **Atomic Removal Operations**: Thread-safe k-mer count zeroing with STATUS_REMOVE coordination
  - **Hash Forest Cleanup**: Recursive binary tree traversal for collision victim removal
  - **Table Regeneration**: Ownership clearing and hash table compaction after removal operations

**Performance Optimization Features**:

- **shaveVFast Acceleration**: Optional branch pre-filtering reducing unnecessary exploration by analyzing left-side extensions before full path traversal
- **Count-Based Filtering**: Dual threshold processing (processCell_low/processCell_high) optimizing exploration for different coverage ranges
- **Atomic Work Distribution**: Lock-free table assignment preventing thread contention during parallel processing (lines 320, 349, 565)

**Graph Integrity Management**:

- **Ownership Claiming System**: Thread-safe k-mer ownership with atomic claiming operations preventing race conditions during exploration
- **Count Validation**: Continuous threshold checking ensuring processed k-mers remain within assembly quality bounds
- **Path Validation**: ByteBuilder-based sequence construction with k-mer consistency verification throughout exploration

**Memory Management**:
- **KmerTableSetU Integration**: Direct access to unlimited k-mer hash tables with collision handling (lines 169, 214, 322, 533)
- **Thread-Local Storage**: Per-thread count arrays and sequence builders preventing synchronization overhead
- **Efficient K-mer Operations**: Direct table access methods for count retrieval, ownership management, and nucleotide frequency analysis

**Usage**: Integrated into assembly pipelines for advanced graph cleaning after initial k-mer table construction:
- Input: Populated KmerTableSetU with assembled k-mer graph structure
- Configuration: Threshold parameters (minCount, maxCount, minSeed), removal flags (removeHair, removeBubbles), and performance limits (maxLengthToDiscard, maxDistanceToExplore)
- Processing: Parallel exploration phase followed by coordinated removal phase
- Output: Cleaned k-mer tables with dead-ends and bubbles eliminated

**Implementation Notes**: Uses unlimited k-mer encoding through KmerTableSetU for complex assemblies requiring k-mers >31 bases. The bidirectional exploration algorithm ensures complete path characterization before removal decisions. Advanced branch detection with configurable stringency balances assembly contiguity against false positive removal. Thread safety achieved through atomic operations and per-thread data structures without explicit locking.

### Tadpole2.java - Long K-mer Assembly Engine with Advanced Error Correction

**Purpose**: High-performance long k-mer assembler implementation extending Tadpole with sophisticated error correction algorithms and multithreaded contig construction for complex genome assembly tasks.

**Core Function**: Implements De Bruijn graph-based assembly using long k-mers (>31bp) stored in KmerTableSetU hash structures with comprehensive error correction capabilities including pincer-based detection, tail extension analysis, and reassembly validation algorithms (lines 70-76, 1367-1387).

**K-mer Table Architecture**:
- **Long K-mer Support**: Uses KmerTableSetU for k-mers exceeding 31bp limit with configurable extra bytes per k-mer for metadata storage (lines 61-68, 70)
- **Ownership System**: Thread-safe k-mer claiming mechanism prevents duplicate processing during multithreaded assembly operations (lines 64, 89-91, 122-132)
- **Dual Hash Structure**: Processes both main hash arrays and collision victim trees (HashForestU) for complete k-mer coverage (lines 194-219, 244-275)
- **Progressive Depth Parameters**: Multi-pass shaving algorithm with increasing depth parameters (a=1, b=maxShaveDepth, c=i+1) for iterative graph cleaning

**Assembly Construction Algorithm**:
- **Multi-Pass Seeding**: Implements progressive minimum count reduction across multiple passes:
  - Initial high-stringency passes with elevated minimum counts calculated via exponential decay formula
  - Final pass using base minCountSeed for comprehensive coverage
  - Each pass processes both hash tables and victim structures completely (lines 173-174, 179-180)

- **Bidirectional Extension**: Core contig building through two-phase extension:
  - **Phase 1**: Forward extension from seed k-mer with ownership claiming and junction detection
  - **Phase 2**: Reverse complement extension after sequence reversal for full contig construction
  - **Termination Handling**: Multiple exit conditions including DEAD_END, LOOP, BAD_SEED, BAD_OWNER, and branch codes (lines 391-424, 438-460)

- **Coverage-Based Filtering**: Applies minimum/maximum coverage thresholds with dynamic coverage calculation for quality control

**Advanced Error Correction Framework**:
- **Multi-Algorithm Approach**: Implements five distinct error detection and correction methods:
  - **Pincer Correction**: Bidirectional validation using flanking k-mer analysis for single-base substitutions (lines 1368, 1434-1484)
  - **Tail Extension**: Unidirectional correction with progressive extension validation (lines 1371-1381, 1486-1532)
  - **Reassembly Validation**: Complete local reassembly for complex error patterns (lines 1383-1387, 1534-1620)

- **Error Detection Criteria** (lines 1451, 1499, 1576): Sophisticated multi-factor validation:
  - Coverage similarity analysis between flanking regions (isSimilar function calls)
  - Quality score integration for probabilistic error assessment
  - Junction detection to avoid correcting legitimate biological variation

- **Rollback Protection**: Comprehensive quality control preventing overcorrection:
  - Tracks correction count vs expected error rate ratios
  - Monitors k-mer count degradation after correction attempts
  - Implements complete rollback with original sequence restoration when corrections appear harmful

**Contig Processing and Graph Connectivity**:
- **Graph Initialization**: Establishes k-mer ownership mapping for assembled contigs with bidirectional k-mer claiming for both left and right termini

- **Edge Detection Algorithm**: Identifies connections between contigs through k-mer overlap analysis:
  - **Junction Analysis**: Uses fillLeftCounts/fillRightCounts to identify branching k-mers (lines 645-648, 682-685)
  - **Path Exploration**: exploreRight method traces connections up to 500bp with termination condition tracking
  - **Orientation Detection**: Determines forward/reverse complement relationships for scaffold construction

**Performance Optimization Features**:
- **Thread-Local K-mers**: Each BuildThread maintains private Kmer instances (myKmer, myKmer2) to eliminate synchronization overhead
- **Atomic Work Distribution**: Uses AtomicInteger arrays for lock-free table and victim processing coordination
- **Memory-Efficient Extension**: extendToRight2_inner implements bounded extension without ownership claiming for read extension operations
- **Fast Error Detection**: hasErrorsFast provides rapid error screening using sampling approach to avoid expensive full correction on clean reads

**Junk Detection and Quality Control**:
- **Read Quality Assessment**: Multi-metric evaluation including:
  - K-mer coverage depth analysis across entire read length
  - Left/right terminus connectivity validation (lines 1233-1234, 1253-1254)
  - Paired-read aware thresholds for improved accuracy (lines 1237, 1246)

- **Low-Quality K-mer Identification**: Configurable fraction-based filtering:
  - Counts valid vs invalid k-mers based on depth thresholds
  - Compensates for N-base effects on expected k-mer counts
  - Returns binary classification for downstream processing decisions

**Insert Size Estimation**: Specialized algorithm for paired-read insert size calculation:
- Identifies rightmost k-mers from both reads in pair
- Traces connection path between k-mers with maximum distance limits
- Accounts for k-mer overlap in final size calculation

**Usage**: Primary assembly engine for complex genomes requiring long k-mer resolution and sophisticated error correction. Typical workflow:
1. K-mer table construction via loadKmers() with coverage-based filtering
2. Iterative graph cleaning through shave() operations 
3. Multi-pass contig construction with progressive stringency
4. Error correction application with rollback protection
5. Graph connectivity analysis for scaffold construction

**Implementation Notes**: Extends Tadpole base class while replacing hash table implementation with KmerTableSetU for long k-mer support. The ownership system ensures thread safety during parallel processing, while the multi-algorithm error correction framework provides robust quality control for complex assembly scenarios. All k-mer operations use long arrays for k-mers exceeding 32bp, with automatic fallback to integer operations when possible for performance optimization.

### Rollback.java - Assembly Operation State Management and Recovery

**Purpose**: Provides transactional rollback functionality for reverting failed assembly operations by storing original read state before correction or extension attempts

**Core Function**: Creates immutable snapshots of Read objects with optional k-mer count data to enable recovery from unsuccessful assembly operations. Implements deep copying for critical sequence data while handling nullable quality scores and count lists appropriately.

**Key Features**:

- **Immutable Checkpoint Storage**: Final fields ensure checkpoint data cannot be modified after creation:
  - `id0`: Original read identifier string
  - `flags0`: Original processing flags integer  
  - `bases0`: Deep copy of sequence data preventing corruption
  - `quals0`: Deep copy of quality scores (null if unavailable)
  - `counts0`: Deep copy of k-mer frequency statistics (public for algorithm access)

- **Dual Constructor Interface** (lines 20-22, 29-35): Supports reads with or without k-mer count data:
  - Single-parameter constructor for basic read checkpointing
  - Two-parameter constructor including IntList k-mer counts for frequency-based algorithms
  - Automatic deep copying prevents reference-based data corruption

- **Intelligent Restoration Logic**: Optimizes restoration based on sequence length consistency:
  - **Same-length optimization**: Uses System.arraycopy() for efficient byte array restoration when length unchanged
  - **Length-mismatch handling**: Assigns reference copies when sequence length changed during processing
  - **Count list restoration**: Handles both array copying and collection-based restoration for k-mer statistics

- **Memory Safety Design**: Deep copying strategy prevents checkpoint corruption:
  - `bases.clone()` ensures sequence modifications don't affect saved state
  - `quality.clone()` protects quality score checkpoints when available  
  - `counts.copy()` preserves k-mer frequency data independently

- **Null Safety Handling** (lines 33-34, 55-56, 60-63): Comprehensive null checking for optional data:
  - Quality scores may be null in some sequencing formats
  - K-mer counts are algorithm-specific and may not exist
  - Restoration methods handle all null combinations safely

**Restoration Algorithm**:
- **Phase 1**: Restore core read properties (ID and flags) unconditionally
- **Phase 2**: Determine restoration strategy based on sequence length comparison
- **Phase 3a** (efficient path): Use System.arraycopy for same-length arrays
- **Phase 3b** (reference path): Assign saved references for length mismatches

**Usage**: Essential for assembly algorithms that attempt risky operations on reads:
- Error correction algorithms can revert unsuccessful fixes
- Extension algorithms can restore original state after failed elongation
- Quality improvement routines can rollback unsuccessful modifications
- Iterative assembly methods can reset to known good states

**Implementation Notes**: The length-based restoration strategy optimizes for common assembly scenarios where sequence length remains constant during processing, while gracefully handling expansion/contraction cases through reference assignment. Public access to `counts0` enables algorithm-specific k-mer frequency analysis during rollback decisions.

### Shaver1.java - Specialized Dead-End and Bubble Removal Algorithm

**Purpose**: Implements specialized algorithms for removal of dead ends (hairs) and bubble structures from De Bruijn graphs using bidirectional path exploration and termination code analysis

**Core Function**: Performs graph cleaning through exhaustive path exploration from k-mer seeds, categorizing termination reasons to identify removable structures. Uses atomic work distribution across multiple hash tables to enable parallel processing while maintaining thread safety through ownership claims (lines 58-113, 125-251).

**Key Features**:

- **Dual Structure Removal** (lines 32, 91-109): Configurable removal of two problematic graph structures:
  - **Dead-End Removal**: Eliminates linear paths terminating in single-branch dead ends (removeHair flag, lines 91-99)
  - **Bubble Removal**: Collapses alternative paths between shared branch points (removeBubbles flag, lines 101-109)

- **Bidirectional Path Exploration**: Core algorithm `exploreAndMark()` performs systematic graph traversal:
  - Forward exploration from initial k-mer seed using `explore()` method
  - Reverse-complement exploration for complete path characterization
  - Termination code analysis to classify path endpoints (lines 71-73, 76-84, 86-88)

- **Termination Classification System**: Seven distinct termination codes for path analysis:
  - `DEAD_END`: No valid extensions available (lines 146-149, 189-192, 288-290, 335-337)
  - `TOO_LONG`: Path exceeds maximum exploration distance (lines 248-250, 389-391)
  - `TOO_DEEP`: K-mer coverage exceeds threshold indicating repetitive regions (lines 239-242, 380-382)
  - `F_BRANCH`: Forward branching detected via multiple valid extensions (lines 234-237, 338-340)
  - `B_BRANCH`: Backward branching identified through count comparison analysis (lines 212-228, 360-376)
  - `LOOP`: Circular path detection via first k-mer comparison (lines 178-181, 311-314)

- **Multithreaded Processing Architecture**:
  - **ExploreThread**: Parallel dead-end detection using atomic table assignment (lines 418-435, 438-457)
  - **ShaveThread**: Concurrent k-mer removal operations with ownership cleanup
  - **High-Count Optimization**: Specialized processing for high-coverage k-mers via neighbor exploration

**Algorithm Implementation Details**:

- **K-mer Extension Logic**: Implements sophisticated path following algorithm:
  - Nucleotide count analysis using `fillRightCounts()` and `fillLeftCounts()` (lines 144, 187, 197)
  - Bit-level k-mer construction with mask operations (lines 132-134, 170-174)
  - Branch detection via second-highest count analysis (lines 155-156, 199-200, 212-228)

- **Coverage-Based Filtering** (lines 142, 158, 239-242): Uses configurable thresholds:
  - `minCount`: Minimum k-mer coverage for valid extensions (constructor line 32)
  - `maxCount`: Maximum coverage threshold to avoid repetitive regions (constructor line 32)
  - `minSeed`: Minimum coverage for exploration seeds (lines 495, 508, 528, 542)

- **Loop Detection Mechanism** (lines 138-139, 178-181, 268-269, 311-314): Maintains first k-mer reference for circular path identification

- **Thread Safety Implementation**:
  - **Atomic Work Distribution**: Uses `nextTable.getAndAdd(1)` for lock-free table assignment (lines 419, 439, 762)
  - **Ownership Claims**: Thread-safe k-mer marking via `claim()` operations (lines 77, 82, 94, 104, 111)
  - **Status Management**: Hierarchical status system (STATUS_UNEXPLORED, STATUS_EXPLORED, STATUS_REMOVE) prevents race conditions

**Performance Optimization Features**:

- **High-Count Neighbor Search**: Optimized algorithm for high-coverage k-mers:
  - Processes neighbors of high-count k-mers rather than the k-mers themselves
  - Uses bit masking for efficient neighbor generation
  - Implements safe fallback for core mask conflicts (lines 632, 708-737)

- **Memory-Efficient Processing**: Thread-local bit manipulation variables:
  - Pre-calculated shift values for k-mer size operations
  - Bit masks for efficient DNA encoding/decoding

**Usage**: Integrated into assembly pipelines for graph cleaning after initial k-mer table construction:
- Input: KmerTableSet with coverage statistics and configurable processing parameters
- Processing: Parallel exploration threads identify removable structures, followed by removal threads
- Output: Cleaned k-mer tables with dead ends and bubbles eliminated

**Implementation Notes**: The algorithm implements conservative defaults with low thresholds and both removal types enabled. The dual `explore()` and `explore2()` methods (lines 125-251, 256-392) suggest algorithm refinement, with `explore2()` providing streamlined termination logic. Thread-local data structures optimize performance by avoiding synchronization overhead during intensive graph traversal operations.

### Tadpole.java - Short K-mer De Bruijn Graph Assembler Framework

**Purpose**: Abstract base class implementing a comprehensive short k-mer assembler based on KmerCountExact with support for assembly, read extension, error correction, and sequence quality improvement.

**Core Function**: Provides unified framework for De Bruijn graph-based genome assembly using short k-mers (≤31bp by default) with automatic delegation to Tadpole2 for longer k-mers (>31bp, lines 72-85). Implements four primary processing modes: contig assembly, read extension, error correction, and read filtering/discarding (lines 230-246, 567-581).

**Architecture Overview**:

- **Polymorphic K-mer Implementation**: Factory pattern creates Tadpole1 for short k-mers (≤31bp) or Tadpole2 for long k-mers via makeTadpole() method. Uses AbstractKmerTableSet.MASK_CORE flag to optimize memory layout for different k-mer sizes.

- **Multi-Mode Processing Pipeline**: Four distinct operational modes configured via `mode` parameter:
  - **Contig Mode** (lines 234-235, 578): Default assembly mode for building contigs from k-mer graphs
  - **Extend Mode** (lines 236-237, 572): Read extension using graph traversal for gap filling
  - **Correct Mode** (lines 238-239, 568, 598): Error correction via k-mer frequency analysis  
  - **Discard Mode** (lines 242-243, 575): Quality filtering and problematic read removal

- **Comprehensive Parameter System**: Extensive command-line argument parsing supporting assembly parameters, error correction settings, quality thresholds, and I/O configuration

**Key Features**:

- **Adaptive K-mer Processing** (lines 88-104, 117): Dynamic k-mer size selection with automatic optimization:
  - K-mer size parsed from command line via preparseK() method
  - Automatic adjustment using Kmer.getMult() and Kmer.getK() for efficient representation
  - Supports range filtering with kmerRangeMin/kmerRangeMax thresholds (lines 292-295, 797-799)

- **Advanced Error Correction System**: Multi-algorithm error correction with configurable aggressiveness:
  - **Pincer Mode**: Bidirectional error correction requiring confirmation from both directions (lines 388-389, 394-395)
  - **Tail Mode**: Single-direction correction for read ends
  - **Reassemble Mode**: Complete sequence reconstruction for severe errors
  - **Conservative/Aggressive Modes**: Tunable parameter sets for different error rates

- **Graph Processing Operations** (lines 907-932, 1052-1102): Comprehensive graph analysis and improvement:
  - **Shave and Rinse**: Dead-end removal and bubble elimination via shaveAndRinse() method
  - **Bubble Popping**: Advanced bubble detection with configurable passes (lines 1061-1067, 201-204)
  - **Contig Graph Generation**: Edge detection and connectivity analysis

- **Read Extension Algorithm**: Sophisticated extension system for gap filling:
  - **Bidirectional Extension**: Separate left/right extension with configurable distances (lines 256-264, 1663-1670)
  - **Junction Handling**: Configurable behavior for branching points via extendThroughLeftJunctions (lines 323-324, 1451)
  - **Extension Rollback**: Quality-based trimming of low-confidence extensions (lines 262-263, 1671-1682)

**Threading Architecture**:

- **Parallel Processing Framework**: Multi-level parallelism with configurable thread counts (lines 301-302, 1014-1050):
  - **Build Threads**: Dedicated threads for contig construction (lines 197-198, 1017-1022)
  - **Extend Threads**: Specialized threads for read extension operations
  - **Process Threads**: Graph analysis and connectivity detection

- **Thread-Local Storage**: Comprehensive per-thread data structures to avoid synchronization overhead:
  - DNA base counting arrays, k-mer objects, sequence builders, and error tracking
  - Automatic initialization via initializeThreadLocals() method

**Error Correction Algorithms**:

- **Multi-Pass Reassembly**: Iterative error correction with bidirectional validation:
  - Builds candidate sequences from both directions
  - Applies sliding window quality filtering
  - Performs consensus calling with rollback protection

- **Quality-Based Filtering**: Sophisticated base marking system:
  - Marks low-coverage bases as unreliable using k-mer frequency thresholds
  - Implements delta-only marking to preserve high-confidence regions
  - Adjusts quality scores based on correction confidence

- **Merge Validation**: Statistical validation of read merging operations:
  - Analyzes k-mer frequency patterns around merge boundaries
  - Detects chimeric assemblies using coverage discontinuities
  - Configurable sensitivity via testMergeWidth and testMergeThresh parameters

**Memory Management**:

- **Efficient K-mer Storage**: Optimized memory layout with optional core masking
- **Buffer Management**: Sorted buffers and fast fill optimizations (lines 486-488, 333-334)
- **Thread Pool Control**: Configurable thread counts with automatic CPU detection

**Input/Output Processing**:

- **Flexible File Support**: Handles multiple input/output files with automatic format detection
- **Compression Integration**: Built-in support for compressed files via pigz/unpigz
- **Statistics Output**: Comprehensive assembly metrics with optional detailed reporting (lines 254-255, 749-875)

**Quality Control Features**:

- **Coverage Filtering**: Min/max coverage thresholds for k-mer inclusion (lines 275-279, 797-799)
- **Length Constraints**: Configurable minimum/maximum contig lengths
- **Junk Detection**: Automatic identification and removal of unassemblable reads (lines 373-383, 1649-1660)

**Usage**: Primary entry point for short k-mer genome assembly workflows:
- **Command Line Interface**: Comprehensive parameter system via main() method
- **Programmatic Access**: Factory creation via makeTadpole() for integration into larger pipelines
- **Mode Selection**: Automatic mode detection based on parameters or explicit mode specification

**Implementation Notes**: 
- Abstract class requiring concrete implementation of k-mer table operations, extension algorithms, and error correction methods
- Extensive parameter validation and automatic optimization based on input characteristics
- Built-in performance monitoring with detailed timing and statistics reporting
- Thread-safe design supporting parallel execution across all operational modes
- Supports both forward and reverse complement processing with canonical k-mer handling
