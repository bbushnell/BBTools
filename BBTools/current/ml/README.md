# ML Package - Neural Network and Machine Learning Infrastructure

BBTools machine learning package providing comprehensive neural network implementation, activation functions, training utilities, and bioinformatics-specific sequence processing capabilities.

## Neural Network Core Components

### CellNet (CellNet.java)
**Purpose**: Multi-layer feedforward neural network implementation supporting dense and sparse network architectures with configurable training methodologies.
- **Core Function**: Provides complete neural network training, forward propagation, backpropagation, and gradient descent learning with support for various connectivity patterns and optimization techniques
- **Key Features**:
  - Supports both dense and sparse network connectivity patterns
  - Configurable layer dimensions, connection densities, and block-aligned sparsity
  - Thread-safe training with parallel gradient accumulation
  - Advanced weight initialization with multi-distribution sampling
  - Built-in simulated annealing and weight normalization
  - Comprehensive performance metrics tracking (error rates, FPR, FNR)
  - Serialization and network state reproduction support
- **Usage**: Machine learning model training across various domains, flexible neural network architecture exploration, supporting research and experimental neural network implementations

### Cell (Cell.java)
**Purpose**: Represents a single neural network cell with advanced computational and training capabilities.
- **Core Function**: Handles neural network cell operations including activation, error computation, edge updates, and weight management for both dense and sparse network architectures.
- **Key Features**:
  - Supports multiple activation functions through dynamic function selection
  - Handles dense and sparse network connectivity modes
  - Implements vectorized weight and gradient computations
  - Supports advanced error calculation with class weighting
  - Provides methods for edge updates in final and hidden layers
- **Usage**: Core computational unit in neural network training, manages individual neuron state, including weights, biases, and activation values

## Activation Functions

### Function (Function.java)
**Purpose**: Abstract base class for defining activation functions in neural network architectures.
- **Core Function**: Provides a standardized interface for implementing various activation function transformations and their derivatives in machine learning models.
- **Key Features**:
  - Abstract method `activate(double x)` for computing activation function output
  - Abstract method `derivativeX(double x)` for calculating function derivatives during backpropagation
  - Static method for converting function types from string representations
  - Supports multiple built-in activation function types (Sigmoid, Tanh, RsLog, etc.)
  - Weighted random function selection via `randomFunction()`
- **Usage**: Serves as a template for defining different activation functions in neural network implementations

### Sigmoid (Sigmoid.java)
**Purpose**: Sigmoid activation function for neural networks, implementing the logistic sigmoid transformation which outputs values between 0 and 1.
- **Core Function**: Provides a standard sigmoid activation implementation, converting input values to a probability-like output range used in classification tasks
- **Key Features**:
  - Implements `Function` abstract class for activation function standardization
  - Delegates activation calculations to `Functions.sigmoid(x)`
  - Provides derivative methods for backpropagation
  - Static singleton instance with predefined name "SIG"
- **Usage**: Used in neural network architectures requiring a standard sigmoid activation function, commonly applied in binary classification problems

### Tanh (Tanh.java)
**Purpose**: Hyperbolic tangent activation function for neural networks with symmetric output range between -1 and 1.
- **Core Function**: Implements the tanh activation function, providing a symmetrically centered activation transformation for neural network layers
- **Key Features**:
  - Extends abstract `Function` class
  - Provides activation via `Functions.tanh(x)` method
  - Symmetric output range from -1 to 1
  - Offers multiple derivative calculation methods for backpropagation
- **Usage**: Used in neural network architectures requiring a symmetrically centered activation function

### Swish (Swish.java)
**Purpose**: Swish activation function implementation for neural network architectures.
- **Core Function**: Implements the swish(x) = x * sigmoid(x) activation function, providing smoother gradients and improved performance in deep neural networks.
- **Key Features**:
  - Extends the abstract `Function` class for standardized activation function implementation
  - Uses `Functions.swish(x)` for activation computation
  - Provides multiple derivative calculation methods for backpropagation
  - Static singleton instance with pre-defined name and type
- **Usage**: Used in neural network layers requiring the Swish activation function, offers an alternative to traditional activation functions like ReLU

### MSig (MSig.java)
**Purpose**: Mirrored sigmoid activation function implementation for neural networks.
- **Core Function**: Provides a symmetric sigmoid activation function centered around zero for specialized neural network architectures
- **Key Features**:
  - Implements the `Function` abstract class
  - Delegates activation and derivative calculations to `Functions` utility class
  - Static singleton implementation with consistent access
  - Uses `Functions.mSig(x)` for symmetric activation
- **Usage**: Used in neural network layers requiring a symmetric activation function with specialized sigmoid behavior

### ExtendedSigmoid (ExtendedSigmoid.java)
**Purpose**: Provides an extended sigmoid activation function implementation in the machine learning function hierarchy.
- **Core Function**: Implements derivative and activation calculations for an extended sigmoid function using utility methods from Functions class
- **Key Features**:
  - Extends abstract Function class
  - Singleton pattern with static instance
  - Supports multiple derivative calculations
- **Usage**: Used in neural network activation calculations, specifically for extended sigmoid transformations

### ExtendedMSig (ExtendedMSig.java)
**Purpose**: Extended mirrored sigmoid activation function for neural networks with symmetric [-1, 1] output range.
- **Core Function**: Provides a specialized activation function implementing a symmetric, extended-range sigmoid for neural network layers
- **Key Features**:
  - Implements `Function` abstract class
  - Generates activation value via `Functions.emSig(x)`
  - Calculates derivative using `Functions.emSigDerivativeX(x)`
  - Static singleton instance for consistent access
- **Usage**: Used in neural network architectures requiring a symmetric activation function with enhanced expressiveness between -1 and 1

### Bell (Bell.java)
**Purpose**: Bell curve activation function for neural networks, generating values between 0 and 1.
- **Core Function**: Implements a Gaussian-like activation function with maximum output at x=0
- **Key Features**:
  - Overrides activation method using `Functions.bell(x)`
  - Provides derivatives for x, function value (fx), and combined derivative
  - Singleton implementation with static `instance` field
- **Usage**: Used in neural network architectures requiring bell-shaped activation responses, particularly where output needs to be normalized between 0 and 1

### RSLog (RSLog.java)
**Purpose**: Rotational Symmetric Logarithm (RSLOG) activation function implementation for neural networks.
- **Core Function**: Provides a specialized activation function from the Function hierarchy, utilizing RSLOG mathematical transformation for neural network layers.
- **Key Features**:
  - Extends abstract `Function` class for activation function standardization
  - Implements `activate(double x)` method using `Functions.rslog(x)`
  - Provides derivative calculations for backpropagation
  - Static singleton instance with pre-defined name and type
- **Usage**: Used in neural network architectures requiring the Rotational Symmetric Logarithm activation function

### Functions (Functions.java)
**Purpose**: Utility class providing mathematical activation functions and derivative calculations for machine learning models.
- **Core Function**: Implements a collection of static mathematical transformation methods commonly used in neural network and machine learning architectures, including sigmoid, tanh, swish, and specialized derivatives.
- **Key Features**:
  - Standard and Extended Sigmoid implementations (`sigmoid()`, `eSigmoid()`)
  - TanH activation function (`tanh()`)
  - Swish activation function (`swish()`)
  - Specialized derivative calculations for each activation function
  - Gaussian bell curve function (`bell()`)
  - Mirrored Sigmoid variants (`mSig()`, `emSig()`)
  - Robust Safe Log function (`rslog()`)
  - Mean Squared Error calculation for model evaluation (`mse()`)
- **Usage**: Provides mathematical transformations for neural network activation layers, supports gradient calculations during backpropagation

## Training Infrastructure

### Trainer (Trainer.java)
**Purpose**: Comprehensive neural network training orchestrator that manages complex machine learning model generation, configuration, and evaluation processes.
- **Core Function**: Coordinates end-to-end neural network training lifecycle, including network initialization, parallel processing, performance tracking, and output generation across multiple training cycles.
- **Key Features**:
  - Supports multi-threaded neural network training with configurable worker threads
  - Handles network generation via random seeding or pre-configured network loading
  - Implements advanced training techniques like simulated annealing and learning rate scheduling
  - Provides comprehensive performance metrics tracking (error rates, FPR, FNR, TPR, TNR)
  - Manages data loading, preprocessing, and subset creation
  - Supports flexible network architecture exploration
- **Usage**: Primary driver for neural network training in machine learning research and applications, enables systematic exploration of neural network architectures

### TrainerThread (TrainerThread.java)
**Purpose**: Manages parallel neural network processing tasks, enabling distributed computation across multiple threads with fine-grained job processing.
- **Core Function**: Orchestrates thread-safe sample processing, network computation, and result collection during distributed machine learning training
- **Key Features**:
  - Implements `Comparable<TrainerThread>` for thread identification and sorting
  - Thread-safe job processing using `ArrayBlockingQueue` for job distribution
  - Supports dynamic job configuration with job metadata tracking
  - Implements built-in performance profiling with nanosecond-level precision
- **Usage**: Parallel neural network training across multiple computational threads, distributed processing of machine learning sample sets

### WorkerThread (WorkerThread.java)
**Purpose**: Manages parallel neural network processing tasks, enabling distributed computation across multiple threads with fine-grained job processing.
- **Core Function**: Orchestrates thread-safe sample processing, network computation, and result collection during distributed machine learning training
- **Key Features**:
  - Implements `Comparable<WorkerThread>` for thread identification and sorting
  - Thread-safe job processing using `ArrayBlockingQueue` for job distribution
  - Supports dynamic job configuration with job metadata tracking
  - Implements built-in performance profiling with nanosecond-level precision
  - Handles sample processing, sorting, and synchronization
  - Accumulates comprehensive performance metrics (TP, TN, FP, FN rates)
- **Usage**: Parallel neural network training across multiple computational threads, implementing thread-safe gradient accumulation strategies

### ScannerThread (ScannerThread.java)
**Purpose**: Parallel neural network architecture generation and evaluation thread that explores candidate network configurations and returns top-performing networks.
- **Core Function**: Generates multiple neural network architectures, trains them, and selectively returns the most promising network configurations for further exploration.
- **Key Features**:
  - Implements multi-threaded network architecture search
  - Uses deterministic random seed generation for reproducible experiments
  - Maintains a priority queue (min-heap) of top-performing network candidates
  - Supports configurable training epochs and sample evaluation limits
  - Handles thread-safe job processing with synchronization mechanisms
- **Usage**: Automated neural network architecture exploration, parallel evaluation of multiple network configurations

## Data Management

### DataLoader (DataLoader.java)
**Purpose**: Loads and preprocesses machine learning datasets from files with sophisticated parsing and transformation capabilities.
- **Core Function**: Handles loading, parsing, shuffling, and splitting of machine learning datasets with support for multiple file formats
- **Key Features**:
  - Supports loading data from single or multiple file paths
  - Performs automatic data dimension inference
  - Handles data shuffling with a reproducible random seed
  - Implements class balancing by cloning underrepresented samples
  - Supports train/test dataset splitting with configurable parameters
- **Usage**: Prepares raw data files for machine learning model training, handling data preprocessing steps like shuffling, balancing, and splitting

### Matrix (Matrix.java)
**Purpose**: Represents a data matrix for machine learning operations, handling input/output data transformation, range detection, and format conversion.
- **Core Function**: Manages data matrix preprocessing, including range detection, binary classification conversion, and linear scaling of output values.
- **Key Features**:
  - Automatic range detection for output data
  - Binary classification conversion with configurable threshold
  - Linear range adjustment for output values
  - Statistical tracking of output data properties (min, max, mean, midpoint)
- **Usage**: Preprocessing step in machine learning pipelines, normalizing and transforming output data for neural network training

### Sample (Sample.java)
**Purpose**: Represents a single training sample in a machine learning dataset with comprehensive error computation and classification capabilities.
- **Core Function**: Manages individual data points for neural network training, including input features, target goals, prediction results, and sophisticated error calculation mechanisms
- **Key Features**:
  - Supports binary classification with a 0.5 threshold
  - Implements custom pivot calculation for sample prioritization
  - Computes error magnitude with optional class weighting
  - Thread-safe methods for pivot and epoch management
  - Supports conversion of sample to string and byte representations
- **Usage**: Represents an individual data point in machine learning training datasets, enables precise tracking of prediction errors and sample characteristics

### SampleSet (SampleSet.java)
**Purpose**: Manages a set of machine learning samples for training and evaluation across various performance metrics and data manipulation techniques.
- **Core Function**: Provides comprehensive functionality for generating, manipulating, and analyzing machine learning sample sets, including subset creation, ROC curve generation, and performance metric calculations.
- **Key Features**:
  - Generates samples from input matrices with automatic positive/negative categorization
  - Supports cross-validation through subset generation and rotation
  - Calculates key performance metrics like False Positive Rate (FPR) and False Negative Rate (FNR)
  - Generates Receiver Operating Characteristic (ROC) curve
  - Implements flexible sample shuffling with reproducible random seed
- **Usage**: Preprocessing and preparing machine learning datasets for training, performance evaluation of classification models, cross-validation setup for machine learning algorithms

### Subset (Subset.java)
**Purpose**: Manages sample subset creation, sorting, and manipulation for machine learning training datasets with advanced sample distribution techniques.
- **Core Function**: Provides methods for organizing, sorting, and shuffling training samples with support for binary classification, error-based prioritization, and reproducible randomization.
- **Key Features**:
  - Supports binary classification by splitting samples into positive and negative sets
  - Implements advanced sorting methods for sample prioritization based on error magnitude
  - Provides methods for multithreaded and single-threaded sorting of samples
  - Supports interleaving of samples from positive and negative classes for balanced training
  - Implements sample triage mechanism with configurable epochs and delay distances
- **Usage**: Preprocessing machine learning datasets with fine-grained sample management, implementing advanced sample selection strategies during neural network training

## Utility Components

### Seed (Seed.java)
**Purpose**: Represents a deterministic seed configuration for neural network training with performance tracking.
- **Core Function**: Manages network seeds and performance pivot for reproducible machine learning experiments, providing a comparable representation of network initialization states.
- **Key Features**:
  - Implements `Comparable<Seed>` interface for seed sorting and comparison
  - Stores network seed and performance pivot
  - Provides custom compareTo method for deterministic sorting
  - Supports equality checks based on network seed
- **Usage**: Used in neural network initialization to ensure reproducible random states, facilitates seed tracking and comparison during model training

### Status (Status.java)
**Purpose**: Represents a tracking object for neural network training state, enabling comparison and hashcode generation for network instances.
- **Core Function**: Encapsulates key training parameters including network reference, epoch counter, learning rate, and annealing factor, providing mechanisms for object comparison and identification.
- **Key Features**:
  - Immutable storage of training hyperparameters (epoch, alpha, anneal)
  - Implements `Comparable<Status>` for network instance sorting
  - Provides `hashCode()` generation based on epoch
  - Supports null-safe comparison in `compareTo()` method
- **Usage**: Tracking and comparing neural network training states, enabling sorting of network instances during performance evaluation

### Source (Source.java)
**Purpose**: Abstract base class defining fundamental properties for machine learning data sources with minimal implementation.
- **Core Function**: Provides a lightweight template for data sources with basic value, ID, and bias tracking capabilities
- **Key Features**:
  - Abstract `check()` method for source validation
  - Abstract `terminal()` method identifying sources with no inputs
  - Protected `value` field for storing computational results
  - Default methods for `id()` and `bias()` returning zero
- **Usage**: Serves as a base class for creating specialized data source implementations in machine learning contexts

### Profiler (Profiler.java)
**Purpose**: Lightweight, configurable performance timing utility for measuring execution phases in precise nanosecond increments.
- **Core Function**: Tracks elapsed time across multiple computational phases, accumulating timing data with minimal performance overhead
- **Key Features**:
  - High-precision nanosecond-level timing using `System.nanoTime()`
  - Configurable profiling with boolean `PROFILING` flag to enable/disable
  - Supports multiple phase timing accumulation in a single array
  - Allows resetting, clearing, and accumulating timing data
- **Usage**: Performance benchmarking of computational methods, tracking time spent in different phases of machine learning algorithms

## Job Management

### JobData (JobData.java)
**Purpose**: Encapsulates configuration and state for parallel machine learning job processing with thread-safe data management.
- **Core Function**: Provides a thread-safe container for neural network training jobs, managing sample subsets, network references, and job-specific parameters during distributed learning
- **Key Features**:
  - Thread-safe job configuration with immutable and optional mutable network references
  - Supports configurable learning rate, backpropagation, and weight regularization
  - Manages job-specific sample subsets within a larger dataset
  - Implements thread synchronization using `ReentrantReadWriteLock`
  - Provides a sentinel `POISON` job for safe thread termination
- **Usage**: Coordinates distributed neural network training across multiple threads, manages job-specific configurations like epoch, learning rate, and sample processing limits

### JobResults (JobResults.java)
**Purpose**: Encapsulates computational results and performance metrics for individual machine learning job executions during neural network training.
- **Core Function**: Tracks and stores detailed performance statistics for a specific training job, including error metrics, classification outcomes, and network reference
- **Key Features**:
  - Stores comprehensive job execution details including epoch, thread ID, and job sequence number
  - Tracks raw and weighted error sums
  - Captures classification performance metrics (true/false positives/negatives)
  - Provides a compact string representation for debugging
  - Implements `Comparable` for sorting job results by epoch and job ID
- **Usage**: Used in parallel neural network training to aggregate and compare job-level performance metrics, enables tracking of per-thread and per-job training statistics

## Comparators

### SampleValueComparator (SampleValueComparator.java)
**Purpose**: Comparator for sorting samples by their neural network output values.
- **Core Function**: Implements a three-way comparison mechanism for Sample objects, sorting based on their primary result value with an ID-based tie-breaker
- **Key Features**:
  - Implements Java's Comparator interface for Sample objects
  - Uses the first result value of a Sample for primary comparison
  - Provides deterministic ordering with ID-based secondary comparison
  - Singleton static instance for consistent access
- **Usage**: Sorting neural network samples by prediction confidence or classification score, ranking samples in machine learning evaluation and selection processes

### SampleErrorComparator (SampleErrorComparator.java)
**Purpose**: Deprecated comparator for comparing Sample objects with experimental error handling.
- **Core Function**: Implements a problematic comparator for Sample objects that always triggers an assertion failure
- **Key Features**:
  - Implements Java's Comparator interface
  - Marked as @Deprecated, indicating non-functional implementation
  - Contains a non-functional compare method with commented-out error value logic
  - Falls back to default Sample comparison via `compareTo()`
- **Usage**: Currently non-functional due to intentional assertion failure, appears to be an experimental or placeholder implementation

### CellNetReverseComparator (CellNetReverseComparator.java)
**Purpose**: Reverse comparator for sorting neural networks in descending order during training.
- **Core Function**: Implements a reverse comparison mechanism for CellNet objects
- **Key Features**:
  - Implements Java's Comparator interface
  - Reverses default comparison order 
  - Used for selecting best-performing networks
- **Usage**: Allows sorting of CellNet objects from highest to lowest performance metric

## Configuration and Parsing

### CellNetParser (CellNetParser.java)
**Purpose**: Parses and loads Cell Network configuration files with support for dense and sparse network representations.
- **Core Function**: Reads network configuration files, parsing header information and edge details for neural network structures
- **Key Features**:
  - Supports parsing both dense and sparse network configurations
  - Handles multiple header metadata including layers, seed, density, epochs
  - Extracts cell functions, biases, and weights from configuration
  - Supports loading from files with optional null-on-failure behavior
- **Usage**: Used to load pre-configured neural network structures from custom file formats, supports loading network configurations for machine learning model initialization

## Bioinformatics Integration

### NetFilter (NetFilter.java)
**Purpose**: Filters genomic sequences using neural network-based scoring and classification.
- **Core Function**: Processes input sequence files through a pre-trained neural network, scoring and optionally filtering sequences based on network predictions.
- **Key Features**:
  - Supports both single-end and paired-end sequence filtering
  - Configurable neural network scoring modes (single sequence, sliding window)
  - Flexible filtering options with adjustable cutoff thresholds
  - Supports forward and reverse complement sequence scoring
  - Can annotate sequences with neural network scores
- **Usage**: Sequence quality control in bioinformatics pipelines, filtering genomic reads based on machine learning predictions

### ScoreSequence (ScoreSequence.java)
**Purpose**: Neural network-based genomic sequence scoring and filtering utility for bioinformatics data processing.
- **Core Function**: Processes genomic sequences through a pre-trained neural network, computing scores, applying optional filtering, and generating comprehensive processing statistics.
- **Key Features**:
  - Supports neural network scoring for single-end and paired-end sequences
  - Configurable scoring modes including forward and reverse complement orientation
  - Dynamic neural network input width detection
  - Supports sequence annotation with network scores
  - Generates score distribution histograms
- **Usage**: Genomic sequence quality control in bioinformatics pipelines, neural network-based sequence classification and filtering

### SequenceToVector (SequenceToVector.java)
**Purpose**: Converts DNA sequences into numerical vectors for machine learning processing.
- **Core Function**: Transforms nucleotide sequences into mathematical representations suitable for neural network input
- **Key Features**:
  - Converts DNA sequences to numerical feature vectors
  - Supports various encoding strategies for nucleotide representation
  - Handles sequence preprocessing and normalization
- **Usage**: Preprocessing genomic data for machine learning models, converting biological sequences into numerical formats

### ProcessBBMergeHeaders (ProcessBBMergeHeaders.java)
**Purpose**: Processes and extracts statistical headers from BBMerge merged read files, transforming read headers into feature vectors for downstream analysis.
- **Core Function**: Parses BBMerge read merge statistics, extracting detailed insert size, overlap, error, and quality information from read headers
- **Key Features**:
  - Converts BBMerge read headers into structured statistical feature vectors
  - Supports parsing complex header information with multiple statistical markers
  - Handles concurrent read input and output streams for efficient processing
  - Validates and filters read headers based on merge quality criteria
- **Usage**: Post-processing BBMerge merged read files to extract comprehensive merge statistics, preparing read merge data for further bioinformatics analysis

## Data Processing Utilities

### ReduceColumns (ReduceColumns.java)
**Purpose**: Command-line utility for extracting specific columns from tab-delimited text files with efficient memory management.
- **Core Function**: Processes input files by selectively extracting user-specified columns, writing results to a new output file with low memory overhead
- **Key Features**:
  - Supports dynamic column selection via command-line arguments
  - Uses streaming file processing for large file handling
  - Preserves file formatting with tab-separated output
  - Tracks processing performance with timing and statistics
  - Skips comment lines (lines starting with '#')
- **Usage**: Column extraction from large tabular data files, data preprocessing for machine learning datasets, reducing dimensionality of input files