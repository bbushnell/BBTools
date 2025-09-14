# BBTools Structures Package

## Overview

The `structures` package provides the foundational data structures and algorithms that power the entire BBTools suite. This package implements high-performance, memory-efficient collections and specialized data structures optimized for bioinformatics applications. These structures handle everything from basic integer lists to sophisticated hash maps, heaps, and specialized collections for genomic data processing.

## Core Design Principles

### Performance Optimization
The structures in this package are designed with performance as the primary concern:
- **Primitive Collections**: Uses primitive arrays instead of boxed objects to reduce memory overhead and eliminate boxing/unboxing costs
- **Memory Efficiency**: Minimal memory footprint through careful design and specialized implementations
- **Cache-Friendly Layouts**: Data structures optimized for modern CPU cache hierarchies
- **Reduced GC Pressure**: Fewer object allocations and better memory reuse

### Bioinformatics-Specific Features
Many structures include specialized functionality for genomic data:
- **Large-Scale Datasets**: Designed to handle millions or billions of genomic elements
- **Statistical Operations**: Built-in support for percentiles, histograms, and statistical analysis
- **Genomic Coordinates**: Specialized handling of genomic positions and ranges
- **k-mer Processing**: Optimized structures for k-mer counting and analysis

## Package Structure

### Core Collection Types

#### List Implementations
- **`IntList`**: High-performance integer list with specialized operations for sorted data
- **`IntList2`**: Enhanced integer list with additional features
- **`IntList3`**: Further optimized integer list implementation
- **`LongList`**: Long integer list for large values and genomic coordinates
- **`LongList2`**: Enhanced long list with extended functionality
- **`LongList3`**: Advanced long list implementation
- **`DoubleList`**: Double-precision floating point list for statistical data
- **`FloatList`**: Single-precision floating point list for performance

#### Hash Map Implementations
- **`IntHashMap`**: Primitive integer-to-integer mapping with high performance
- **`LongHashMap`**: Long-to-integer mapping optimized for genomic coordinates
- **`LongLongHashMap`**: Long-to-long mapping for large value ranges
- **`LongLongHashMapHybrid`**: Hybrid implementation optimizing for different value distributions
- **`IntLongHashMap`**: Integer-to-long mapping for mixed data types
- **`LongArrayListHashMap`**: Maps long keys to lists of values

#### Hash Set Implementations
- **`IntHashSet`**: High-performance integer set for membership testing
- **`LongHashSet`**: Long integer set for large coordinate ranges
- **`LongListSet`**: Set implementation backed by sorted lists for memory efficiency

#### Specialized Hash Maps
- **`IntHashSetList`**: Maps integers to sets of integers
- **`IntListHashMap`**: Maps integers to lists of integers
- **`LongLongListHashMap`**: Maps longs to lists of longs

### Advanced Data Structures

#### Heap Implementations
- **`Heap`**: Generic min-heap with rollover capability for top-N problems
- **`HeapLoc`**: Location-aware heap for spatial algorithms
- **`LongHeap`**: Primitive long heap for performance-critical applications
- **`LongHeapMap`**: Heap-based map for maintaining top-N key-value pairs
- **`LongHeapSet`**: Heap-based set for top-N element selection

#### Bit Set Implementations
- **`AbstractBitSet`**: Base class for bit manipulation structures
- **`RawBitSet`**: Low-level bit manipulation with minimal overhead
- **`AtomicBitSet`**: Thread-safe bit set for concurrent access
- **`MultiBitSet`**: Multi-bit values per position for advanced encoding

#### Specialized Collections
- **`SuperLongList`**: Hybrid structure for histogram data with array + list storage
- **`ByteBuilder`**: High-performance string builder for byte data
- **`RingBuffer`**: Circular buffer implementations for streaming data
- **`RingBufferMask`**: Mask-based ring buffer for power-of-2 sizes
- **`RingBufferMod`**: Modulo-based ring buffer for arbitrary sizes

### Genomic-Specific Structures

#### Coverage Analysis
- **`CoverageArray`**: Basic coverage tracking for genomic positions
- **`CoverageArray2`**: Enhanced coverage with additional statistics
- **`CoverageArray2A`**: Array-optimized coverage implementation
- **`CoverageArray3`**: Advanced coverage with multi-level features
- **`CoverageArray3A`**: Array-optimized advanced coverage

#### Coordinate Handling
- **`Range`**: Genomic coordinate range representation
- **`CRange`**: Compressed range for memory efficiency
- **`Point`**: 2D coordinate representation
- **`Feature`**: Genomic feature with metadata
- **`SeqPos`**: Sequence position with orientation
- **`SeqPosM`**: Multi-sequence position tracking

#### Counting and Statistics
- **`SeqCount`**: Sequence occurrence counting
- **`SeqCountM`**: Multi-sequence count tracking
- **`StringCount`**: String occurrence counting with statistics
- **`Quantizer`**: Value quantization for histogram generation
- **`StandardDeviator`**: Efficient standard deviation calculation

### Utility Structures

#### Text and String Handling
- **`ByteBuilder`**: High-performance byte array builder with extensive formatting support
- **`StringNum`**: String-to-number associations
- **`AtomicStringNum`**: Thread-safe string-number mapping
- **`StringPair`**: Paired string storage

#### Numeric Pairs and Tuples
- **`LongPair`**: Paired long values for coordinates
- **`LongM`**: Extended long value with metadata

#### List Utilities
- **`ListNum`**: List with associated numeric identifier
- **`IntListCompressor`**: Compression utilities for integer lists

## Key Features and Usage Patterns

### High-Performance Collections

The primitive collections provide significant performance advantages over standard Java collections:

```java
// Memory-efficient integer storage
IntList list = new IntList();
for(int i = 0; i < 1000000; i++) {
    list.add(i);  // No boxing overhead
}

// Fast hash-based counting
LongHashMap counter = new LongHashMap();
counter.increment(key);  // Efficient k-mer counting
```

### Statistical Analysis

Many structures include built-in statistical operations:

```java
IntList data = new IntList();
// ... populate with data
data.sort();
double median = data.percentile(0.5);
double mean = data.sum() / data.size();
```

### Heap-Based Operations

Efficient top-N selection and maintenance:

```java
Heap<Integer> topValues = new Heap<>(100, true);  // Keep top 100
for(int value : dataset) {
    topValues.add(value);  // Automatically maintains top N
}
```

### Memory-Efficient Building

ByteBuilder provides high-performance string construction:

```java
ByteBuilder sb = new ByteBuilder();
sb.append("Value: ").append(123).append(" (").append(45.67, 2).append(")");
String result = sb.toString();
```

## Advanced Features

### Coverage Analysis

Sophisticated genomic coverage tracking:

```java
CoverageArray3 coverage = new CoverageArray3(chromosomeLength);
coverage.incrementRange(start, stop, increment);
double avgCoverage = coverage.averageCoverage();
```

### Hybrid Storage

SuperLongList efficiently handles histogram data:

```java
SuperLongList histogram = new SuperLongList();
// Small values stored in array, large values in list
histogram.add(smallValue);  // Fast array access
histogram.add(largeValue);  // List storage for outliers
```

### Multi-Threading Support

Thread-safe structures for concurrent access:

```java
AtomicBitSet threadSafeBits = new AtomicBitSet(size);
// Safe for concurrent access from multiple threads
```

## Performance Characteristics

### Memory Efficiency
- **IntList vs ArrayList<Integer>**: ~75% less memory usage
- **Primitive collections**: 4-16x less memory per element
- **Specialized structures**: Optimized for specific use patterns

### Speed Optimizations
- **Direct array access**: Eliminates indirection overhead
- **Cache-friendly layouts**: Better CPU cache utilization
- **Reduced allocations**: Lower garbage collection pressure
- **Specialized algorithms**: Optimized for biological data patterns

### Scalability Features
- **Large dataset support**: Handles billions of elements efficiently
- **Memory-aware sizing**: Automatic capacity management
- **Streaming-friendly**: Support for incremental processing

## Integration with BBTools Ecosystem

### Core Infrastructure
These structures form the foundation for:
- **Alignment algorithms**: Coordinate tracking and scoring
- **k-mer processing**: High-speed counting and storage
- **Coverage analysis**: Genomic position coverage tracking
- **Statistical analysis**: Quality metrics and data summarization

### Memory Management
- **KillSwitch integration**: Centralized memory allocation tracking
- **Shared memory pools**: Efficient memory reuse across tools
- **Large array support**: Handles datasets exceeding standard limits

### Performance Monitoring
- **Built-in benchmarks**: Performance comparison utilities
- **Memory usage tracking**: Real-time memory consumption monitoring
- **Statistics collection**: Built-in performance metrics

## Best Practices

### Choosing the Right Structure
1. **Use primitive collections** when possible for better performance
2. **Consider memory patterns** when selecting between array and list-based structures
3. **Leverage specialized structures** for domain-specific operations
4. **Use heap structures** for top-N problems and priority queues

### Performance Optimization
1. **Size appropriately**: Initialize with expected capacity when known
2. **Use batch operations** for bulk data processing
3. **Avoid boxing/unboxing** by using primitive-specific methods
4. **Consider memory locality** when accessing data sequentially

### Memory Management
1. **Monitor memory usage** for large datasets
2. **Use shrink methods** to reduce memory footprint after processing
3. **Clear collections** when no longer needed
4. **Consider hybrid structures** for mixed data distributions

## Design Patterns

### Builder Pattern
ByteBuilder demonstrates efficient incremental construction with method chaining.

### Factory Pattern
Specialized constructors and static methods provide optimized instances for different use cases.

### Template Method Pattern
Abstract base classes define common interfaces while allowing specialized implementations.

### Strategy Pattern
Multiple implementations of similar structures allow selection based on performance requirements.

This structures package provides the high-performance foundation that enables BBTools to efficiently process massive genomic datasets while maintaining memory efficiency and computational speed. The careful optimization of these fundamental data structures directly impacts the performance of all higher-level BBTools algorithms.