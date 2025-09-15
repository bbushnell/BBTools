# BBMin Package

This package provides minimizer generation functionality for genomic sequence analysis.

## Minimizer (Minimizer.java)
**Purpose**: Generates an array of minimal hash codes for a sequence's k-mer windows
**Core Function**: Generates canonical hash codes for overlapping k-mers, tracking minimal values
**Key Features**:
- Computes hash codes for forward and reverse complement k-mers
- Tracks minimal hash code within sliding windows
- Eliminates duplicate minimizers 
- Supports flexible k-mer and window size configurations
**Usage**: Used for genomic sequence fingerprinting and minimizer-based sequence mapping