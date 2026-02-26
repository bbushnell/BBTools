# BBTools

**Official BBTools repository maintained by Brian Bushnell**

**BBTools** is a suite of fast, multithreaded bioinformatics tools designed for analysis of DNA and RNA sequence data. BBTools can handle common sequencing and bioinformatics file formats such as fastq, fasta, sam, bam, scarf, fasta+qual, gff, gtf, vcf, and gzip/bgzip/bzip2 compressed files.

## 👥 Authors

- **Brian Bushnell** - Primary Developer
- **Jon Rood** - Contributor
- **Shijie Yao** - Contributor
- **Isla** - Contributor

## 📊 Version

Current Version: **39.77**

## 🚀 Features

- **High Performance**: Multithreaded, efficient memory usage, and optimized algorithms
- **Comprehensive Suite**: Over 100 tools for various bioinformatics tasks
- **Format Flexibility**: Handles multiple file formats with automatic format detection and proper support for twin or interleaved files
- **Platform Independent**: Pure Java implementation runs on any system with Java 8+
- **Production Ready**: Used in production at JGI, and cited in thousands of publications

## 📦 Main Components

### Some Core Tools
- **BBMap**: Short read aligner for DNA and RNA-seq data
- **BBDuk**: Adapter trimming, quality filtering, and contaminant removal  
- **BBMerge**: Paired read merging with error correction
- **BBNorm**: Error correction and normalization
- **Tadpole**: Fast assembler for small genomes or metagenomes
- **Clumpify**: Reorder reads to maximize compression and speed up analysis
- **CallVariants**: Variant calling from aligned reads
- **BBCMS**: Error correction via conditional median filtering

### Quality Control
- **BBDuk**: Quality trimming and adapter removal
- **BBSplit**: Binning reads by mapping to multiple references
- **Clumpify**: Reorder reads to maximize compression

### Assembly Tools
- **Tadpole**: Extremely fast micro-assembler
- **BBTools Assembly Pipeline**: Complete assembly workflow

### Sketch Tools
- **SendSketch**: Identify organisms rapidly using a remote server
- **CompareSketch**: Compare organisms or custom databases locally using MinHash

## 🔧 Installation

### Quick Start
```bash
# Download the latest version
wget https://sourceforge.net/projects/bbmap/files/latest/download -O BBTools.tar.gz
tar -xzf BBTools.tar.gz
cd bbmap

# Test installation
./bbversion.sh
```

### Requirements
- Java 8 or higher (Java 25+ recommended)
- Bash recommended but not required

## 📖 Documentation

Comprehensive documentation is available in the `/docs` directory:
- `readme.txt` - General information
- `UsageGuide.txt` - Detailed usage instructions
- `ToolDescriptions.txt` - Description of all tools

## 💻 Usage Examples

### Adapter Trimming
```bash
bbduk.sh in=reads.fq out=clean.fq ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo
```

### Read Mapping
```bash
bbmap.sh ref=reference.fasta in=reads.fq out=mapped.sam
```

### Error Correction
```bash
tadpole.sh in=reads.fq out=corrected.fq mode=correct k=62
```

## 🤝 Contributing

We welcome contributions! Contact Brian if there is a tool or functionality you want to add.

## 📄 License

BBTools is free for unlimited use.

## 📚 Citation

If you use BBTools in your research, please cite:
```
Bushnell, B. (2014) BBMap: A Fast, Accurate, Splice-Aware Aligner.

BBTools official website: https://bbmap.org
GitHub: https://github.com/bbushnell/BBTools
```


## 🔗 Links

- **Official Website**: [bbmap.org](https://bbmap.org)
- **SourceForge**: [sourceforge.net/projects/bbmap](https://sourceforge.net/projects/bbmap/)

---

*BBTools is developed at the Joint Genome Institute (JGI), part of Lawrence Berkeley National Laboratory (LBNL)*
