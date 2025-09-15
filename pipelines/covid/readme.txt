Written by Brian Bushnell
Last modified May 8, 2020

Contents:
This directory contains a collection of scripts written for calling variants and generating consensus genomes from SARS-CoV-2 ("Covid") Illumina data using BBTools.  They were designed and optimized for several sets of libraries (from NextSeq and NovaSeq platforms) containing a mix of shotgun and amplicon data using various primers, and may need adjustment for your specific experimental design.
These scripts, and those in BBTools, should work fine in Bash on Linux, MacOS, or Windows 10.  Other shells generally work fine.

Usage:
These scripts are just guidelines showing how I am processing Covid data.
If you want to use them without modification, follow these steps:

1) Install BBTools (you need the latest version - 38.85+ - due to some new features I added for Covid!).
   a) If you don't have Java, install Java.
      i) If you run java --version on the command line and it reports version 8+, you're fine.
      ii) Otherwise, you can get it from https://openjdk.java.net/install/index.html
   b) Download BBTools 38.85 or higher from https://sourceforge.net/projects/bbmap/
   c) Unzip the archive: tar -xzf BBMap_38.85.tar.gz
   d) Add it to the path: export PATH=/path/to/bbmap/:$PATH
      i) Now you can run the shell scripts from wherever.
   e) Samtools is not necessary, but recommended if you want to make the bam files for visualization.
2) Rename and interleave the files if they are paired, so libraries are in the format "prefix.fq.gz" with one file per library.
3) Copy the Covid reference to the directory.
4) Modify the template processCoronaWrapper.sh:
   a) Pick a library to use for quality-score calibration if desired (line 16).
   b) Add lines, or make a loop, so that processCorona.sh is called on all of the libraries (lines 20-21).
   c) Delete lines 8-9 to let the script run.
5) Run processCoronaWrapper.sh.

Note:
If you need to perform dehosting (removal of human reads), I would suggest doing that first prior to any other data manipulations such as trimming or filtering to allow the resultant dehosted data to be as close to the raw data as possible.
You can remove human reads while maintaining fastq format and keeping pairing and read headers intact with a command like this (where reads.fq can be single-ended or paired/interleaved):
   bbmap.sh ref=hg19.fa in=reads.fq outm=human.fq outu=nonhuman.fq bloom
...where hg19.fa would be a human reference.  This will remove both reads in a pair if either of them map to human.


