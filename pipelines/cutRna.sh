#!/bin/bash
set -e

#Written by Brian Bushnell
#Last updated August 8, 2019

#This script is designed to cut out the ribosomal and tRNAs from annotated genomes,
#which have paired fna.gz and gff.gz files.  The output is then used to create
#sets of kmers such that every sequence will contain at least one of the kmers in the set.
#In other words, a sequence not sharing a kmer with the set is probably not
#a sequence of that type.  These are used by CallGenes, but Silva is used for the
#ssu and lsu (16S and 23S).


cutgff.sh fastawrap=10k types=rRNA out=23S.fa attributes=23S */*.fna.gz
cutgff.sh fastawrap=10k types=rRNA out=16S.fa attributes=16S */*.fna.gz
cutgff.sh fastawrap=10k types=rRNA out=5S.fa attributes=5S */*.fna.gz
cutgff.sh fastawrap=10k types=tRNA out=tRNA.fa */*.fna.gz

kmerfilterset.sh in=23S.fa k=15 out=23S_15mers.fa ow minkpp=1 maxkpp=1 rcomp=f
kmerfilterset.sh in=16S.fa k=15 out=16S_15mers.fa ow minkpp=1 maxkpp=1 rcomp=f
kmerfilterset.sh in=5S.fa k=9 out=5S_9mers.fa ow minkpp=1 maxkpp=1 rcomp=f
kmerfilterset.sh in=tRNA.fa k=9 out=tRNA_9mers.fa ow minkpp=1 maxkpp=1 rcomp=f
