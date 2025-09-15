#!/bin/bash
set -e

#Written by Brian Bushnell
#Last updated August 7, 2019

#Fetches and renames RefSeq.
#Be sure the taxonomy server is updated first, or run with local taxonomy data!
#To use this script outside of NERSC, modify $TAXPATH to point to your directory with the BBTools taxonomy data,
#e.g. TAXPATH="/path/to/taxonomy_directory/"

TAXPATH="auto"

#Ensure necessary executables are in your path
#module load pigz

#Fetch RefSeq
#time wget -nv ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/*genomic.fna.gz
#The line below requires pigz!
wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/*genomic.fna.gz | gi2taxid.sh -Xmx1g in=stdin.fa.gz out=renamed.fa.gz pigz=16 unpigz zl=9 server ow maxbadheaders=5000 badheaders=badHeaders.txt bgzip

