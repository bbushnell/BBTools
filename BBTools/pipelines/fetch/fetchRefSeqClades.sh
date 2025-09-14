#!/bin/bash
set -e

#Written by Brian Bushnell
#Last updated August 7, 2019

#Fetches specific RefSeq clades, plus nt and nr.
#Be sure taxonomy is updated first!
#To use this script outside of NERSC, modify $TAXPATH to point to your directory with the BBTools taxonomy data,
#e.g. TAXPATH="/path/to/taxonomy_directory/"

TAXPATH="auto"

#Ensure necessary executables are in your path
#module load pigz

#These commands require pigz!

time wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/invertebrate/*.genomic.fna.gz | gi2taxid.sh -Xmx1g in=stdin.fna.gz out=refseq.invertebrate.fna.gz pigz=6 unpigz=t bgzip=t preferbgzip=t zl=9 server ow maxbadheaders=5000 badheaders=badHeaders_invertebrate.txt taxpath=$TAXPATH
time wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/*.genomic.fna.gz | gi2taxid.sh -Xmx1g in=stdin.fna.gz out=refseq.plant.fna.gz pigz=6 unpigz=t bgzip=t preferbgzip=t zl=9 server ow maxbadheaders=5000 badheaders=badHeaders_plant.txt taxpath=$TAXPATH
time wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/*.genomic.fna.gz | gi2taxid.sh -Xmx1g in=stdin.fna.gz out=refseq.plasmid.fna.gz pigz=6 unpigz=t bgzip=t preferbgzip=t zl=9 server ow maxbadheaders=5000 badheaders=badHeaders_plasmid.txt taxpath=$TAXPATH

time wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/*.genomic.fna.gz | gi2taxid.sh -Xmx1g in=stdin.fna.gz out=refseq.plastid.fna.gz pigz=6 unpigz=t bgzip=t preferbgzip=t zl=9 server ow maxbadheaders=5000 badheaders=badHeaders_plastid.txt taxpath=$TAXPATH
time wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/protozoa/*.genomic.fna.gz | gi2taxid.sh -Xmx1g in=stdin.fna.gz out=refseq.protozoa.fna.gz pigz=6 unpigz=t bgzip=t preferbgzip=t zl=9 server ow maxbadheaders=5000 badheaders=badHeaders_protozoa.txt taxpath=$TAXPATH
time wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/*.genomic.fna.gz | gi2taxid.sh -Xmx1g in=stdin.fna.gz out=refseq.mitochondrion.fna.gz pigz=6 unpigz=t bgzip=t preferbgzip=t zl=9 server ow maxbadheaders=5000 badheaders=badHeaders_mitochondrion.txt taxpath=$TAXPATH

time wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/*.genomic.fna.gz | gi2taxid.sh -Xmx1g in=stdin.fna.gz out=refseq.vertebrate_mammalian.fna.gz pigz=6 unpigz=t bgzip=t preferbgzip=t zl=9 server ow maxbadheaders=5000 badheaders=badHeaders_vertebrate_mammalian.txt taxpath=$TAXPATH
time wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_other/*.genomic.fna.gz | gi2taxid.sh -Xmx1g in=stdin.fna.gz out=refseq.vertebrate_other.fna.gz pigz=6 unpigz=t bgzip=t preferbgzip=t zl=9 server ow maxbadheaders=5000 badheaders=badHeaders_vertebrate_other.txt taxpath=$TAXPATH
time wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/archaea/*.genomic.fna.gz | gi2taxid.sh -Xmx1g in=stdin.fna.gz out=refseq.archaea.fna.gz pigz=6 unpigz=t bgzip=t preferbgzip=t zl=9 server ow maxbadheaders=5000 badheaders=badHeaders_archaea.txt taxpath=$TAXPATH

time wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/*.genomic.fna.gz | gi2taxid.sh -Xmx1g in=stdin.fna.gz out=refseq.bacteria.fna.gz pigz=6 unpigz=t bgzip=t preferbgzip=t zl=9 server ow maxbadheaders=5000 badheaders=badHeaders_bacteria.txt taxpath=$TAXPATH
time wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/*.genomic.fna.gz | gi2taxid.sh -Xmx1g in=stdin.fna.gz out=refseq.fungi.fna.gz pigz=6 unpigz=t bgzip=t preferbgzip=t zl=9 server ow maxbadheaders=5000 badheaders=badHeaders_fungi.txt taxpath=$TAXPATH
time wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/*.genomic.fna.gz | gi2taxid.sh -Xmx1g in=stdin.fna.gz out=refseq.viral.fna.gz pigz=6 unpigz=t bgzip=t preferbgzip=t zl=9 server ow maxbadheaders=5000 badheaders=badHeaders_viral.txt taxpath=$TAXPATH

time wget -q -O - ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz | gi2taxid.sh -Xmx1g in=stdin.fna.gz out=nt.fna.gz pigz=6 unpigz=t bgzip=t preferbgzip=t zl=9 server ow maxbadheaders=5000 badheaders=badHeaders_nt.txt taxpath=$TAXPATH
time wget -q -O - ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz | gi2taxid.sh -Xmx1g in=stdin.faa.gz out=nr.faa.gz pigz=6 unpigz=t bgzip=t preferbgzip=t zl=9 server ow maxbadheaders=-1 badheaders=badHeaders_nr.txt taxpath=$TAXPATH

touch done
