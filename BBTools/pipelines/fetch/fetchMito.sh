#Fetches mitochondrial genomes and annotations from NCBI.

wget -q -O - ftp://ftp.ncbi.nih.gov/genomes/refseq/mitochondrion/*genomic.gbff.gz > mito.genomic.gbff.gz
wget -q -O - ftp://ftp.ncbi.nih.gov/genomes/refseq/mitochondrion/*genomic.fna.gz > mito.genomic.fna.gz
gbff2gff.sh mito.genomic.gbff.gz mito.genomic.gff.gz

