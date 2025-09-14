#Fetches plastid genomes and annotations from NCBI.

wget -q -O - ftp://ftp.ncbi.nih.gov/genomes/refseq/plastid/*genomic.gbff.gz > plastid.genomic.gbff.gz
wget -q -O - ftp://ftp.ncbi.nih.gov/genomes/refseq/plastid/*genomic.fna.gz > plastid.genomic.fna.gz
gbff2gff.sh plastid.genomic.gbff.gz plastid.genomic.gff.gz

