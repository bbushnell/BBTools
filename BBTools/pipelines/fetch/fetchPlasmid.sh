#Fetches plasmid genomes and annotations from NCBI.

wget -q -O - ftp://ftp.ncbi.nih.gov/genomes/refseq/plasmid/*genomic.gbff.gz > plasmid.genomic.gbff.gz
wget -q -O - ftp://ftp.ncbi.nih.gov/genomes/refseq/plasmid/*genomic.fna.gz > plasmid.genomic.fna.gz
gbff2gff.sh plasmid.genomic.gbff.gz plasmid.genomic.gff.gz

