time fetchproks.sh ftp://ftp.ncbi.nih.gov:21/genomes/refseq/archaea archaea.sh 1>fetchA.o 2>&1
time fetchproks.sh ftp://ftp.ncbi.nih.gov:21/genomes/refseq/bacteria bacteria.sh 1>fetchB.o 2>&1

mkdir archea
cp archaea.sh archaea
cd archaea
sh archaea.sh
cd ..
mkdir bacteria
cp bacteria.sh bacteria
cd bacteria
sh bacteria.sh
cd ..

time nice analyzegenes.sh archaea/*.fna.gz out=archaea.pgm -Xmx1g
time nice analyzegenes.sh bacteria/*.fna.gz out=bacteria.pgm -Xmx1g
time nice analyzegenes.sh */*.fna.gz out=model.pgm -Xmx1g

cutRna.sh
