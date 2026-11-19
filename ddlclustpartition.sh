#!/usr/bin/env bash
# Balanced seeded protein partitions for independent clustering jobs.
usage(){
cat <<'HELP'
Written by Keqing and Fischl.

Splits protein FASTA records into balanced files without clustering or filtering.
Loads all inputs into memory. Partition sizes differ by at most one record;
the first partitions receive any remainder, so divisibility is not required.

Usage: ddlclustpartition.sh in=tails.faa out=part#.faa map=partition.tsv
                           ways=4 order=shuffle seed=1

in=             One disk FASTA or an ordered comma-separated list of FASTA files.
                Paired input and reads limits are rejected.
out=            FASTA pattern with exactly one #, replaced by zero-based shard ID.
map=            TSV: parent_input_index, shard, shard_input_index, input_id.
ways=4          Positive number of balanced partitions; some may be empty.
order=shuffle   Seeded shuffle, then balanced contiguous blocks.
                hydro sorts the shuffled records by ascending hydrophobicity,
                preserving seeded ties; undefined hydrophobicity sorts last.
seed=1          Reproducible shuffle seed, also used to order hydrophobicity ties.
tailcrop=0      No cropping by default. A positive value retains that many final
                amino acids before sorting; a terminal * is preserved but not
                counted. Any shorter record causes failure instead of filtering.
ready=          Optional completion-file pattern with exactly one #. Each new
                file signals its closed FASTA shard, not completion of the map
                or whole partition. Existing completion files are always rejected.
ow=f            Existing output files require explicit overwrite permission.

All records are retained, including short proteins and duplicate text IDs.
Headers and retained residue bytes, including case, ambiguity symbols and
terminal stops, are preserved; line wrapping may change. Parent ordinals follow
the declared input-file order and distinguish duplicate text IDs.
All output parent directories must already exist. Input/output aliases,
including symbolic and hard links, are rejected before any output is opened.
No alignment, consensus construction or assignment is performed.
config=file     Standard BBTools configuration file.
-Xmx            Override Java heap allocation.
HELP
}
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then usage; exit 0; fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")" && pwd)
if [[ -f $DIR/bbtools.jar ]]; then CP=$DIR/bbtools.jar; else CP=$DIR/current; fi
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
parseJavaArgs --xmx=4g --xms=32m --mode=auto "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -ea -cp "$CP" prot.cluster.DdlClustPartition "$@"
