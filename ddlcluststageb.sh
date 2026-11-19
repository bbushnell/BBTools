#!/usr/bin/env bash
# Collects shard consensuses and rebuilds merged seeds from original members.
usage(){
cat <<'HELP'
Written by Yelan and Keqing.

Collects consensuses from independently clustered shards, or rebuilds merged
consensuses from their original input members. Loads the complete input in memory.

Usage: ddlcluststageb.sh mode=collect codecs=part0.pddl.gz,part1.pddl.gz
                       out=consensuses.faa roster=roster.tsv ow=f
       ddlcluststageb.sh mode=union in=originals.faa partitionmap=map.tsv.gz
                       stageacodecs=part0.pddl.gz,part1.pddl.gz
                       stageamembers=part0.members.tsv.gz,part1.members.tsv.gz
                       roster=roster.tsv stagebcodec=merged.pddl.gz
                       stagebmembers=merged.members.tsv.gz
                       out=seeds.pddl.gz report=union.tsv t=8 ow=f

mode=collect    Emits each shard consensus and its stable shard/family identity.
codecs=         Ordered comma-separated list of shard cluster-sketch files.
roster=         Collection output; union input. Preserve this exact mapping.
mode=union      Attempts one consensus rebuild from each group's original members.
memberranking=t  Default for mutual, layercake, probe, localprobe, localplan and
                liveprobe. Ranks admitted candidates using retained-member length
                and KD with the existing cosine; memberranking=f restores legacy
                ordering. Local selected-baseline modes require serialized member
                descriptors, including last-nonempty values for empty templates.
in=             Complete original protein FASTA, in the partitioner's input order.
partitionmap=   Mapping emitted by ddlclustpartition.sh.
stageacodecs=   Same ordered shard codec list supplied during collection.
stageamembers=  Matching ordered shard membership files, including unassigned rows.
stagebcodec=    Codec from clustering the collected consensuses.
stagebmembers=  Membership file assigning every collected consensus exactly once.
out=            Collected FASTA, or union seed codec, according to mode.
report=         Required for union mode; records rebuild status and source counts.
t=              Workers for independent union rebuilds.
minspan=70      Minimum aligned span when validating a rebuilt model.
minid=0.35      Minimum aligned-query identity fraction.
minblosumratio=0.14  Minimum alignment score ratio under the selected matrix.
matrix=POSSUM55  Loads resources/bbgap/POSSUM55.tsv unless matrixfile is supplied.
                blosum62 uses built-in scores.
matrixfile=     Optional file for POSSUM55; required without an installed default.
                The file's matrix name must match the selector. Do not supply a
                matrixfile with blosum62.
denominator=aligned  Identity denominator: aligned or full query.
maxalignlen=    Optional alignment length limit.
pinright=t     Required stop-anchored alignment.
trimleft=t     Allow left overhangs during alignment.
consensustrimleft=t  Trim unsupported left consensus columns.
trimdepth=0    Consensus depth-trimming fraction.
blacklist=     Required when the codecs declare a blacklist; bytes must match.
blacklistsha80= Optional additional assertion of the codec's blacklist identity.
ow=f           Use explicitly to protect existing outputs.
config=file    Standard BBTools configuration file.
-Xmx           Override Java heap allocation.

Cluster the collected consensuses with passes=0 mincluster=1 topn=0 minlen=70
minspan=70 to retain every representative. Rejected union proposals retain the
merged seed and ALL source members; the report identifies this provisional state.
Use ddlclust.sh seeds=seeds.pddl.gz with passes>=1 to reassign the complete original
input. Serialized source counts are provenance, not new input assignments.
All output directories must exist. Input/output aliases, including hard links,
are rejected. Paired inputs, reads limits and append output are unsupported.
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
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -ea -cp "$CP" prot.cluster.DdlClustStageB "$@"
