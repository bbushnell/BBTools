#!/bin/bash
# Generates the combined ribosomal DDL blacklist for FindSSU.
# Builds separate blacklists for 16S, 18S, and ITS at multiple
# taxonomic levels, then merges them into one file.
#
# Cutoffs tuned May 24, 2026:
#   family/200, order/70, class/24, phylum/8
# These were validated on ITS (35k fungi) and give a good balance
# of genus-level signal retention vs deep-level noise removal.
#
# Usage: bash dev/make_ribo_blacklist.sh
# Requires: resources/all_prok_16S_best_taxsorted.fa.gz
#           resources/all_euk_18S_best_taxsorted.fa.gz
#           resources/all_ITS_best_taxsorted.fa.gz
#
# Author: Brian Bushnell, Noire
# Date: May 24, 2026

set -e
DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$DIR"

TMPDIR="${TMPDIR:-/tmp}/ribo_blacklist_$$"
mkdir -p "$TMPDIR"
echo "Working in $TMPDIR"

LEVELS="family/200 order/70 class/24 phylum/8"

for TYPE in 16S 18S ITS; do
    case $TYPE in
        16S) REF="resources/all_prok_16S_best_taxsorted.fa.gz" ;;
        18S) REF="resources/all_euk_18S_best_taxsorted.fa.gz" ;;
        ITS) REF="resources/all_ITS_best_taxsorted.fa.gz" ;;
    esac

    echo ""
    echo "===== $TYPE: Building 256-bucket sketches with kmers ====="
    bash ddlwriter.sh in="$REF" out="$TMPDIR/${TYPE}_sketches.tsv" \
        k=19 buckets=256 exponent=4 mode=persequence kmers=t

    BLFILES=""
    for LEVELCUT in $LEVELS; do
        LEVEL="${LEVELCUT%%/*}"
        CUT="${LEVELCUT#*/}"
        OUTFA="$TMPDIR/bl_${TYPE}_${LEVEL}${CUT}.fa"
        echo "--- $TYPE: $LEVEL/$CUT ---"
        bash ddlblacklist.sh in="$TMPDIR/${TYPE}_sketches.tsv" \
            out="$OUTFA" k=19 exponent=4 \
            mintaxcount="$CUT" level="$LEVEL"
        BLFILES="$BLFILES,$OUTFA"
    done
    BLFILES="${BLFILES#,}"

    echo "--- $TYPE: merging ---"
    java -ea -cp current/ --add-modules jdk.incubator.vector \
        ddl.MergeDDLBlacklists in="$BLFILES" \
        out="$TMPDIR/bl_${TYPE}_merged.fa" k=19 ow=t
done

echo ""
echo "===== Merging all three types ====="
java -ea -cp current/ --add-modules jdk.incubator.vector \
    ddl.MergeDDLBlacklists \
    in="$TMPDIR/bl_16S_merged.fa,$TMPDIR/bl_18S_merged.fa,$TMPDIR/bl_ITS_merged.fa" \
    out="$TMPDIR/riboDDLBlacklist.fa" k=19 ow=t

echo ""
echo "===== Validating combined blacklist ====="
for TYPE in 16S 18S ITS; do
    case $TYPE in
        16S) REF="resources/all_prok_16S_best_taxsorted.fa.gz" ;;
        18S) REF="resources/all_euk_18S_best_taxsorted.fa.gz" ;;
        ITS) REF="resources/all_ITS_best_taxsorted.fa.gz" ;;
    esac
    echo ""
    echo "--- $TYPE baseline ---"
    bash ddlwriter.sh in="$REF" out="$TMPDIR/${TYPE}_baseline.tsv" \
        k=19 buckets=128 exponent=4 mode=persequence
    bash ddlblacklist.sh in="$TMPDIR/${TYPE}_baseline.tsv" out=/dev/null \
        k=19 exponent=4 validate=t samples=200000

    echo "--- $TYPE with blacklist ---"
    bash ddlwriter.sh in="$REF" out="$TMPDIR/${TYPE}_filtered.tsv" \
        k=19 buckets=128 exponent=4 mode=persequence \
        blacklist="$TMPDIR/riboDDLBlacklist.fa"
    bash ddlblacklist.sh in="$TMPDIR/${TYPE}_filtered.tsv" out=/dev/null \
        k=19 exponent=4 validate=t samples=200000
done

echo ""
echo "===== Done ====="
echo "Combined blacklist: $TMPDIR/riboDDLBlacklist.fa"
echo "Copy to resources with: cp $TMPDIR/riboDDLBlacklist.fa resources/"
