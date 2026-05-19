#!/bin/bash
set -e

# Assembly Preprocessing Pipeline for Bacterial Isolates
# Brian Bushnell, DOE Joint Genome Institute
# Presented at SFAF 2026, Santa Fe, NM
#
# This pipeline preprocesses 2x151bp Illumina paired-end reads for
# bacterial isolate assembly. It performs quality score recalibration,
# filtering, error correction, and read merging, then assembles with
# either Tadpole (fast, included) or SPAdes (external).
#
# Usage: assemblyPipeline.sh <SRA_accession>
#   e.g.: assemblyPipeline.sh SRR36826327
#
# Requirements:
#   - BBTools 39.85+  (https://github.com/bbushnell/BBTools)
#   - SRA Toolkit     (https://github.com/ncbi/sra-tools)
#   - SPAdes 4.x      (optional, for SPAdes assembly)
#
# Reference files needed for contamination filtering:
#   - Illumina.artifacts.fa.gz    (Illumina spike-ins and artifacts)
#   - nextera_LMP_linker.fa.gz    (Nextera linker sequences)
#   - lambda.fa.gz                (Lambda phage)
#   - pJET1.2.fa.gz               (pJET cloning vector)
#   - short.fa.gz                 (Short artifact sequences)
#   - commonMicrobes fused ref    (Common lab contaminant microbes)
#   - mousecatdoghuman ref        (Mammalian host references)
# These are available at: https://portal.nersc.gov/dna/RQC/RQCFilterData/
#
# The pipeline is designed for isolate genomes. For metagenomes,
# some steps (especially the taxID-based microbe filtering) would
# need to be modified or removed.

ACC=$1
if [ -z "$ACC" ]; then
  echo "Usage: assemblyPipeline.sh <SRA_accession>"
  echo "  e.g.: assemblyPipeline.sh SRR36826327"
  exit 1
fi

# Adjust these paths for your system
BB=$(dirname $0)           # Assumes script is in BBTools directory
RQCDATA=/path/to/RQCFilterData
THREADS=$(nproc)
ARGS="t=$THREADS zl=4 ow"
XMX="-Xmx450g"            # Adjust for your system's RAM

DIR=./$ACC
mkdir -p $DIR/ref
cd $DIR

echo "========================================"
echo "Assembly Preprocessing Pipeline"
echo "Accession: $ACC"
echo "Started: $(date)"
echo "Threads: $THREADS"
echo "========================================"


###############################################
# STEP 0: Download from SRA
###############################################
# Restore original Illumina headers with --defline-seq.
# Interleave with reformat.sh (ain=add interleave info, vpair=verify pairing).

echo "=== Step 0: Download from SRA ==="

fastq-dump $ACC --split-files \
  --defline-seq '@$sn $ri:N:0:$sg' --defline-qual '+' \
  --outdir .

$BB/reformat.sh in1=${ACC}_1.fastq in2=${ACC}_2.fastq \
  out=00_raw.fq.gz ain vpair zl=4
rm -f ${ACC}_1.fastq ${ACC}_2.fastq


###############################################
# PHASE I: Build quality-score recalibration
#          matrices from a quick self-assembly
###############################################
# We subsample, trim adapters, assemble with SPAdes, map back,
# and call variants to identify true polymorphisms. Then
# CalcTrueQuality builds per-tile, per-cycle, per-base correction
# matrices, excluding variant sites.
# This allows PhiX-level quality accuracy from any library,
# even when no spike-in reference is available.

echo "=== Phase I: Build recalibration matrices ==="

echo "--- I.1: Subsample 1M reads ---"
$BB/reformat.sh $XMX in=00_raw.fq.gz out=sub.fq.gz reads=1m

echo "--- I.2: Adapter-trim subsample ---"
$BB/bbduk.sh $XMX in=sub.fq.gz out=trimmed.fq.gz \
  ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo ftm=5 ordered

echo "--- I.3: Quick assembly (SPAdes) ---"
spades.py --isolate -k 25,55,95,127 --phred-offset 33 \
  --12 trimmed.fq.gz -o spades_out -t $THREADS -m 450 \
  1>spades_phase1.log 2>&1
cp spades_out/scaffolds.fasta ref/genome.fa

echo "--- I.4: Map reads to assembly ---"
$BB/bbmap.sh $XMX in=trimmed.fq.gz out=mapped.sam.gz \
  ref=ref/genome.fa nodisk ordered t=$THREADS

echo "--- I.5: Call variants ---"
$BB/callvariants.sh $XMX in=mapped.sam.gz ref=ref/genome.fa \
  out=vars.vcf.gz ploidy=1

echo "--- I.6: Build recalibration matrices ---"
$BB/calctruequality.sh $XMX in=mapped.sam.gz ref=ref/genome.fa \
  vcf=vars.vcf.gz usetiles


###############################################
# Detect taxID for microbe exclusion
###############################################
# SendSketch identifies the organism so we can exclude it from
# the microbe contamination reference (otherwise the organism's
# own reads would be filtered out).

echo "=== Detect organism taxID ==="

TAXID=$($BB/sendsketch.sh $XMX in=00_raw.fq.gz reads=400k \
  records=1 json 2>/dev/null | grep TaxID | head -1 | grep -o "[0-9]*")
echo "Detected taxID: $TAXID"

echo "--- Build filtered microbe reference ---"
$BB/filterbytaxa.sh $XMX names=$TAXID tree=$RQCDATA/tree.taxtree.gz \
  level=order \
  in=$RQCDATA/commonMicrobes/fusedERPBBmasked2.fa.gz \
  out=ref/filteredMicrobes.fa.gz ow=t include=f besteffort=f \
  results=microbesUsed.txt


###############################################
# PHASE II: Read preprocessing
###############################################

echo "=== Phase II: Read preprocessing ==="

echo "--- II.1: Clumpify (deduplicate optical/PCR duplicates) ---"
$BB/clumpify.sh $XMX $ARGS in=00_raw.fq.gz out=01_clumped.fq.gz \
  dedupe passes=2 reorder entryfilter groups=1

echo "--- II.2: Quality-score recalibration ---"
$BB/bbduk.sh $XMX $ARGS in=01_clumped.fq.gz out=02_recal.fq.gz \
  recalibrate ordered

echo "--- II.3: FilterByTile (remove low-quality tile regions) ---"
$BB/filterbytile.sh $XMX $ARGS in=02_recal.fq.gz out=03_fbt.fq.gz

echo "--- II.4: Adapter trimming ---"
$BB/bbduk.sh $XMX $ARGS in=03_fbt.fq.gz out=04_trimmed.fq.gz \
  ktrim=r k=23 mink=11 hdist=1 tpe tbo rcomp=f ftm=5 ordered \
  minlen=25 minlenfraction=0.333 ref=adapters

echo "--- II.5: Homopolymer filter ---"
$BB/polyfilter.sh $XMX $ARGS in=04_trimmed.fq.gz out=05_polyfiltered.fq.gz \
  polymers=GC ordered

echo "--- II.6: Artifact filter 1 (spike-ins, vectors, PhiX) ---"
$BB/bbduk.sh $XMX $ARGS in=05_polyfiltered.fq.gz out=06_filtered1.fq.gz \
  k=31 hdist=1 \
  ref=$RQCDATA/Illumina.artifacts.fa.gz,$RQCDATA/nextera_LMP_linker.fa.gz,$BB/resources/phix2.fa.gz,$RQCDATA/lambda.fa.gz,$RQCDATA/pJET1.2.fa.gz \
  ordered maq=5 maxns=0 minlen=25 minlenfraction=0.333 \
  trimpolygleft=6 trimpolygright=6 maxnonpoly=2

echo "--- II.7: Artifact filter 2 (short artifacts) ---"
$BB/bbduk.sh $XMX $ARGS in=06_filtered1.fq.gz out=07_filtered2.fq.gz \
  k=20 hdist=1 ref=$RQCDATA/short.fa.gz ordered

echo "--- II.8: Microbe filter (target organism excluded) ---"
$BB/bbmap.sh $XMX $ARGS in=07_filtered2.fq.gz \
  outu=08_nomicrobe.fq.gz outm=microbes.fq.gz \
  ref=ref/filteredMicrobes.fa.gz nodisk \
  deterministic ordered quickmatch k=13 idtag=t printunmappedcount \
  qtrim=rl trimq=10 untrim ef=0.001 minid=.95 idfilter=.95 \
  maxindel=3 minhits=2 bw=12 bwr=0.16 fast=true \
  maxsites2=10 tipsearch=0 scafstats=commonMicrobes.txt

echo "--- II.9: Mammal filter (cat, dog, human, mouse) ---"
$BB/bbsplit.sh $XMX $ARGS in=08_nomicrobe.fq.gz \
  outu=09_clean.fq.gz outm=human.fq.gz \
  path=$RQCDATA/mousecatdoghuman/ \
  deterministic ordered=false k=14 idtag=t usemodulo \
  printunmappedcount qtrim=rl trimq=10 untrim kfilter=25 \
  maxsites=1 tipsearch=0 minratio=.9 maxindel=3 minhits=2 \
  bw=12 bwr=0.16 fast=true maxsites2=10 build=1 ef=0.03 \
  forcereadonly refstats=refStats.txt


###############################################
# PHASE III: Error correction, merging, assembly
###############################################

echo "=== Phase III: Error correction + merging + assembly ==="

echo "--- III.1: Overlap-based error correction (BBMerge ecco) ---"
$BB/bbmerge.sh $XMX $ARGS in=09_clean.fq.gz out=10_ecco.fq.gz \
  ecco mix adapters=default kfilter=1 k=31 prefilter=1 quantize=3

echo "--- III.2: Clump-based error correction ---"
$BB/clumpify.sh $XMX $ARGS in=10_ecco.fq.gz out=11_eccc.fq.gz \
  ecc passes=4

echo "--- III.3: Kmer-based error correction (Tadpole) ---"
$BB/tadpole.sh $XMX $ARGS in=11_eccc.fq.gz out=12_ecct.fq.gz \
  ecc k=62 wash prefilter=1 tossjunk

echo "--- III.4: Merge overlapping pairs ---"
$BB/bbmerge.sh $XMX $ARGS in=12_ecct.fq.gz \
  outm=13_merged0.fq.gz outu=13_unmerged0.fq.gz \
  kfilter=1 adapters=default prefilter=1 ihist=ihist_merge.txt

echo "--- III.5: Merge with kmer extension (k=93) ---"
$BB/bbmerge.sh $XMX $ARGS in=13_unmerged0.fq.gz extra=13_merged0.fq.gz \
  out=13_merged_ext.fq.gz outu=13_unmerged_ext.fq.gz \
  rem k=93 extend2=100 strict

echo "--- III.6: Combine all merged reads ---"
cat 13_merged0.fq.gz 13_merged_ext.fq.gz > 14_merged.fq.gz

echo "--- III.7: Quality-trim unmerged reads ---"
$BB/bbduk.sh $XMX $ARGS in=13_unmerged_ext.fq.gz out=14_qtrimmed.fq.gz \
  qtrim=r trimq=15 maq=14 minlen=90 ftr=149 maxns=1 ordered

echo "--- III.8a: Assemble with Tadpole (fast, included) ---"
$BB/tadpole.sh $XMX $ARGS \
  in=14_merged.fq.gz,14_qtrimmed.fq.gz \
  out=contigs_tadpole.fa k=145

echo "--- III.8b: Assemble with SPAdes (optional, external) ---"
spades.py --isolate -k 25,55,95,127 --phred-offset 33 --only-assembler \
  --pe-m 1 14_merged.fq.gz --pe-12 1 14_qtrimmed.fq.gz \
  -o spades_final -t $THREADS -m 450 1>spades_final.log 2>&1
cp spades_final/scaffolds.fasta contigs_spades.fa

echo ""
echo "=== Assembly stats ==="
$BB/stats.sh in=contigs_tadpole.fa
echo ""
$BB/stats.sh in=contigs_spades.fa
echo ""
echo "Pipeline complete: $(date)"
