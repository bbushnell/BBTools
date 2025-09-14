#Isolate assembly script for producing poly-G-artifact-free assemblies.
#For Perlmutter
#Version 1.1
#Brian Bushnell
#October 30, 2024

#link dev version
export PATH=$PATH:/global/cfs/cdirs/bbtools/jgi-bbtools/

#Link mousecatdoghuman index (this is the Perlmutter location)
MCDH=/global/cfs/cdirs/bbtools/RQCFilterData_Local/mousecatdoghuman/

#Link your RAW (not filtered) reads
ln -s path/to/reads raw.fq.gz


#Set variables
CORES=64
ZL=6
HIGHRAM=31g
#MAXRAM here is for login nodes; scheduled jobs should set this to 85% of physical memory requested
MAXRAM=48g

ARGS="t=$CORES zl=$ZL ow"

#Max memory commands
MAX="-Xmx$MAXRAM"

#High memory commands
HIGH="-Xmx$HIGHRAM"

#Low memory commands
LOW="-Xmx4g"


#Find adapter sequence
bbmerge.sh "$LOW" "$ARGS" in=raw.fq.gz outa=adapters.fa

#Trim adapters
bbduk.sh "$LOW" "$ARGS" in=raw.fq.gz out=raw_trimmed.fq.gz tbo tpe hdist=2 k=23 mink=9 hdist2=1 ref=adapters.fa minlen=135 ktrim=r

#Deduplicate
clumpify.sh "$MAX" "$ARGS" in=raw_trimmed.fq.gz out=deduped.fq.gz passes=4 groups=1 dedupe optical dist=50

#Filter artifacts
bbduk.sh "$LOW" "$ARGS" ref=artifacts,phix literal=AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA in=deduped.fq.gz k=31 hdist=1 out=filtered.fq.gz

#Filter human
bbsplit.sh "$MAX" "$ARGS" deterministic ordered=false k=14 usemodulo printunmappedcount kfilter=25 maxsites=1 tipsearch=0 minratio=.9 maxindel=3 minhits=2 bw=12 bwr=0.16 fast maxsites2=10 build=1 ef=0.03 bloomfilter bloomk=29 bloomhashes=1 bloomminhits=6 bloomserial path="$MCDH" refstats=refStats.txt forcereadonly in=filtered.fq.gz out=clean.fq.gz outm=human.fq.gz

#Filter microbes: Disabled due to potential false positives
###bbmap.sh "$HIGH" "$ARGS" deterministic quickmatch k=13 idtag=t printunmappedcount qtrim=rl trimq=10 untrim ef=0.001 path=/global/cfs/cdirs/bbtools/RQCFilterData_Local/commonMicrobes/ pigz=t unpigz=t "$ZL" minid=.95 idfilter=.95 maxindel=3 minhits=2 bw=12 bwr=0.16 fast=true maxsites2=10 build=1 tipsearch=0 in=filtered.fq.gz out=nonmicrobe.fq.gz outm=microbes.fq.gz scafstats=commonMicrobes.txt

#make quick assembly
tadpole.sh "$HIGH" "$ARGS" in=clean.fq.gz out=qecc.fq.gz k=62 merge wash ecc tossjunk tu ldf=0.4 tossdepth=1 aecc
tadpole.sh "$HIGH" "$ARGS" in=qecc.fq.gz out=quick.fa k=124 mcs=5 mce=4 merge

#map trimmed reads to quick assembly
bbmap.sh "$HIGH" "$ARGS" in=clean.fq.gz out=clean.sam.gz vslow maxindel=40 ref=quick.fa

#Make recal matrices
calctruequality.sh "$HIGH" "$ARGS" in=clean.sam.gz usetiles ref=quick.fa callvars

#Recalibrate quality
bbduk.sh "$LOW" "$ARGS" in=clean.fq.gz out=clean_recal_tile.fq.gz recalibrate usetiles

#Filter low-quality flowcell areas
filterbytile.sh "$HIGH" "$ARGS" in=clean_recal_tile.fq.gz out=fbt_recal_tile.fq.gz lowqualityonly=t

#Primary poly-G filter
polyfilter.sh "$HIGH" "$ARGS" in=fbt_recal_tile.fq.gz out=polyfilter_fbt_recal_tile.fq.gz

#Trim residual poly-G reads; no longer needed since polyfilter does trimming and kmer matching now.
#bbduk.sh "$LOW" "$ARGS" in=polyfilter_fbt_recal_tile.fq.gz trimpolyg=6 trimpolyc=6 maxnonpoly=2 minlen=135 literal=GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG k=29 hdist=2 out=hdist2.fq.gz


#Poly-G removal and other generic steps are above this line
#----------------------------------------------------------
#Isolate-assembly-specific preprocessing is below this line


#Error-correct phase 1
bbmerge.sh "$HIGH" "$ARGS" in=hdist2.fq.gz out=ecco.fq.gz ecco mix adapters=adapters.fa kfilter=1 k=31

#Error-correct phase 3; note k has to be less than half of read length here
tadpole.sh "$HIGH" "$ARGS" in=ecco.fq.gz out=ecct.fq.gz ecc k=62 wash tu tossdepth=1 ldf=0.15

#Read merging pass 1
bbmerge.sh "$HIGH" "$ARGS" in=ecct.fq.gz outm=merged0.fq.gz outu=unmerged0.fq.gz kfilter=1 adapters=adapters.fa

#Read merging with kmer extension
bbmerge.sh "$HIGH" "$ARGS" in=unmerged0.fq.gz extra=merged0.fq.gz out=merged_rem.fq.gz outu=unmerged_rem.fq.gz rem k=124 extend2=120
bbmerge.sh "$HIGH" "$ARGS" in=unmerged_rem.fq.gz extra=merged0.fq.gz,merged_rem.fq.gz out=merged_rem2.fq.gz outu=unmerged_rem2.fq.gz rem k=145 extend2=140
bbmerge.sh "$HIGH" "$ARGS" in=unmerged_rem2.fq.gz extra=merged0.fq.gz,merged_rem.fq.gz,merged_rem2.fq.gz out=merged_rem3.fq.gz outu=unmerged_rem3.fq.gz rem k=93 extend2=100 strict

#Combine merged read files
zcat merged0.fq.gz merged_rem.fq.gz merged_rem2.fq.gz merged_rem3.fq.gz | reformat.sh "$LOW" "$ARGS" in=stdin.fq int=f out=merged_both.fq.gz

#Q-trim unmerged reads
bbduk.sh "$LOW" "$ARGS" in=unmerged_rem3.fq.gz out=qtrimmed.fq.gz qtrim=r trimq=15 cardinality cardinalityout maq=14 minlen=90 ftr=149 maxns=1

#Extend reads
tadpole.sh "$HIGH" "$ARGS" in=merged_both.fq.gz extra=qtrimmed.fq.gz out=merged_ext1.fq.gz mode=extend mce=4 er=10 el=10 k=145
tadpole.sh "$HIGH" "$ARGS" in=qtrimmed.fq.gz extra=merged_both.fq.gz out=unmerged_ext1.fq.gz mode=extend mce=4 el=10 k=124

#Assemble with Spades
shifter --image=staphb/spades:4.0.0 spades.py -t "$CORES" -k 25,55,95,127 --phred-offset 33 --only-assembler --isolate --pe-m 1 merged_ext1.fq.gz --pe-12 1 unmerged_ext1.fq.gz -o spades_out 1>spades.o 2>&1

#Test for residual poly-G contam
bbduk.sh "$LOW" "$ARGS" in=spades_out/contigs.fasta literal=GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG hdist=2 k=25

#Test contiguity
stats.sh in=spades_out/contigs.fasta

