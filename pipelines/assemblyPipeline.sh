#!/bin/bash
set -e

#Written by Brian Bushnell
#Last updated March 4, 2019

#This script is designed to preprocess data for assembly of overlapping 2x150bp reads from Illumina HiSeq 2500.
#Some numbers and steps may need adjustment for different data types or file paths.
#For large genomes, tadpole and bbmerge (during the "Merge" phase) may need the flag "prefilter=2" to avoid running out of memory.
#"prefilter" makes these take twice as long though so don't use it if you have enough memory.
#The "rm temp.fq.gz; ln -s reads.fq.gz temp.fq.gz" is not necessary but added so that any pipeline stage can be easily disabled,
#without affecting the input file name of the next stage.


# --- Setup ---

#Load dependencies.
#These module load commands are for Genepool; getting the correct executables in your path will vary by system.
if [[ $NERSC_HOST == genepool ]]; then
	#module load bbtools
	module load spades/3.9.0
	module load megahit
	module load pigz
	module load quast
elif [[ $NERSC_HOST == denovo ]]; then
	#TODO
elif [[ $NERSC_HOST == cori ]]; then
	#TODO
fi


#Link the interleaved input file as "temp.fq.gz"
rm temp.fq.gz; ln -s reads.fq.gz temp.fq.gz

# --- Preprocessing ---

#Remove optical duplicates
clumpify.sh in=temp.fq.gz out=clumped.fq.gz dedupe optical
rm temp.fq.gz; ln -s clumped.fq.gz temp.fq.gz

#Remove low-quality regions
filterbytile.sh in=temp.fq.gz out=filtered_by_tile.fq.gz
rm temp.fq.gz; ln -s filtered_by_tile.fq.gz temp.fq.gz

#Trim adapters.  Optionally, reads with Ns can be discarded by adding "maxns=0" and reads with really low average quality can be discarded with "maq=8".
bbduk.sh in=temp.fq.gz out=trimmed.fq.gz ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=70 ref=adapters ftm=5 ordered
rm temp.fq.gz; ln -s trimmed.fq.gz temp.fq.gz

#Remove synthetic artifacts and spike-ins by kmer-matching.
bbduk.sh in=temp.fq.gz out=filtered.fq.gz k=31 ref=artifacts,phix ordered cardinality
rm temp.fq.gz; ln -s filtered.fq.gz temp.fq.gz

#Decontamination by mapping can be done here.
#JGI removes these in two phases:
#1) common microbial contaminants (E.coli, Pseudomonas, Delftia, others)
#2) common animal contaminants (Human, cat, dog, mouse)

#Error-correct phase 1
bbmerge.sh in=temp.fq.gz out=ecco.fq.gz ecco mix vstrict ordered ihist=ihist_merge1.txt
rm temp.fq.gz; ln -s ecco.fq.gz temp.fq.gz

#Error-correct phase 2
clumpify.sh in=temp.fq.gz out=eccc.fq.gz ecc passes=4 reorder
rm temp.fq.gz; ln -s eccc.fq.gz temp.fq.gz

#Error-correct phase 3
#Low-depth reads can be discarded here with the "tossjunk", "tossdepth", or "tossuncorrectable" flags.
#For very large datasets, "prefilter=1" or "prefilter=2" can be added to conserve memory.
#Alternatively, bbcms.sh can be used if Tadpole still runs out of memory.
tadpole.sh in=temp.fq.gz out=ecct.fq.gz ecc k=62 ordered
rm temp.fq.gz; ln -s ecct.fq.gz temp.fq.gz

#Normalize
#This phase can be very beneficial for data with uneven coverage like metagenomes, MDA-amplified single cells, and RNA-seq, but is not usually recommended for isolate DNA.
#So normally, this stage should be commented out, as it is here.
#bbnorm.sh in=temp.fq.gz out=normalized.fq.gz target=100 hist=khist.txt peaks=peaks.txt
#rm temp.fq.gz; ln -s normalized.fq.gz temp.fq.gz

#Merge
#This phase handles overlapping reads,
#and also nonoverlapping reads, if there is sufficient coverage and sufficiently short inter-read gaps
#For very large datasets, "prefilter=1" or "prefilter=2" can be added to conserve memory.
bbmerge-auto.sh in=temp.fq.gz out=merged.fq.gz outu=unmerged.fq.gz strict k=93 extend2=80 rem ordered ihist=ihist_merge.txt

#Quality-trim the unmerged reads.
bbduk.sh in=unmerged.fq.gz out=qtrimmed.fq.gz qtrim=r trimq=10 minlen=70 ordered

# --- Assembly ---

#You do not need to assemble with all assemblers, but I have listed the commands for the 3 I use most often

#Assemble with Tadpole
#For very large datasets, "prefilter=1" or "prefilter=2" can be added to conserve memory.
tadpole.sh in=merged.fq.gz,qtrimmed.fq.gz out=tadpole_contigs.fa k=124

#Or assemble with TadWrapper (which automatically finds the best value of K but takes longer)
tadwrapper.sh in=merged.fq.gz,qtrimmed.fq.gz out=tadwrapper_contigs_%.fa outfinal=tadwrapper_contigs k=40,124,217 bisect

#Assemble with Spades
spades.py -k 25,55,95,125 --phred-offset 33 -s merged.fq.gz --12 qtrimmed.fq.gz -o spades_out

#Assemble with Megahit
#Note that the above error correction phases tend to not be beneficial for Megahit
megahit --k-min 45 --k-max 225 --k-step 26 --min-count 2 -r merged.fq.gz --12 qtrimmed.fq.gz -o megahit_out

# --- Evaluation ---

#Evaluate assemblies with AssemblyStats
statswrapper.sh contigs.fa spades_out/scaffolds.fasta megahit_out/contigs.fa format=3 out=

#Evaluate assemblies with Quast (leave out "-R ref.fa if you don't have a reference)
quast.py -f -o quast -R ref.fa tadpole_contigs.fa spades_out/scaffolds.fasta megahit_out/contigs.fa

#Pick which assembly you like best

#Determine the overall taxonomic makeup of the assembly
sendsketch.sh in=tadpole_contigs.fa

#Or, try to determine taxonomy on a per-contig basis.  If this is not sensitive enough, try BLAST instead.
sendsketch.sh in=tadpole_contigs.fa persequence minhits=1 records=4

#Calculate the coverage distribution, and capture reads that did not make it into the assembly
bbmap.sh in=filtered.fq.gz ref=tadpole_contigs.fa nodisk covhist=covhist.txt covstats=covstats.txt outm=assembled.fq.gz outu=unassembled.fq.gz maxindel=200 minid=90 qtrim=10 untrim ambig=all

