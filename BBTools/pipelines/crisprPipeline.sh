#!/bin/bash
set -e

#Written by Brian Bushnell
#Last updated September 20, 2023

#This script is designed to preprocess data for CRISPR detection in Illumina 2x150bp libraries.
#The "rm temp.fq.gz; ln -s reads.fq.gz temp.fq.gz" is not necessary but added so that any pipeline stage can be easily disabled,
#without affecting the input file name of the next stage.


# --- Setup ---

#Link the interleaved input file as "temp.fq.gz"
#rm temp.fq.gz #Only if already present
ln -s reads.fq.gz temp.fq.gz

# --- Preprocessing (for raw data only) ---

#Remove reads from low-quality regions of the flowcell
filterbytile.sh in=temp.fq.gz out=filtered_by_tile.fq.gz
rm temp.fq.gz; ln -s filtered_by_tile.fq.gz temp.fq.gz

#Trim adapters.  Optionally, reads with Ns can be discarded by adding "maxns=0" and reads with really low average quality can be discarded with "maq=8".
bbduk.sh in=temp.fq.gz out=trimmed.fq.gz ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=90 ref=adapters ftm=5
rm temp.fq.gz; ln -s trimmed.fq.gz temp.fq.gz

#Remove synthetic artifacts and spike-ins by kmer-matching.
bbduk.sh in=temp.fq.gz out=filtered.fq.gz k=31 ref=artifacts,phix cardinality
rm temp.fq.gz; ln -s filtered.fq.gz temp.fq.gz

#It would be nice to do quality-score correction here but that requires a reference so we will skip it.

# --- Merging and Error Correction ---

#Merge the reads first since they will be merged at some point anyway. 'mix' puts the merged and unmerged reads in the same file.
bbmerge.sh in=temp.fq.gz out=merged.fq.gz mix strict
rm temp.fq.gz; ln -s merged.fq.gz temp.fq.gz

#Error-correct phase 1.  Note the int=f flag because the file is no longer interleaved.
#This is optional but seems safe and effective.
clumpify.sh in=temp.fq.gz int=f out=eccc.fq.gz ecc conservative passes=9
rm temp.fq.gz; ln -s eccc.fq.gz temp.fq.gz

#Error-correct phase 2
#If this runs out of memory just skip it.
#This is optional but seems safe and effective.
tadpole.sh in=temp.fq.gz int=f out=ecct.fq.gz ecc k=72 conservative
rm temp.fq.gz; ln -s ecct.fq.gz temp.fq.gz

# --- CRISPR-Finding ---

#With a reference of known repeats:
bbcrisprfinder.sh in=temp.fq.gz int=f outc=crisprs.fq outr=repeats.fa outs=spacers.fa chist=chist.txt phist=phist.txt ref=knownRepeats.fa outref=uses.fa

#- or -

#With no reference:
bbcrisprfinder.sh in=temp.fq.gz int=f outc=crisprs.fq outr=repeats.fa outs=spacers.fa chist=chist.txt phist=phist.txt ow


#In either case, you can subsequently run additional passes using output from the first pass to refine the results.
#Without a reference, 3 total passes are advised (including the initial pass); with a reference, 1 or 2 total passes are advised.
#Here, mincount=3 will restrict the tool to use repeats encountered at least 3 times in the first pass as references,
#to prevent polluting the reference pool with spurious repeats.  In practice, lower seems better (even 1 is fine),
#but it can become much slower for pass 2 if mincount is set low.

bbcrisprfinder.sh in=temp.fq.gz int=f outc=crisprs2.fq outr=repeats2.fa outs=spacers2.fa chist=chist2.txt phist=phist2.txt ref=repeats.fa mincount=3
bbcrisprfinder.sh in=temp.fq.gz int=f outc=crisprs3.fq outr=repeats3.fa outs=spacers3.fa chist=chist3.txt phist=phist3.txt ref=repeats2.fa mincount=3
