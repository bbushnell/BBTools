#!/bin/bash

#Title:       NovaSeq full-lane analysis script.
#Description: Gathers full-lane data for later per-library processing; does not modify any data.
#             Should be run on every raw Illumina lane.
#             The output files PHIX, TILEDUMP, QHIST, and COUNTS should be added to Jamo, linked to the lane.
#             For easiest downstream processing, expected.txt (list of barcodes in the lane)
#             should also be present in Jamo for the lane.
#             Retaining qhist.txt would also be helpful.
#             This script requires PhiX to be spiked in; 0.1% is sufficient but 1%+ is better.
#Author:      Brian Bushnell
#Version:     0.07
#Date:        Aug 17, 2024


#link dev version, or use a container with at least v39.09
#Perlmutter:
##export PATH=$PATH:/global/cfs/cdirs/bbtools/jgi-bbtools/
#Dori:
##export PATH=$PATH:/clusterfs/jgi/groups/gentech/homes/bbushnell/BBTools/


###Set environment variables


#Lane ID (any lane-specific identifier)
LANEID=ABXYZ
#Lane fastq input name
RAW=ABXYZ.1.fq.gz
#Lane input path
RAWPATH=/foo/bar/"$RAW"
#Output path for large fastqs (temporary)
OUT="$PSCRATCH"/"$LANEID"
#Path for recalibrated reads file (temporary)
RECAL="$OUT"/"$LANEID"_recal.fq.gz
#Path for mapped phix (retained)
PHIX=phix.sam.gz
#Path for tile quality dump (retained)
TILEDUMP=tiledump.txt.gz
#Path for barcode counts (retained)
COUNTS=barcodecounts.txt.gz
#Quality histogram (retained)
QHIST=qhist.txt


### BBTools variables

#Physical CPUs
CORES=64
#Compression level.  If you don't have a working version of bgzip in the path set this to 4.
ZL=9
#Shared arguments for all steps
ARGS="t=$CORES zl=$ZL ow"

#Set this to 85% of physical RAM; 48 is for Perlmutter login nodes
MAXRAM=48g
#Keep this at 31
HIGHRAM=31g
#Keep this at 4
LOWRAM=4g

#Max memory flag
MAX="-Xmx$MAXRAM"
#High memory flag
HIGH="-Xmx$HIGHRAM"
#Low memory flag
LOW="-Xmx$LOWRAM"


###Processing starts here

rm finished

#Filter phix reads
bbduk.sh "$LOW" "$ARGS" ref=phix k=25 hdist=2 in="$RAWPATH" outm=phix.fq.gz

#Clumpify reduces phix file size but can make the pipeline less stable when run incorrectly.
#clumpify.sh "$MAX" "$ARGS" in=phix.fq.gz out=clumped.fq.gz

#Adapter-trim phix reads
bbduk.sh "$LOW" "$ARGS" in=phix.fq.gz out=phix_trimmed.fq.gz ref=adapters k=23 mink=11 hdist=2 hdist2=0 tbo tpe ktrim=r minlen=100 ordered

#Align phix reads
#This sam file should be retained for future recalibration if desired
bbmap.sh "$HIGH" "$ARGS" ref=phix nodisk vslow maxindel=100 in=phix_trimmed.fq.gz outm=phix.sam.gz qhist="$QHIST" qahist=qahist.txt mhist=mhist.txt bhist=bhist.txt ordered

#Make recalibration matrices
calctruequality.sh "$HIGH" "$ARGS" in="$PHIX" usetiles callvars ref=phix

#Recalibrate quality
bbduk.sh "$LOW" "$ARGS" in="$RAWPATH" out="$RECAL" recalibrate usetiles

#Filter by tile - this step is best done on the whole lane, after quality recalibration
filterbytile.sh "$MAX" "$ARGS" in="$RECAL" dump="$TILEDUMP"

#Count barcodes
countbarcodes2.sh "$HIGH" "$ARGS" in="$RAWPATH" counts="$COUNTS"

touch finished

#At this point PHIX, DUMP, QHIST, and COUNTS files need to be put in Jamo, associated with the lane.
#Additionally, expected.txt (the file with a 1-per-line list of barcodes in the lane) should also be in Jamo associated with the lane.  Or samplemap.txt (barcode tab sampleid per line).
