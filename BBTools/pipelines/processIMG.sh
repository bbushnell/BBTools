#!/bin/bash
set -e

#Written by Brian Bushnell
#Last updated February 21, 2018

#Combines IMG (a JGI genome set) into a single flat file, renamed by taxonomy, and sketches it.
#This script only works on Genepool or other systems connected to NERSC.
#"time" before each command is optional.


#Rename all contigs by prefixing them with the IMG ID, and put them in a single file.
#"in=auto" reads the IMG ID, taxID, and file location from a file at /global/projectb/sandbox/gaag/bbtools/tax/img2/IMG_taxonID_ncbiID_fna.txt
#This is a 3-column tsv file.
#The "imghq" flag uses a different file, "IMG_taxonID_ncbiID_fna_HQ.txt", which contains only high-quality genomes.
time renameimg.sh in=auto imghq out=renamed.fa.gz fastawrap=255 zl=6

#Make the IMG blacklist of uninformative kmers occuring in over 300 different species.
#This is optional, but tends to increase query speed and reduce false positives.
time sketchblacklist.sh -Xmx31g in=renamed.fa.gz prepasses=1 tree=auto taxa taxlevel=species ow out=blacklist_img_species_300.sketch mincount=300 k=31,24 imghq

#Sketch the reference genomes, creating one sketch per IMG ID.
#They are written to 31 files, img0.sketch through img30.sketch.
#The only reason for this is to allow parallel loading by CompareSketch.
time sketch.sh -Xmx31g in=renamed.fa.gz out=img#.sketch files=31 mode=img tree=auto img=auto gi=null ow blacklist=blacklist_img_species_300.sketch k=31,24 imghq

#A query such as contigs.fa can now be compared to the new reference sketches like this:
#comparesketch.sh in=contigs.fa k=31,24 tree=auto img*.sketch blacklist=blacklist_img_species_300.sketch printimg

#On NERSC systems, you can then set the default path to img by pointing /global/projectb/sandbox/gaag/bbtools/img/current at the path to the new sketches.
#Then you can use the default set of img sketches like this:
#comparesketch.sh in=contigs.fa img tree=auto printimg
#That command automatically adds the default path to the sketches, the blacklist, and the correct values for K.
