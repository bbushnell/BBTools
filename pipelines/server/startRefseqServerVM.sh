#!/bin/bash
#Run this on jgi-web-4

LOG=refseqlogVM_32.txt
PASS=xxxxx
DOMAIN=https://refseq-sketch.jgi.doe.gov
KILL=https://refseq-sketch.jgi.doe.gov/kill/
PORT=3072
DB=RefSeq

nohup /global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -da -Xmx28g prealloc=0.9 port=$PORT verbose tree=auto sizemult=2 sketchonly index domain=$DOMAIN killcode=$PASS oldcode=$PASS oldaddress=$KILL $DB k=32,24 1>>$LOG 2>&1 &

#simple mode, for testing:
#nohup /global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -ea -Xmx28g port=3072 verbose tree=auto sizemult=2 sketchonly RefSeq k=32,24 index=t
