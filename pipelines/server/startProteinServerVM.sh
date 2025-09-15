#!/bin/bash
#Run this on jgi-web-5

LOG=proteinlogVM_32.txt
PASS=xxxxx
DOMAIN=https://protein-sketch.jgi.doe.gov
KILL=https://protein-sketch.jgi.doe.gov/kill/
PORT=3074
DB=ProkProt

nohup /global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -da -Xmx16g prealloc=0.9 port=$PORT verbose tree=auto sizemult=2 sketchonly index amino domain=$DOMAIN killcode=$PASS oldcode=$PASS oldaddress=$KILL $DB k=12,9 1>>$LOG 2>&1 &

#simple mode, for testing:
#nohup /global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -ea -Xmx28g port=$PORT verbose tree=auto sizemult=2 sketchonly $DB k=12,9 index=f
