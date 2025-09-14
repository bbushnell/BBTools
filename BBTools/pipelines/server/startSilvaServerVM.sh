#!/bin/bash
#Run this on jgi-web-4

LOG=ribologVM_32.txt
PASS=xxxxx
DOMAIN=https://ribo-sketch.jgi.doe.gov
KILL=https://ribo-sketch.jgi.doe.gov/kill/
PORT=3073
REF=/global/projectb/sandbox/gaag/bbtools/silva/latest/both_seq#.sketch
DB=silva

nohup /global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -da -Xmx10g port=$PORT verbose tree=auto sketchonly index whitelist domain=$DOMAIN killcode=$PASS oldcode=$PASS oldaddress=$KILL ref=$REF dbname=Silva blacklist=silva k=32,24 1>>$LOG 2>&1 &

#simple mode, for testing:
#/global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -ea -Xmx10g port=3073 verbose tree=auto sketchonly silva k=32,24 index=f domain=https://ribo-sketch.jgi-psf.org
