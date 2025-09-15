#!/bin/bash
#Run this on jgi-web-4

LOG=ntlogVM_32.txt
PASS=xxxxx
DOMAIN=https://nt-sketch.jgi.doe.gov
KILL=https://nt-sketch.jgi.doe.gov/kill/
PORT=3071
DB=nt

nohup /global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -da -Xmx9g port=$PORT verbose tree=auto sketchonly index domain=$DOMAIN killcode=$PASS oldcode=$PASS oldaddress=$KILL $DB k=32,24 1>>$LOG 2>&1 &

#simple mode, for testing:
#/global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -ea -Xmx9g port=3071 verbose tree=auto sketchonly nt k=32,24 index=f
