#!/bin/bash
#Run this on jgi-web-1

LOG=taxlogVM_55.txt
PASS=xxxxx
DOMAIN=https://taxonomy.jgi.doe.gov
KILL=https://taxonomy.jgi.doe.gov/kill/
PORT=3068

nohup /global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -da -Xmx31g port=$PORT verbose accession=auto tree=auto table=auto size=auto img=auto pattern=auto prealloc domain=$DOMAIN killcode=$PASS oldcode=$PASS oldaddress=$KILL html 1>>$LOG 2>&1 &

#simple mode, for testing:
#/global/projectb/sandbox/gaag/bbtools/jgi-bbtools/taxserver.sh -ea -Xmx8g port=$PORT verbose accession=null tree=auto table=null
