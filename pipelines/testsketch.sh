#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 17, 2019

Description:  Runs SendSketch at various taxonomic levels 
for sensitivity analysis.

Usage: sketchtest.sh <file> <taxid> <db> <optional extra flag>
"
}

if [ -z "$1" ] || [ -z "$2" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

QUERY=$1
TID=$2
DB=$3
EXTRA=$4
EXTRA2=$5

echo ""

echo "********** Strain **********"

sendsketch.sh $QUERY $DB -Xmx1g colors=f records=1 minhits=1 silent banunclassified printcommonancestorlevel tid=$TID include=$TID $EXTRA $EXTRA2

echo "********** Species **********"

sendsketch.sh $QUERY $DB -Xmx1g colors=f records=1 minhits=1 silent banunclassified printcommonancestorlevel tid=$TID include=$TID includelevel=species exclude=$TID excludelevel=strain $EXTRA $EXTRA2

echo "********** Genus **********"

sendsketch.sh $QUERY $DB -Xmx1g colors=f records=1 minhits=1 silent banunclassified printcommonancestorlevel tid=$TID include=$TID includelevel=genus exclude=$TID excludelevel=species $EXTRA $EXTRA2

echo "********** Family **********"

sendsketch.sh $QUERY $DB -Xmx1g colors=f records=1 minhits=1 silent banunclassified printcommonancestorlevel tid=$TID include=$TID includelevel=family exclude=$TID excludelevel=genus $EXTRA $EXTRA2

echo "********** Order **********"

sendsketch.sh $QUERY $DB -Xmx1g colors=f records=1 minhits=1 silent banunclassified printcommonancestorlevel tid=$TID include=$TID includelevel=order exclude=$TID excludelevel=family $EXTRA $EXTRA2

echo "********** Class **********"

sendsketch.sh $QUERY $DB -Xmx1g colors=f records=1 minhits=1 silent banunclassified printcommonancestorlevel tid=$TID include=$TID includelevel=class exclude=$TID excludelevel=order $EXTRA $EXTRA2

echo "********** Phylum **********"

sendsketch.sh $QUERY $DB -Xmx1g colors=f records=1 minhits=1 silent banunclassified printcommonancestorlevel tid=$TID include=$TID includelevel=phylum exclude=$TID excludelevel=class $EXTRA $EXTRA2

echo "********** Superkingdom **********"

sendsketch.sh $QUERY $DB -Xmx1g colors=f records=1 minhits=1 silent banunclassified printcommonancestorlevel tid=$TID include=$TID includelevel=superkingdom exclude=$TID excludelevel=phylum $EXTRA $EXTRA2

echo "********** Life **********"

sendsketch.sh $QUERY $DB -Xmx1g colors=f records=1 minhits=1 silent banunclassified printcommonancestorlevel tid=$TID include=$TID includelevel=life exclude=$TID excludelevel=superkingdom $EXTRA $EXTRA2

echo ""
