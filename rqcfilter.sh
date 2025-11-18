#!/bin/bash

#For more information, please see rqcfilter2.sh
#RQCFilter was deprecated and replaced with RQCFilter2

resolveSymlinks(){
	SCRIPT="$0"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
}

resolveSymlinks
"$DIR/rqcfilter2.sh" "$@"
