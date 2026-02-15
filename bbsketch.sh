#!/bin/bash

#For more information, please see sketch.sh
#This exists for people who type bbsketch.sh instead of sketch.sh
#I haven't decided which one will be the canonical version.

resolveSymlinks(){
	SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
}

resolveSymlinks
"$DIR/sketch.sh" "$@"
