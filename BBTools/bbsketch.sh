#!/bin/bash

#For more information, please see sketch.sh
#This exists for people who type bbsketch.sh instead of sketch.sh
#I haven't decided which one will be the canonical version.

pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"

"$DIR"sketch.sh $@
