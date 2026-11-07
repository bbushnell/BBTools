#!/bin/bash

usage(){
echo "
StegoFast: embed arbitrary file bytes or UTF-8 text in synthetic DNA reads.

Usage:
stegofast.sh in=metadata.tar,publication.pdf out=encoded.fq.gz len=150
stegofast.sh literal='This is some metadata' appendto=reads.fq.gz out=combined.fq.gz
stegofast.sh in=mystery.fq.gz decode out=/path/to/extracted/

in=             Files to embed, comma-separated; in= may be repeated.
literal=        UTF-8 text to embed, optionally alongside files.
name=metadata.txt  Filename for literal text.
out=            FASTQ, FASTA, SAM or BAM output; directory in decode mode.
len=150         Bases per metadata read, including envelope (128..1000000).
paired=auto     Emit paired metadata; auto follows append input or out2.
out2=           Optional separate mate output.
appendto=       Existing sequence file to copy before encoded metadata.
appendto2=      Optional second existing mate file.
decode=f        Extract metadata using sequences alone. Ignores ordinary reads.
in2=            Optional second sequence file to scan when decoding.
maxfiles=128    Maximum simultaneously assembled files during decoding.
ow=f            Allow replacing sequence output. Extraction never overwrites.
tin=1 tout=1    Reader and writer threads. t= caps shared thread settings.

Read headers and quality scores are not needed for recovery. Reordering,
reverse complements and identical duplicate parts are supported. All bases
must be preserved. Missing/conflicting parts or failed checksums cause errors.
Arbitrary input files are embedded verbatim, including compression bytes.
See docs/guides/StegoFastGuide.txt for the version-1 sequence format.
"
}
if [[ $# == 0 || "$1" == -h || "$1" == --help ]]; then
	usage
	exit 0
fi

SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
while [[ -h "$SCRIPT" ]]; do
	DIR="$(dirname "$SCRIPT")"
	SCRIPT="$(readlink "$SCRIPT")"
	[[ "$SCRIPT" != /* ]] && SCRIPT="$DIR/$SCRIPT"
done
DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
# New development classes are compiled in current; release jars must include them.
if [[ -f "$DIR/bbtools.jar" ]]; then CP="$DIR/bbtools.jar"; else CP="$DIR/current/"; fi
. "$DIR/javasetup.sh"
. "$DIR/memdetect.sh"
parseJavaArgs "--xmx=1g" "--xms=128m" "--mode=fixed" "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" stego.StegoFast "$@"
