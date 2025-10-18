#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 14, 2025

Description:  Sends taxonomic queries to a remote QuickClade server for
classification.  This client-server architecture allows users to classify
sequences without loading the reference database locally, dramatically reducing
memory requirements and improving performance for multiple queries.  The client
sends sequence data to a remote server running CladeServer with a preloaded
reference database, receives taxonomic classifications, and displays the results.

The client-server design offers several key advantages:
- No local database loading (saves gigabytes of memory)
- Faster startup time (no database initialization)
- Consistent results across multiple users
- Centralized database maintenance and updates
- Ideal for batch processing of many samples

SendClade mirrors the SendSketch architecture and provides the same taxonomic
classification capabilities as QuickClade but with reduced local resource
requirements.  It is particularly useful in compute environments where memory
is limited or when processing many samples sequentially.

Usage Examples:
sendclade.sh in=sequences.fasta
sendclade.sh in=sequences.fasta address=http://myserver.com:3069
sendclade.sh in=sequences.fasta hits=10 oneline out=results.tsv
sendclade.sh in=sequences.fasta local=t mode=perseq minlen=1000
sendclade.sh in=bin1.fa,bin2.fa,bin3.fa hits=5 heap=10

File Parameters:
in=<file,file>  Query files or directories. Input can be fasta, fastq, .clade,
                or .spectra format. Pre-computed .clade/.spectra files are
                sent directly without sequence processing. Multiple files can be
                specified comma-separated, or loose file names are permitted as
                additional arguments.
out=stdout      Output file for results.  If not specified, results are written
                to standard output.  Progress messages always go to stderr.
local=f         Use local server at localhost:5002 instead of the default remote
                server.  Useful for testing or when running your own CladeServer.
address=<url>   Specify custom server address.  Should include full URL with
                protocol and port, e.g., http://myserver.com:3069/clade.
                If protocol is omitted, http:// is assumed.

Basic Parameters:
hits=1          Number of top taxonomic hits to return per query.  More hits
                provide alternative classifications but increase output size.
oneline=f       Print results in tab-delimited format with one line per query.
                Default format is human-readable with detailed information.
                Oneline format includes: QueryName, Q_GC, Q_Bases, Q_Contigs,
                RefName, R_TaxID, R_GC, R_Bases, R_Contigs, R_Level, GCdif,
                STRdif, k3dif, k4dif, k5dif, lineage.
percontig=f     Process each contig/sequence separately instead of combining
                all sequences from each file into a single query.  When true,
                each contig gets its own taxonomic classification.  When false,
                all sequences in a file are combined for classification.
minlen=0        Minimum contig length in percontig mode.  Contigs shorter than
                this threshold are ignored.  Only applies when percontig=true.

Advanced Parameters:
heap=1          Number of intermediate comparison results to store during
                processing.  Higher values may improve accuracy for complex
                queries at the cost of increased memory usage on the server.
printqtid=f     Print query TaxID if present in sequence headers.  Useful for
                benchmarking when query sequences have known taxonomic labels
                in the format 'tid_1234' or similar.
banself=f       Ban self-matches by ignoring records with the same TaxID as
                the query.  Makes the program behave as if that organism is
                not in the reference database.  Useful for testing accuracy.
verbose=f       Enable detailed progress reporting and timing information.
                Shows batch processing, server communication details, and
                performance metrics.

Standard BBTools Parameters:
overwrite=f     Allow overwriting of existing output files.
append=f        Append to existing output files instead of overwriting.

Server Communication:
The default server is: https://bbmapservers.jgi.doe.gov/quickclade
Sequences are sent in batches of up to 100 clades for efficient processing.
The server responds with taxonomic classifications in either human-readable
or tab-delimited format depending on the oneline parameter.

Performance Notes:
SendClade is designed for high-throughput processing.  It batches sequences
efficiently and provides detailed timing information when verbose=true.
Memory usage on the client is minimal as no reference database is loaded.
Server-side processing benefits from preloaded databases and optimized
comparison algorithms.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
"
}

#This block allows symlinked shellscripts to correctly set classpath.
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
CP="$DIR""current/"

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
    # Source the new scripts
    source "$DIR""/memdetect.sh"
    source "$DIR""/javasetup.sh"

    # Parse Java arguments with sendclade-specific defaults
    # For sendclade: fixed 2g memory for client operations
    parseJavaArgs "--mem=2g" "--mode=fixed" "$@"

    # Set environment paths
    setEnvironment
}
calcXmx "$@"

sendclade() {
	local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP clade.SendClade $@"
	echo $CMD >&2
	eval $CMD
}

sendclade "$@"