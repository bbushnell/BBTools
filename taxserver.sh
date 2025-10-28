#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Shijie Yao
Last modified Jan 7, 2020

Description:   Starts a server that translates NCBI taxonomy.

Usage:  taxserver.sh tree=<taxtree file> table=<gitable file> port=<number>

Usage examples:
taxserver.sh tree=tree.taxtree.gz table=gitable.int1d.gz port=1234

At LBL:
taxserver.sh tree=auto table=auto port=1234

For accession number support, add accession=<file,file>  E.g.:

External:
taxserver.sh -Xmx45g tree=tree.taxtree.gz table=gitable.int1d.gz accession=prot.accession2taxid.gz,nucl_wgs.accession2taxid.gz port=1234

At LBL:
taxserver.sh tree=auto table=auto accession=auto port=1234

If all expected files are in some specific location, you can also do this:
taxserver.sh -Xmx45g tree=auto table=auto accession=auto port=1234 taxpath=/path/to/files

To kill remotely, launch with the flag kill=password, then access /kill/password

Parameters:

tree=auto           taxtree path.  Always necessary.
table=auto          gitable path.  Necessary for gi number support.
accession=null      Comma-delimited paths of accession files.
                    Necessary for accession support.
img=null            IMG dump file.
pattern=null        Pattern file, for storing accessions more efficiently.
port=3068           Port number.
domain=             Domain to be displayed in the help message.
                    Default is taxonomy.jgi-psf.org.
dbname=             Set the name of the database in the help message.
sketchcomparethreads=16    Limit compare threads per connection.
sketchloadthreads=4 Limit load threads (for local queries of fastq).
sketchonly=f        Don't hash taxa names.
k=31                Kmer length, 1-32.  To maximize sensitivity and 
                    specificity, dual kmer lengths may be used:  k=31,24
prealloc=f          Preallocate some data structures for faster loading.

Security parameters:

killcode=           Set a password to allow remote killing.
oldcode=            Set the password of a prior instance.
oldaddress=         Attempt to kill a prior instance after initialization,
                    by sending the old code to this address.  For example,
                    taxonomy.jgi-psf.org/kill/
allowremotefileaccess=f   Allow non-internal queries to use internal files
                    for sketching in local mode.
allowlocalhost=f    Consider a query internal if it originates from localhost
                    without being proxied.
addressprefix=128.  Queries originating from this IP address prefix will be
                    considered internal.


Unrecognized parameters with no = symbol will be treated as sketch files.
Other sketch parameters such as index and k are also allowed.
Please consult bbmap/docs/guides/TaxonomyGuide.txt and BBSketchGuide.txt for more information.

Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.

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

z="-Xmx45g"
z2="-Xms45g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
}
calcXmx "$@"

taxserver() {
	local CMD="java $EA $SIMD $EOOM $z $z2 -cp $CP tax.TaxServer $@"
	echo $CMD >&2
	eval $CMD
}

taxserver "$@"
