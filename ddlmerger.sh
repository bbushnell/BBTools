#!/bin/bash

usage(){
echo "
Written by Ady
Last modified July 26, 2026

Description:  Loads multiple DDL files, optionally condenses to fewer buckets,
merges records sharing a TID (e.g., combining mitochondrial/plastid/plasmid
with host genomes), sorts, renumbers, and writes a single combined DDL file.

Source categories are detected from filenames:
  *mito* = mitochondrial, *plastid* = plastid, *plasmid* = plasmid

Usage:  ddlmerger.sh in=a.ddl.gz,b.ddl.gz out=combined.ddl.gz
    or: ddlmerger.sh in=big.tsv.gz out=small.tsv.gz buckets=2048 merge=f
    or: ddlmerger.sh in=all.ddl.gz out=filtered.ddl.gz tidsout=keep.txt \\
          minsize=refseq_bacteria.fa.gz,100000+refseq_invertebrate.fa.gz,1000000

Parameters:
in=<file>           Input DDL file(s), comma-delimited.
                    Can also be given as bare arguments.
out=<file>          Output combined DDL file.
k=31                K-mer length (for validation).
buckets=-1          Target bucket count.  Input DDLs with more buckets
                    will be condensed by folding (taking max per bucket).
                    -1 disables condensing.
merge=t             Merge all categories with matching TIDs.
mergemito=t         Merge mitochondrial records into host genome.
mergeplastid=t      Merge plastid records into host genome.
mergeplasmid=t      Merge plasmid records into host genome.
exponent=5          Exponent bit width within each 16-bit bucket value; the
                    mantissa gets the remaining 16-exponent bits.  An input
                    file's #exponent header overrides this, and inputs that
                    disagree with each other are a fatal error.
blacklist=<file>    Kmer blacklist; condensing prefers non-blacklisted kmers,
                    baking the blacklist into the output.
overwrite=f         Overwrite existing output.
verbose=f           Print per-file load counts and categories.

Size filtering:
minsize=<list>      Per-source-file minimum genome size, as key,bases pairs
                    joined by '+':
                      minsize=refseq_bacteria.fa.gz,100000+refseq_fungi.fa.gz,1m
                    A record whose source file matches a key and whose genome is
                    smaller than that many bases is dropped at load.  Sizes accept
                    k/m/g suffixes.  The key is matched against each record's
                    origin (its source filename), ignoring directory and the
                    .gz/.tsv/.ddl/.fa/.fna extensions.
                    A source file named by NO rule is kept whole.  That is how
                    organelles, plasmids and viruses stay exempt: do not list
                    them.  No clade logic is involved, so 'plastid' vs
                    'chloroplast' naming does not matter.
                    The run reports, per rule, how many records it matched and
                    dropped; a rule matching NOTHING is called out, since that
                    usually means an input file was renamed.  Origins covered by
                    no rule are listed too.
                    NOTE: the separator is '+', not ';'.  BBTools launchers run
                    'eval', so a ';' inside any argument is read by bash as a
                    command separator and everything after it is silently lost
                    from the java command line.  A '~' is also unsafe (bash
                    expands it).  Repeating the flag is always safe and needs no
                    separator: minsize=a,100000 minsize=b,1000000
tidsout=<file>      Write the taxIDs surviving to the output, sorted and unique,
                    one per line.  Feed this to cladeloader.sh whitelist= to make
                    a spectra database congruent with this sketch database.
"
}

if [ -z "$1" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	usage
	exit
fi

resolveSymlinks(){
	SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=4g" "--xms=400m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ddl.DDLMerger $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
