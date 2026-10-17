#!/usr/bin/env bash
# Prepares downloaded NCBI PGAP HMMs and Swiss-Prot for BBGAP public searches.
set -euo pipefail

sources=
manifest=
out=
hmmer_sif=
diamond_sif=
threads=8

usage(){
	cat >&2 <<'HELP'
Usage: bbgap_public_sources_prepare.sh sources=/path/public_sources out=/fresh/prepared \
  hmmer_sif=/path/hmmer.sif diamond_sif=/path/diamond.sif [threads=8] \
  [manifest=/path/sources_manifest.tsv]

The sources directory is the raw snapshot produced by
dev/bbgap_labels_fetch_sources.sh. This script validates its recorded pins,
concatenates and presses NCBI PGAP HMMs, builds a Swiss-Prot DIAMOND database,
copies source metadata, and writes prepare_manifest.tsv.
HELP
	exit 2
}

for arg in "$@"; do
	case "$arg" in
		sources=*) sources=${arg#*=} ;;
		manifest=*) manifest=${arg#*=} ;;
		out=*) out=${arg#*=} ;;
		hmmer_sif=*) hmmer_sif=${arg#*=} ;;
		diamond_sif=*) diamond_sif=${arg#*=} ;;
		threads=*) threads=${arg#*=} ;;
		-h|--help) usage ;;
		*) printf 'Unknown argument: %s\n' "$arg" >&2; usage ;;
	esac
done

[[ -n $sources && -n $out && -n $hmmer_sif && -n $diamond_sif ]] || usage
[[ -n $manifest ]] || manifest=$sources/sources_manifest.tsv
[[ $threads =~ ^[1-9][0-9]*$ ]] ||
	{ printf 'threads must be a positive integer: %s\n' "$threads" >&2; exit 2; }
[[ ! -e $out ]] || { printf 'Output already exists: %s\n' "$out" >&2; exit 1; }
for file in "$manifest" "$sources/hmm_PGAP.HMM.tgz" "$sources/hmm_PGAP.tsv" \
		"$sources/uniprot_sprot.fasta.gz" "$hmmer_sif" "$diamond_sif"; do
	[[ -s $file ]] || { printf 'Missing required file: %s\n' "$file" >&2; exit 1; }
done
for cmd in apptainer sha80 gzip tar find sort xargs wc stat; do
	command -v "$cmd" >/dev/null ||
		{ printf 'Missing required command: %s\n' "$cmd" >&2; exit 1; }
done

manifest_cell(){
	local file=$1 column=$2
	awk -F '\t' -v file="$file" -v column="$column" '
		NR==1{
			for(i=1; i<=NF; i++){if($i==column){c=i}}
			if($1!="source" || $2!="file" || !c){exit 2}
			next
		}
		$2==file{if(seen++){exit 3}; print $c; found=1}
		END{if(!found){exit 1}}
	' "$manifest"
}

check_source(){
	local file=$1 path=$sources/$1 expected_bytes expected_pin
	expected_bytes=$(manifest_cell "$file" bytes)
	expected_pin=$(manifest_cell "$file" sha80)
	[[ $(stat -c '%s' "$path") == "$expected_bytes" && $(sha80 "$path") == "$expected_pin" ]] ||
		{ printf 'Raw source pin mismatch: %s\n' "$path" >&2; exit 1; }
}

check_source hmm_PGAP.HMM.tgz
check_source hmm_PGAP.tsv
check_source uniprot_sprot.fasta.gz

mkdir -p "$out"
tmp=$(mktemp -d "$out/.tmp.prepare.XXXXXX")
trap 'rm -rf "$tmp"' EXIT

tar -xzf "$sources/hmm_PGAP.HMM.tgz" -C "$tmp"
find "$tmp/hmm_PGAP" -type f -name '*.HMM' | sort > "$tmp/hmm.list"
hmm_count=$(wc -l < "$tmp/hmm.list" | tr -d ' ')
if [[ ! $hmm_count =~ ^[0-9]+$ || $hmm_count -lt 1000 ]]; then
	printf 'Too few HMM files extracted from %s: %s\n' \
		"$sources/hmm_PGAP.HMM.tgz" "$hmm_count" >&2
	exit 1
fi
xargs -a "$tmp/hmm.list" cat > "$tmp/hmm_PGAP.HMM"
apptainer --silent exec "$hmmer_sif" hmmpress "$tmp/hmm_PGAP.HMM" >/dev/null

apptainer --silent exec "$diamond_sif" diamond makedb \
	--in "$sources/uniprot_sprot.fasta.gz" \
	--db "$tmp/swissprot" \
	--threads "$threads" >/dev/null

cp "$sources/hmm_PGAP.tsv" "$tmp/hmm_PGAP.tsv"
cp "$manifest" "$tmp/sources_manifest.tsv"

{
	printf 'file\tbytes\tsha80\n'
	for file in sources_manifest.tsv hmm_PGAP.tsv hmm_PGAP.HMM \
			hmm_PGAP.HMM.h3f hmm_PGAP.HMM.h3i hmm_PGAP.HMM.h3m \
			hmm_PGAP.HMM.h3p swissprot.dmnd; do
		printf '%s\t%s\t%s\n' "$file" \
			"$(stat -c '%s' "$tmp/$file")" \
			"$(sha80 "$tmp/$file")"
	done
} > "$tmp/prepare_manifest.tsv"

for file in sources_manifest.tsv hmm_PGAP.tsv hmm_PGAP.HMM \
		hmm_PGAP.HMM.h3f hmm_PGAP.HMM.h3i hmm_PGAP.HMM.h3m \
		hmm_PGAP.HMM.h3p swissprot.dmnd prepare_manifest.tsv; do
	mv "$tmp/$file" "$out/$file"
done
