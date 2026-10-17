#!/usr/bin/env bash
set -euo pipefail

out=""
force=f

for arg in "$@"; do
	case "$arg" in
		out=*) out="${arg#out=}" ;;
		force=*) force="${arg#force=}" ;;
		*)
			echo "Unknown argument: $arg" >&2
			echo "Usage: $0 out=/path/to/public_sources [force=t]" >&2
			exit 1
			;;
	esac
done

if [[ -z "$out" ]]; then
	echo "Missing required out=/path/to/public_sources" >&2
	exit 1
fi

case "$force" in
	t|true|1) force=t ;;
	f|false|0) force=f ;;
	*)
		echo "force must be t or f; found $force" >&2
		exit 1
		;;
esac

if ! command -v curl >/dev/null 2>&1; then
	echo "curl not found in PATH" >&2
	exit 1
fi
if ! command -v sha80 >/dev/null 2>&1; then
	echo "sha80 not found in PATH" >&2
	exit 1
fi

mkdir -p "$out"

declare -a rows=(
	"ncbi_pgap_hmm	hmm_PGAP.HMM.tgz	https://ftp.ncbi.nih.gov/hmm/current/hmm_PGAP.HMM.tgz"
	"ncbi_pgap_hmm	hmm_PGAP.tsv	https://ftp.ncbi.nih.gov/hmm/current/hmm_PGAP.tsv"
	"ncbi_pgap_hmm	RELEASE_NOTES.txt	https://ftp.ncbi.nih.gov/hmm/current/RELEASE_NOTES.txt"
	"uniprot_swissprot	uniprot_sprot.fasta.gz	https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
	"uniprot_swissprot	reldate.txt	https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt"
	"uniprot_swissprot	README	https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/README"
	"uniprot_swissprot	LICENSE	https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/LICENSE"
)

for row in "${rows[@]}"; do
	IFS=$'\t' read -r source file url <<< "$row"
	path="$out/$file"
	if [[ -e "$path" && "$force" != "t" ]]; then
		echo "Refusing to overwrite existing file: $path" >&2
		exit 1
	fi
done

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

printf "source\tfile\tbytes\tsha80\turl\n" > "$tmp/manifest.tsv"
for row in "${rows[@]}"; do
	IFS=$'\t' read -r source file url <<< "$row"
	curl -fL --retry 3 --retry-delay 10 -o "$tmp/$file" "$url"
	bytes=$(wc -c < "$tmp/$file" | tr -d ' ')
	hash=$(sha80 "$tmp/$file")
	printf "%s\t%s\t%s\t%s\t%s\n" "$source" "$file" "$bytes" "$hash" "$url" >> "$tmp/manifest.tsv"
done

for row in "${rows[@]}"; do
	IFS=$'\t' read -r source file url <<< "$row"
	mv -f "$tmp/$file" "$out/$file"
done
mv -f "$tmp/manifest.tsv" "$out/sources_manifest.tsv"
