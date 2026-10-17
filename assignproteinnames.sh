#!/usr/bin/env bash
# Assigns evidence-supported names to protein families from public search results.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
	cat <<'HELP'
Usage: assignproteinnames.sh families=family_identity.tsv ncbi=ncbi.tsv \
    swiss=swissprot.tsv out=family_labels.tsv evidence=decisions.tsv \
    unresolved=unresolved.tsv summary=summary.tsv labelversion=RELEASE \
    identitypin=SHA80 consensuspin=SHA80

Assigns protein-family names from existing public-database search evidence.
No database download or sequence search is performed by this tool.

families=FILE     Authoritative family identity TSV, including family_id,
                  active_index and rep_id. Compressed input is supported.
ncbi=FILE         NCBI model-candidate TSV from BbgapFullEvidenceTriage.
swiss=FILE        Swiss-Prot compact evidence TSV from BbgapLabelEvidence.
reviewed=FILE     Optional reviewed decisions, including explicit abstentions.
reviewedevidence=FILE  Optional reviewed-evidence sidecar for reviewed= rows.
reviewedevidenceout=FILE  Validated copy of the reviewed-evidence sidecar that
                  was projected into schema-2 reviewed labels.
                  eggNOG rows require selected-seed geometry; Swiss-Prot rows
                  require explicit dots in eggNOG-only payload fields.
eggnoghomology=FILE  Optional opt-in eggNOG same-seed Pfam/OG homology leads
                  for otherwise unnamed families; requires provisional=t and
                  eggnoghomologyout=. Rows are weak provisional domain-containing
                  protein hints only, not whole-protein or function assertions.
eggnoghomologyout=FILE  Validated copy of the eggNOG homology candidates with
                  selected=N/Y and ineligibility reasons. This is intentionally
                  separate from the mixed NCBIfam/Swiss-Prot decision evidence.
eggnogdomainrelaxed=f  Extend eggnoghomology= to two reviewed relaxed domain
                  cases: matching single-root rows with an unresolved symbol
                  and no KO/EC, and matching multi-root architecture rows with
                  no Preferred_name/KO/EC functional lead.
synonyms=FILE     Optional explicit product-name equivalences.
mode=exact        Require agreement between the eligible NCBI and Swiss-Prot names.
                  mode=hierarchy uses eligible NCBI model precedence; Swiss-Prot
                  alternatives are retained for review without a mandatory vote.
swissonly=f       Experimental Swiss-Prot-only fallback; requires mode=hierarchy.
                  For otherwise unnamed families with no eligible NCBI naming
                  product, require exactly one coverage-qualified Swiss product
                  comparison key and reject uninformative placeholder products.
minswissidentity=90  Minimum percent identity of the highest-bit-score hit for
                  that product (0-100). Used only when swissonly=t. Identity and
                  accession break score ties; a weaker hit cannot rescue it.
ncbinearquery=f   Experimental near-query NCBI fallback; requires mode=hierarchy.
                  For otherwise insufficient families, retain naming-eligible
                  equivalog/exception rows with native sequence/domain scores,
                  hmm_cov >= mincov and 0.75 <= consensus_cov < mincov, and
                  require exactly one NCBI key matching one full-coverage Swiss key.
ncbinearquerydomains=FILE  Corrected BbgapCoverageAudit coverage.domains.tsv
                  with all retained domain alternatives for ncbinearquery=t.
provisional=f    Emit schema-2 sidecar labels for otherwise unresolved families.
                  Provisional rows stay visibly weak and carry source, coverage,
                  score, product-text and review tags for GFF provenance.
broadconflicts=FILE  Curated exact conflict-set broad-name rules; requires
                  mode=hierarchy. This opt-in pilot can emit
                  domain_supported_broad_name for NCBI top-tier conflicts or
                  homology_supported_broad_name for Swiss-Prot-only conflicts.
                  It does not change minswissidentity or the specific
                  Swiss-Prot-only fallback.
broadparents=FILE  Curated source-scoped parent broad-name rules; requires
                  mode=hierarchy. Every competing product must be explicitly
                  assigned to the same parent or the family remains unresolved.
                  Exact broadconflicts rules take precedence when both are set.
ignorecase=f      Ignore ASCII capitalization when comparing product names and
                  looking up synonyms; preserve original evidence and display names.
out=FILE          Versioned family labels for annotategenes.sh labels=FILE.
evidence=FILE     Per-family decisions and evidence supporting each assignment.
unresolved=FILE   Families whose evidence does not support an assigned name.
summary=FILE      Assignment and unresolved counts.
audit=FILE        Deterministic sample of decisions for independent inspection.
conflictgroups=FILE  Optional mode=hierarchy grouped unresolved-conflict
                  diagnostic, sorted by largest recurring effective NCBI-top-tier
                  or Swiss keyset.
membercoverage=FILE  Optional reference-package member-coverage diagnostic;
                  groups final decisions by named/unresolved bucket and reports
                  members/raw_residues fractions from family_identity.tsv. The
                  raw-residue fraction is . when total raw residues are zero.
memberunresolved=FILE  Optional reference-package unresolved-family diagnostic;
                  lists unresolved families sorted by family member count.
maxaudit=20       Maximum sampled rows per decision stratum.
seed=1            Seed for deterministic audit sampling.
mincov=0.8        Minimum Swiss-Prot query and subject coverage. NCBI eligibility
                  is supplied by the upstream model-candidate evidence table.
ow=t             Allow existing output files to be overwritten. Inputs and
                  outputs must refer to distinct files.
labelversion=ID   Name-assignment release identifier.
identitypin=HEX   Decoded family-identity SHA-256 suffix (20 hexadecimal digits).
consensuspin=HEX  Decoded consensus SHA-256 suffix (20 hexadecimal digits).

Preserves reviewed exceptions. In hierarchy mode, incompatible products at the
highest NCBI precedence remain unresolved; lower-precedence hits do not veto the
selected tier. Specific and broad names have different decision classes. The
experimental hierarchy does not resolve strong Swiss-Prot contradictions.
Swiss-only fallback preserves reviewed decisions and NCBI selections, applies
mincov to both sequence coverages, and records additions as swissprot_full_length.
Its identity threshold is experimental, not a calibrated guarantee of function.
Near-query NCBI fallback is default-off, preserves reviewed decisions, existing
names and conflicts, and records additions separately from ordinary NCBI calls.
Curated broad conflict rules are default-off, preserve the exact source/keyset
provenance, and are not a generic ligand, substrate or regulon deletion rule.
eggNOG homology leads are default-off, require no existing accepted, reviewed or
provisional label, require one coherent same-seed Pfam/OG root, exclude
functional KO/EC/preferred-name rows, and record all candidates separately.
eggnogdomainrelaxed=t remains default-off and only adds symbol-ambiguity or
multi-root architecture cases that otherwise pass the same numeric gates.
Name equivalences must be explicit; similar wording alone does not establish
the same biochemical function.
HELP
	exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")" && pwd)
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
ARGS=(); JVM_ARGS=()
for arg in "$@"; do
	case ${arg,,} in
		--xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
		-ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
		simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("$arg");;
		*) ARGS+=("$arg");;
	esac
done
parseJavaArgs --xmx=1g --xms=256m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$DIR/current" prot.AssignProteinNames "${ARGS[@]}"
