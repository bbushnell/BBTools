#!/bin/bash
# Experimental multi-K: in=reads out=NEW.fa k=96,124,64,32 assemblek=96
# Optional livenet=/livecutoff= and postnet=/postcutoff= pairs. No truth input.
# graphoverlapguard=f: experimental final-walk veto for exact overlaps from
# graphK through assembleK-1. Requires an already requested final graph phase.
# It never repairs a join, adds a count-table load, or changes earlier phases.
# initialmergecoverageratio=0: disabled by default; a finite ratio >=1 applies
# the existing direct-merge coverage veto to initial cleanup only (try 1.75).
# It does not alter contig growth, the postpass NN, later-K policy, or table loads.
# fuseshortevidence=f: experimental veto of selected long-contig fusions when an
# eligible sub-mincontig endpoint has an equal or longer overlap. Short ends are
# evidence only; they cannot create pairs or break cycles to admit new pairs.
# graphbranchguard=f: experimental final-only veto of implicit graph walks when
# either consumed end is an outward kmer branch, including zero/short paths.
# Explicit overlap edges retain their policies. Requires an existing final
# graph phase; no new table load, initial cleanup change, or default change.
# initialbranchguard=f: opt-in use of the same branch veto during initial
# ordinary graph cleanup only, after the postpass. It restores the caller's
# policy before later-K phases and on failure; no extra table is loaded.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=4g --xms=4g --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.MultiKBranchRescue "$@"
