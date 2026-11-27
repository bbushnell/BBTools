#!/usr/bin/env bash
# Learn signed CM state residual bounds from completed paths, never realign.
set -eo pipefail
if [[ $# == 0 || ${1:-} == --help ]];then
 printf '%s\n' 'Usage: dev/buildncrnacmcorridors.sh in=parts.tsv calls=training.calls.tsv.gz audit=resolved_assembly_audit.tsv models=training.models.tsv cms=first.cm,second.cm expectedcalls=N expectedsearches=N expectedsteps=N expectedmodels=N expectedassemblies=N out=fresh_directory mincalls=100 mingroups=10 folds=5 uniform=32 trimpercent=1 t=1 -Xmx4g -eoom'
 printf '%s\n' 'parts.tsv has calls and paths columns naming each complete call-summary and aggregate state-path file, including the retained pilot. trimpercent=0..99 removes at most that percentage of whole paths, rounding down and stopping before either support minimum would break; omitted means1percent. Reports source-group cross-validation and visited-state interval widths; these are not caller accuracy, all DP-cell counts or measured speed.'
 printf '%s\n' 'Optional context resources: refs=original_consensus_fastas anchorpadding=original_teacher_padding. FASTAs must cover every census model by exact name, including unrefined members; member exports follow their FASTA order. This writes version2 resources binding actual consensus bases and the anchor context. Both options are required together; omitted preserves version1 output. Saved window geometry is checked against the declared context.'
 printf '%s\n' 'margin=0 (default) adds a selected nonnegative number of bases once to each learned state interval after whole-path selection. Unseen states retain uniform= fallback. For the explicit wide-context experiment use anchorpadding=100 margin=2 with the original FASTAs; the context and bounds must also be honored by the serving caller. Visited-width reports still describe the saved teacher windows, not a different serving input or measured speed. No alignment, genome search or caller is run.'
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=4g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.BuildNcrnaCmCorridors "$@"
