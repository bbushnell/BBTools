#!/bin/bash
# Plan whole-genome search batches; write reversible input-qualified FASTAs.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden
Last modified November 9, 2026

Usage: dev/preparecmsearchbatches.sh mode=plan in=genomes.tsv out=new_plan tasks=16
       dev/preparecmsearchbatches.sh mode=write in=plan_directory out=new_task task=0
       dev/preparecmsearchbatches.sh mode=schedule in=plan_directory out=new_schedule
           modelcosts=paired.model_costs.tsv samplebases=303039386 workers=8

Plan mode reads the five-column whole-genome census only: input_id, fasta,
file_bytes, nucleotide_bases, fasta_records. It balances bases across tasks and
packs whole genomes into task-local batches. It opens no sequence files.
Write mode reads a complete saved plan, then writes only the selected task.
No searches or biological classifications occur. Output directories must be new.

mode=plan       plan, write or schedule.
in=             Census file (plan), or completed plan directory (write).
out=            Fresh output directory. Overwrite/append are unsupported.
tasks=16        Planning task count; cannot exceed the selected genome count.
maxgenomes=64   Maximum whole genomes per planned batch.
maxbases=400m   Maximum bases per planned batch; an oversized genome fails.
include=        Optional planning subset: first TSV field is an absolute FASTA
                path, no header. Paths must be unique and all match the census.
task=           Required zero-based task ID in write mode. Saved plan limits apply.
modelcosts=     Schedule-only TSV: accession and candidate_CPU_seconds with the
                shared samplebases= below, OR accession, CPU_seconds, sample_bases,
                status with a base count per model. The latter requires status
                MEASURED or MEASURED_REUSED and forbids a global samplebases=.
                Extra columns are ignored. Each baseline must be >=10 CPU seconds.
samplebases=   Required only for legacy candidate_CPU_seconds one-worker timings.
workers=8      Independent one-CPU search slots per task in schedule mode.
taskperbatch=f Schedule each existing batch as its own task when true. Preserves
                batch IDs and saved preparation; does not change the input plan.
concurrenttasks=0  Modeled concurrent task slots; zero permits all tasks at once.
extraseconds=0 Additional unmeasured CPU seconds per search, for sensitivity only.
                Pilot CPU already includes startup; this is not a second measured
                startup charge or an estimate of fixed startup extracted from it.

Plan outputs: members.tsv, batches.tsv, tasks.tsv and final summary.tsv.
Schedule outputs: jobs.tsv (expensive-first queues), tasks.tsv, workers.tsv and
summary.tsv, plus task_starts.tsv with each task's modeled start/end and slot.
Jobs scale pilot CPU linearly by batch bases using each model's measured sample,
then assign each next job to the earliest available modeled worker. Tasks start
in ID order on the earliest free task slot, honoring concurrenttasks=. The actual
scheduler may choose another order. Summary sample_bases=0 means mixed sample sizes.
Reports distinguish CPU, makespan, idle slot-seconds, longest point and tail above
ideal division. These are projections, not measured timings: queue wait, variable
genome difficulty, I/O, contention, preparation and collection are excluded.
No genome files are opened and no searches are run by schedule mode.
Write outputs: batch_N/batch.fa, targets.tsv, members.tsv and final summary.tsv;
plus a final task summary. FASTAs are seekable and shared across model searches.
Each contig stays separate, with unchanged sequence bytes/orientation. Headers
become i<positive_input_id>c<one_based_contig_ordinal>. The literal TSV sidecar maps
that ID to input, ordinal, original target token and length. Percent signs remain
literal. Description text after native ASCII whitespace is not part of the token.
The current native reader requires printable ASCII target tokens; encoding loss
fails. Repeated tokens across inputs are valid; repeated tokens within one input
fail. Source file size and actual contig/base totals must match the census.
The sequence reader/writer uses one consumer; external tasks supply parallelism.
Partial output from a failed write has no completion summary and is not reusable
as a successful batch. The tool never deletes source or prepared FASTAs.
HELP
}
if [ $# -eq 0 ] || [ "$1" = -h ] || [ "$1" = --help ]; then
	usage
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ]; do
	cd "$(dirname "$SCRIPT")"
	RESOLVED_SCRIPT="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$RESOLVED_SCRIPT"
done
DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)/"
popd > /dev/null
if [ -f "${DIR}bbtools.jar" ]; then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=512m --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.PrepareCmSearchBatches "$@"
