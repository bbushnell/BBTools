#!/usr/bin/env bash
# One-pair no-adapter output oracle for the existing RemoveAdapters2/3 split finding.
# Baseline captures RA2 duplication and RA3's single-end assertion separately.
# Fixed mode requires both programs to reject paired splitting before CLI output construction.
# Forced interleaving uses a configured-main test entry; these legacy CLIs lack an interleaved flag.
set -eo pipefail
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
bbtools=$(cd "$(dirname "$SCRIPT")/.." && pwd)
out=; phase=baseline; program=both
for arg in "$@"; do
	case "$arg" in
	bbtools=*) bbtools=${arg#bbtools=};;
	outdir=*) out=${arg#outdir=};;
	phase=baseline|phase=fixed) phase=${arg#phase=};;
	program=both|program=2|program=3) program=${arg#program=};;
	*) echo 'Use outdir=new-directory [bbtools=absolute-root] [phase=baseline|fixed] [program=both|2|3]' >&2; exit 1;;
	esac
done
[[ $bbtools = /* && -f $bbtools/javasetup.sh && -d $bbtools/current && -n $out && ! -e $out && -d ${out%/*} ]]
mkdir -- "$out"
out=$(cd "$out" && pwd)
# Required when this helper is copied beside results rather than run from dev/.
DIR=$bbtools
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed -ea -eoom simd=f
setEnvironment

left=AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
right=TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
quality=IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
adapter=CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
[[ ${#left} = 64 && ${#right} = 64 && ${#quality} = 64 && ${#adapter} = 64 ]]
printf '@pair/1\n%s\n+\n%s\n' "$left" "$quality" > "$out/read1.fq"
printf '@pair/2\n%s\n+\n%s\n' "$right" "$quality" > "$out/read2.fq"
cat "$out/read1.fq" "$out/read2.fq" > "$out/paired.expected.fq"
cat "$out/read1.fq" "$out/read2.fq" "$out/read2.fq" "$out/read1.fq" > "$out/duplicated.expected.fq"
programs=(2 3)
[[ $program = both ]] || programs=("$program")
failures=0
printf 'program\tlayout\tsplit\texpected_records\tobserved_records\texpected_bases\tobserved_bases\tchild_exit\tresult\n'
for version in "${programs[@]}"; do
	for layout in unpaired paired detected forced-main separate; do
		for split in f t; do
			[[ $layout != separate || $split = f ]] || continue
			caseout=$out/ra$version-$layout-$split
			mkdir "$caseout"
			inputflags=("in=$out/read1.fq")
			entry=("pacbio.RemoveAdapters$version")
			output=$caseout/output.fq
			outputs=("$output")
			expected_files=("$out/read1.fq")
			expected=$out/read1.fq; expected_records=1; expected_bases=64
			if [[ $layout != unpaired ]]; then
				inputflags=("in=$out/read1.fq" "in2=$out/read2.fq")
				expected=$out/paired.expected.fq; expected_records=2; expected_bases=128
				expected_files=("$expected")
			fi
			if [[ $layout = detected || $layout = forced-main ]]; then
				inputflags=("in=$out/paired.expected.fq")
			fi
			if [[ $layout = forced-main ]]; then
				entry=(test.pacbio.RemoveAdaptersForcedMainTest "program=$version")
			elif [[ $layout = separate ]]; then
				output=$caseout/output_#.fq
				outputs=("$caseout/output_1.fq" "$caseout/output_2.fq")
				expected_files=("$out/read1.fq" "$out/read2.fq")
			fi
			if [[ $phase = baseline && $version = 3 && $layout != unpaired && $split = t ]]; then
				if [[ $layout != paired ]]; then
					# The same helper probe cannot establish a second input-layout baseline.
					printf '3\t%s\tt\t2\tnot-run\t128\tnot-run\tnot-run\tSKIPPED_UNSAFE_BASELINE_CLI\n' "$layout"
					continue
				fi
				# Never enter the old worker/writer path that can swallow a paired-input assertion.
				probe=(java $EA $EOOM $SIMD $XMX $XMS -XX:ActiveProcessorCount=1 -cp "$bbtools/current"
					test.pacbio.RemoveAdapters3PairedSplitTest)
				printf '%q ' "${probe[@]}" > "$caseout/command.txt"
				printf '\n' >> "$caseout/command.txt"
				probe_code=0
				"${probe[@]}" > "$caseout/stdout.txt" 2> "$caseout/stderr.txt" || probe_code=$?
				printf '%s\n' "$probe_code" > "$caseout/exit.txt"
				if [[ $probe_code = 0 ]] && grep -Fq $'paired-split-guard\tAssertionError\tAssertionError\tCAPTURED' "$caseout/stdout.txt"; then
					printf '3\tpaired\tt\t2\tnot-written\t128\tnot-written\t0\tCAPTURED_LIBRARY_ASSERTION\n'
				else
					printf '3\tpaired\tt\t2\tnot-written\t128\tnot-written\t%s\tFAIL_LIBRARY_PROBE\n' "$probe_code"
					failures=$((failures+1))
				fi
				continue
			fi
			command=(java $EA $EOOM $SIMD $XMX $XMS -XX:ActiveProcessorCount=1 -cp "$bbtools/current"
				"${entry[@]}" "${inputflags[@]}"
				"out=$output" "adapter=$adapter" "split=$split" minratio=1
				threads=1 pigz=f unpigz=f qin=33 qout=33 ow=f)
			printf '%q ' "${command[@]}" > "$caseout/command.txt"
			printf '\n' >> "$caseout/command.txt"
			code=0
			"${command[@]}" > "$caseout/stdout.txt" 2> "$caseout/stderr.txt" || code=$?
			printf '%s\n' "$code" > "$caseout/exit.txt"
			if [[ $phase = fixed && $layout != unpaired && $split = t ]]; then
				# Match the main exception, not the echoed command's input path or split flag.
				result=FAIL_EARLY_REJECTION
				actual_records=not-written; actual_bases=not-written
				if compgen -G "$caseout/output*.fq" > /dev/null; then
					actual_records=unexpected-output; actual_bases=unexpected-output
				fi
				if [[ $code = 1 ]] &&
						grep -Eq '^Exception in thread "main" java\.lang\.(IllegalArgumentException|RuntimeException): .*([Pp]aired.*split|split.*[Pp]aired)' "$caseout/stderr.txt" &&
						grep -Fq 'split=f' "$caseout/stderr.txt" &&
						[[ $actual_records = not-written ]] &&
						! grep -Eq '^Started (read stream|output threads|[0-9]+ processing thread)' "$caseout/stdout.txt" "$caseout/stderr.txt"; then
					result=PASS_EARLY_REJECTION
				fi
				printf '%s\t%s\tt\t0\t%s\t0\t%s\t%s\t%s\n' "$version" "$layout" "$actual_records" "$actual_bases" "$code" "$result"
				[[ $result != FAIL* ]] || failures=$((failures+1))
				continue
			fi
			actual_records=missing; actual_bases=missing
			bytes_match=t; files_present=t
			for index in "${!outputs[@]}"; do
				cmp -s "${expected_files[$index]}" "${outputs[$index]}" || bytes_match=f
				if [[ ! -f ${outputs[$index]} ]]; then files_present=f; fi
			done
			if [[ $files_present = t ]]; then
				read -r actual_records actual_bases < <(awk 'FNR%4==2{bases+=length($0)} {lines++} END{print lines/4, bases+0}' "${outputs[@]}")
			fi
			result=FAIL
			# Both strands of a C-only query miss the A-only and T-only inputs.
			# Require the real detector's no-hit counters, rather than assuming that path.
			if ! awk '/^Plus adapters:|^Minus adapters:/{seen++; if($3!=0){bad=1}} END{exit(bad || seen!=2)}' "$caseout/stdout.txt"; then
				result=FAIL_ADAPTER_CONTROL
			elif [[ $phase = baseline && $layout != unpaired && $split = t && $version = 2 ]]; then
				if [[ $code = 0 && $actual_records = 4 && $actual_bases = 256 ]] &&
						cmp -s "$out/duplicated.expected.fq" "$output"; then
					result=CAPTURED_DUPLICATES
				fi
			elif [[ $code = 0 && $actual_records = "$expected_records" && $actual_bases = "$expected_bases" ]] &&
					[[ $bytes_match = t ]] &&
					! grep -Eq 'Exception|Error' "$caseout/stderr.txt"; then
				result=PASS
			fi
			# RA2 reports output roots plus reachable mates. Check that report against actual FASTQ records/bases.
			if [[ $version = 2 && $result != FAIL* ]]; then
				if ! awk -v reads="$actual_records" -v bases="$actual_bases" '/^Reads Out:/{gsub(/[()]/, "", $4); ok=($3==reads && $4==bases)} END{exit(!ok)}' "$caseout/stdout.txt"; then
					result=FAIL_OUTPUT_STATISTICS
				fi
			fi
			printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$version" "$layout" "$split" \
				"$expected_records" "$actual_records" "$expected_bases" "$actual_bases" "$code" "$result"
			[[ $result != FAIL* ]] || failures=$((failures+1))
		done
	done
done
[[ $failures = 0 ]]
