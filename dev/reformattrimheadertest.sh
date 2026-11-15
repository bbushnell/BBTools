#!/usr/bin/env bash
# ReformatPacBio outer-pipeline trim/header oracles, with exact untrimmed controls.
# Baseline captures the repeated coordinate update; fixed mode requires the correct trimmed interval.
set -eo pipefail
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
bbtools=$(cd "$(dirname "$SCRIPT")/.." && pwd)
out=; phase=baseline
for arg in "$@"; do
	case "$arg" in
	bbtools=*) bbtools=${arg#bbtools=};;
	outdir=*) out=${arg#outdir=};;
	phase=baseline|phase=fixed) phase=${arg#phase=};;
	*) echo 'Use outdir=fresh-absolute-directory [phase=baseline|fixed] [bbtools=absolute-root]' >&2; exit 1;;
	esac
done
[[ $bbtools = /* && -f $bbtools/javasetup.sh && -d $bbtools/current ]]
[[ -n $out && $out = /* && ! -e $out && -d ${out%/*} ]]
mkdir -- "$out"
out=$(cd "$out" && pwd)
DIR=$bbtools
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed -ea -eoom simd=f
setEnvironment

printf -v retained '%80s' ''
retained=${retained// /C}
printf -v removed '%120s' ''
removed=${removed// /N}
[[ ${#retained} = 80 && ${#removed} = 120 ]]
printf '>movie/11/100_300\n%s%s\n' "$retained" "$removed" > "$out/input-right.fa"
printf '>movie/11/100_300\n%s%s\n' "$removed" "$retained" > "$out/input-left.fa"
cp "$out/input-right.fa" "$out/untrimmed-right.expected.fa"
cp "$out/input-left.fa" "$out/untrimmed-left.expected.fa"
right_header=movie/11/100_180; left_header=movie/11/220_300
if [[ $phase = baseline ]]; then
	right_header=movie/11/60_140; left_header=movie/11/180_260
fi
printf '>%s\n%s\n' "$right_header" "$retained" > "$out/trimmed-right.expected.fa"
printf '>%s\n%s\n' "$left_header" "$retained" > "$out/trimmed-left.expected.fa"

failures=0
printf 'case\texpected_header\tobserved_header\texpected_bases\tobserved_bases\tchild_exit\tresult\n'
for name in untrimmed-right trimmed-right untrimmed-left trimmed-left; do
	layout=${name##*-}
	trim=f; expected_header=movie/11/100_300; expected_bases=200
	if [[ $name = trimmed-* ]]; then
		trim=t; expected_header=$right_header; expected_bases=80
		[[ $layout != left ]] || expected_header=$left_header
	fi
	caseout=$out/$name
	mkdir -- "$caseout"
	command=(java $EA $EOOM $SIMD $XMX $XMS -XX:ActiveProcessorCount=1 -cp "$bbtools/current"
		icecream.ReformatPacBio "in=$out/input-$layout.fa" "out=$caseout/output.fa"
		"trim=$trim" trimpolya=f minlen=0 flaglongreads=f entropy=f
		threads=1 jni=f simd=f pigz=f unpigz=f fastawrap=1000 ow=f)
	printf '%q ' "${command[@]}" > "$caseout/command.txt"
	printf '\n' >> "$caseout/command.txt"
	code=0
	"${command[@]}" > "$caseout/stdout.txt" 2> "$caseout/stderr.txt" || code=$?
	printf '%s\n' "$code" > "$caseout/exit.txt"
	observed_header=missing; observed_bases=missing
	if [[ -f $caseout/output.fa ]]; then
		IFS= read -r observed_header < "$caseout/output.fa" || true
		observed_header=${observed_header#>}
		observed_bases=$(awk '!/^>/{bases+=length($0)} END{print bases+0}' "$caseout/output.fa")
	fi
	result=FAIL
	if [[ $code = 0 ]] && cmp -s "$out/$name.expected.fa" "$caseout/output.fa" &&
			! grep -Eq 'Exception|Error' "$caseout/stderr.txt"; then
		result=PASS
		if [[ $phase = baseline && $trim = t ]]; then result=CAPTURED_WRONG_COORDINATES; fi
	fi
	printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$expected_header" "$observed_header" \
		"$expected_bases" "$observed_bases" "$code" "$result"
	[[ $result != FAIL* ]] || failures=$((failures+1))
done
[[ $failures = 0 ]]
