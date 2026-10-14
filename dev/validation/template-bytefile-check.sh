#!/usr/bin/env bash
# Prepare or run A_SampleByteFile template fixtures from a reusable class root.
set -euo pipefail

phase=${1:-}
case "$phase" in
	prep)
		test "$#" = 3
		repo=$2
		work=$3
		test ! -e "$work"
		mkdir -p "$work/src/template" "$work/classes"
		sed 's/		assert(false) : "TODO";/		\/\/ Template customization hook filled by TemplateByteFileCase./' \
			"$repo/current/template/A_SampleByteFile.java" > "$work/src/template/A_SampleByteFile.java"
		javac -cp "$repo/current" -sourcepath "$work/src" -d "$work/classes" \
			"$work/src/template/A_SampleByteFile.java" "$repo/current/test/template/TemplateByteFileCase.java"
		;;
	run)
		test "$#" = 5
		repo=$2
		work=$3
		mode=$4
		expect=$5
		test -d "$work/classes"
		casework="$work/cases/$mode-$expect"
		test ! -e "$casework"
		exec java -ea --add-modules jdk.incubator.vector -Xmx512m -XX:ActiveProcessorCount=1 \
			-cp "$work/classes:$repo/current" test.template.TemplateByteFileCase \
			"classpath=$work/classes:$repo/current" "work=$casework" "mode=$mode" "expect=$expect"
		;;
	all)
		test "$#" = 5
		repo=$2
		work=$3
		mode=$4
		expect=$5
		"$0" prep "$repo" "$work"
		exec "$0" run "$repo" "$work" "$mode" "$expect"
		;;
	*)
		echo "Usage: $0 prep <repo> <work> | run <repo> <work> <mode> <expect> | all <repo> <work> <mode> <expect>" >&2
		exit 2
		;;
esac
