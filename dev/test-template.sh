#!/usr/bin/env bash
# Build and run the finite template regression fixtures. Requires JDK 17+ with Vector API.
# Usage: bash dev/test-template.sh prepare /absolute/fresh/build
#        bash dev/test-template.sh run /absolute/prepared/build
# JAVA, JAVAC and JAVADOC may select the JDK tools. Preparation compiles locally;
# the complete build directory can be transferred to another host for run mode.
set -euo pipefail
mode=${1:?prepare or run}
work=${2:?absolute build directory}
[[ $work == /* ]]
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
java=${JAVA:-java}
javac=${JAVAC:-javac}
javadoc=${JAVADOC:-javadoc}
launch=(nice -n 10 "$java" -ea --add-modules jdk.incubator.vector -Xmx256m -XX:ActiveProcessorCount=1)
prepare_case(){
    local name=$1 kind=$2 source=$3
    mkdir -p "$work/$name/classes"
    "${launch[@]}" -cp "$work/classes" "test.template.Template${kind}Case" \
        mode=prepare "source=$source" "out=$work/$name/src" "${@:4}" \
        > "$work/$name/prepare.stdout" 2> "$work/$name/prepare.stderr"
    "$javac" --release 8 -cp "$work/classes" -sourcepath "$work/empty" \
        -d "$work/$name/classes" "$work/$name/src/template/"*.java \
        > "$work/$name/compile.stdout" 2> "$work/$name/compile.stderr"
}
run_case(){
    local name=$1 rows=$2
    shift 2
    printf '%s\n' "$name" > "$results/current.txt"
    printf '%q ' "${launch[@]}" -cp "$work/classes" "$@" > "$results/$name.command"
    local status=0
    "${launch[@]}" -cp "$work/classes" "$@" > "$results/$name.stdout" 2> "$results/$name.stderr" &
    local child=$!
    printf '%s\n' "$child" > "$results/current.pid"
    wait "$child" || status=$?
    printf '%s\n' "$status" > "$results/$name.exit"
    test "$status" = 0
    test "$(grep -c '^PASS' "$results/$name.stdout")" = "$rows"
    printf '%s\t%s\t0\n' "$name" "$rows" >> "$results/cases.tsv"
}
if [[ $mode == prepare ]]; then
    test ! -e "$work"
    mkdir -p "$work/classes" "$work/empty" "$work/java8"
    "$java" -version > "$work/jdk.txt" 2>&1
    sources=("$repo/current/template/"*.java "$repo/current/test/template/"*.java)
    test "${#sources[@]}" = 34
    # Empty classpath prevents stale shipped classes from satisfying source dependencies.
    "$javac" -J-Xmx2g --add-modules jdk.incubator.vector -cp "$work/empty" \
        -sourcepath "$repo/current" -d "$work/classes" "${sources[@]}" \
        "$repo/current/fileIO/"*.java "$repo/current/parse/"*.java \
        "$repo/current/stream/"*.java "$repo/current/stream/bam/"*.java \
        > "$work/compile.stdout" 2> "$work/compile.stderr"
    "$javac" --release 8 -cp "$work/classes" -sourcepath "$work/empty" \
        -d "$work/java8" "${sources[@]}" > "$work/java8.stdout" 2> "$work/java8.stderr"
    "$javadoc" -quiet -private -Xdoclint:all,-missing -cp "$work/classes" \
        -sourcepath "$work/empty" -tag 'date:a:Date:' -tag 'contributor:a:Contributor:' \
        -d "$work/javadoc" "${sources[@]}" > "$work/javadoc.stdout" 2> "$work/javadoc.stderr"
    prepare_case args-default ArgSplit "$repo/current/template" group=default
    prepare_case args-remaining ArgSplit "$repo/current/template" group=remaining
    prepare_case generator Generator "$repo/current/template/A_Sample_Generator.java"
    prepare_case filter StreamerFilter "$repo/current/template/A_SampleStreamerMT.java"
    prepare_case completion StreamerCompletion "$repo/current"
    prepare_case multiple-completion MultipleInputCompletion "$repo/current"
    prepare_case sam-completion SamStreamerCompletion "$repo/current"
    prepare_case policy Sample2DCompletionPolicy "$repo/current"
    mkdir -p "$work/bytefile/src/template" "$work/bytefile/classes"
    test "$(grep -c 'assert(false) : "TODO";' "$repo/current/template/A_SampleByteFile.java")" = 1
    sed 's/\t\tassert(false) : "TODO";/\t\t\/\/ Finite test supplies the validation hook./' \
        "$repo/current/template/A_SampleByteFile.java" > "$work/bytefile/src/template/A_SampleByteFile.java"
    "$javac" --release 8 -cp "$work/classes" -sourcepath "$work/empty" \
        -d "$work/bytefile/classes" "$work/bytefile/src/template/A_SampleByteFile.java" \
        > "$work/bytefile/compile.stdout" 2> "$work/bytefile/compile.stderr"
    printf 'PASS: all template and test sources compiled; Java8 API/doclint and prepared hooks compiled.\n' > "$work/build.pass"
    cat "$work/build.pass"
elif [[ $mode == run ]]; then
    test -s "$work/build.pass"
    results=$work/results
    test ! -e "$results"
    mkdir "$results"
    trap 'code=$?; printf "%s\n" "$code" > "$results/exit"' EXIT
    printf '%s\n' "$$" > "$results/controller.pid"
    hostname > "$results/host.txt"
    "$java" -version > "$results/jdk.txt" 2>&1
    printf 'case\tpassed_rows\texit\n' > "$results/cases.tsv"
    for group in default remaining; do
        rows=32; [[ $group != remaining ]] || rows=36
        run_case "args-$group" "$rows" test.template.TemplateArgSplitCase mode=run \
            "classes=$work/args-$group/classes" "work=$results/args-$group" expect=fixed "group=$group"
    done
    run_case threadwaiter 9 test.template.TemplateThreadWaiterCase expect=fixed
    run_case multiple 6 test.template.TemplateMultipleInputCase "work=$results/multiple" expect=fixed
    run_case ordering 8 test.template.TemplateStreamerOrderedCase "work=$results/ordering" expect=fixed
    run_case completion 6 test.template.TemplateStreamerCompletionCase mode=run \
        "classes=$work/completion/classes" "work=$results/completion" expect=fixed
    run_case multiple-completion 3 test.template.TemplateMultipleInputCompletionCase mode=run \
        "classes=$work/multiple-completion/classes" "work=$results/multiple-completion" expect=fixed
    run_case sam-completion 3 test.template.TemplateSamStreamerCompletionCase mode=run \
        "classes=$work/sam-completion/classes" "work=$results/sam-completion" expect=fixed
    run_case policy 10 test.template.TemplateSample2DCompletionPolicyCase mode=run \
        "classes=$work/policy/classes" "work=$results/policy" expect=fixed
    for scenario in noout ordinary preexisting alias; do
        expect=pass; [[ $scenario != preexisting && $scenario != alias ]] || expect=reject
        run_case "bytefile-$scenario" 1 test.template.TemplateByteFileCase \
            "classpath=$work/bytefile/classes:$work/classes" "work=$results/bytefile-$scenario" "mode=$scenario" "expect=$expect"
    done
    for scenario in zero one unlimited empty; do
        run_case "lines-$scenario" 1 test.template.TemplateByteFileMaxLinesCase \
            "classpath=$work/bytefile/classes:$work/classes" "work=$results/lines-$scenario" "mode=$scenario" expect=fixed
    done
    for fragments in 1 2; do
        run_case "generator-$fragments" 1 test.template.TemplateGeneratorCase mode=run \
            "classes=$work/generator/classes" "fragments=$fragments" expect=fixed
    done
    for format in fastq sam; do
        for scenario in allkeep drop; do
            run_case "$format-$scenario" 1 test.template.TemplateStreamerFilterCase mode=run \
                "classes=$work/filter/classes" "work=$results/$format-$scenario" "format=$format" "scenario=$scenario" expect=fixed
        done
    done
    run_case validation 4 template.TemplateStreamerValidationCase expect=fixed
    for empty in f t; do
        name=worker-$empty
        "${launch[@]}" -cp "$work/classes" template.TemplateStreamerWorkerCase "empty=$empty" \
            > "$results/$name.stdout" 2> "$results/$name.stderr"
        if [[ $empty == f ]]; then printf 'false\t1\t4\t[0]\t2\ttrue\n';
        else printf 'true\t1\t4\t[0, 1]\t3\ttrue\n'; fi | cmp - "$results/$name.stdout"
        printf '%s\t1\t0\n' "$name" >> "$results/cases.tsv"
    done
    awk -F '\t' 'NR>1{n+=$2} END{exit !(n==133)}' "$results/cases.tsv"
    printf 'PASS:133 finite checks across all14 fixture classes.\n' > "$results/pass.txt"
    cat "$results/pass.txt"
else
    echo 'Expected prepare or run' >&2
    exit 2
fi
