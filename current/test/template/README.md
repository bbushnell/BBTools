# Template regression fixtures

These finite checks accompany the template maintenance changes. Twelve fixture
classes live in `test.template`. `template.TemplateStreamerWorkerCase` and
`template.TemplateStreamerValidationCase` stay in `template` because they directly
construct the package-private `A_SampleStreamerMT.ProcessThread` and inspect its
package-private counters. Production visibility is not widened for testing.

From a BBTools source checkout, with a JDK containing the incubating Vector API
(JDK 17 or newer):

```bash
bash dev/test-template.sh prepare /tmp/template-tests
bash dev/test-template.sh run /tmp/template-tests
```

The work directory must be new for preparation. Run mode requires a prepared build
without an existing `results` directory. The complete prepared directory can be
copied to a compute host and run there; no compiler is needed on that host. Select
JDK tools through `JAVA`, `JAVAC`, and `JAVADOC` when they are not on PATH. Prepare
mode compiles the full template and test packages plus the I/O dependencies from
source, checks the changed packages against Java 8 APIs, and checks Javadocs. This
is not a claim that the runner supports a Java 8 runtime.

Preparation creates separate source copies for examples with intentional hooks.
It fills only the literal hooks or redirects named I/O factories/close callbacks;
the original source files and their assertions are preserved. The prepared copies
remain in package `template`; fixture callbacks are explicitly qualified with
`test.template`. Each fixture checks its expected records, counters, exceptions,
or metadata and fails on a mismatch. The runner requires 133 successful checks
across all fourteen classes, with assertions enabled.

Coverage includes local argument decoding, output descriptors, optional text
output, line limits, generated-read accounting, filtering, normal/empty batches,
validation before writer handoff, I/O completion status and accumulator lock
release. Constructor tests do not prove actual file-writing policy; finite fake
I/O checks do not certify arbitrary failures or cancellation. Performance is not
measured by these tiny fixtures.

Outputs and diagnostics are retained under the work directory. `results/cases.tsv`
records successful case counts; `results/exit` records the suite exit status.
The individual fixtures retain historical baseline expectations for reproducing
older defects. The portable runner uses the fixed expectations for current code.
Historical source revisions and their original commands remain in project Git
history; old source paths should be used with those historical revisions.
