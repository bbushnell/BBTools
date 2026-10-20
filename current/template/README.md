# BBTools templates

These files are starting points for new tools. Most are copied and customized;
`BBTool_ST` instead supplies an abstract base for subclasses. An intentional
`assert(false)` or unimplemented hook identifies a decision the new tool's author
must make. Implement that hook with assertions enabled; disabling assertions does
not complete the template.

For new sequence-processing tools, start with `A_SampleStreamerMT` when several
workers process reads, or `A_SampleStreamer` for one processing thread. Both use
the current `Streamer` and `Writer` interfaces. For line-oriented text, use a
`ByteFile` template and a reused `LineParser1` for parsing. The ByteFileMT example
includes a serial first-input pass and a threaded second-input pass; adapt that
shape if the new tool needs only one input.

## Choose a starting point

| Template | Intended use | What to customize or check |
| --- | --- | --- |
| [A_SampleStreamerMT](A_SampleStreamerMT.java) | Parallel sequence processing with shared `Streamer`/`Writer` backends. | Edit the private worker `processReadPair` hook. Keep original batch IDs, including batches emptied by filtering, and do not mutate submitted reads. Workers validate deferred roots and mates before processing; the default hook retains everything. |
| [A_SampleStreamer](A_SampleStreamer.java) | One sequence-processing consumer with configurable I/O backends. | Customize both the Read-pair and SamLine hooks if supporting both representations. SAM/BAM input with SAM/BAM or no output uses SamLine batches. Its Read-mode output-count assertions assume unchanged counts; revisit them for filtering or length changes. |
| [A_SampleByteFile](A_SampleByteFile.java) | One text input with optional valid-line and invalid-line outputs. | Implement validation and the line-processing policy. The example skips empty lines, treats `lines=0` as a zero-line limit and negative `lines` values as unlimited, routes `#` lines to the invalid sink, and emits other lines only up to the first tab. |
| [A_SampleByteFileMT](A_SampleByteFileMT.java) | Serial text prepass followed by threaded processing of a second input. | Implement validation and the worker `processLine` hook. `in1` feeds the prepass and its writers; `in2` feeds shared-ByteFile workers. The private invalid-output path is dormant, not a working CLI option. |
| [A_SampleSamStreamer](A_SampleSamStreamer.java) | Specialized SAM/BAM processing with reference scaffolds and optional worker-local realigners. | Resolve parameter, reference, output-format and processing TODOs. The supplied hook does not invoke the realigner. Output currently uses a generated header rather than forwarding the saved input header; choose the header/reference policy before processing alignments. |
| [A_Sample](A_Sample.java) | Single processing thread using legacy concurrent read input/output streams, including paired reads. | Customize read-pair processing while preserving input/output setup, list return and close-status handling. |
| [A_SampleMT](A_SampleMT.java) | Worker-based paired-read processing using legacy concurrent streams. | Implement the worker's read-pair hook and review validation choices. Keep counters local to workers and merge them after completion. |
| [A_Sample_Unpaired](A_Sample_Unpaired.java) | A single-input, single-output scaffold whose hook and counters operate on read roots. | Implement `process(Read)` and establish an unpaired-input policy. This is not a paired-read transformation example; its counters do not include linked mates. |
| [A_SampleMultipleInput](A_SampleMultipleInput.java) | Several inputs consumed sequentially into one optional output stream. | Customize the per-read loop. Each input owns a separate read stream; one output is shared across inputs. `in=` replaces the input list, while positional inputs append to it. |
| [A_Sample2](A_Sample2.java) | A compact, older single-pass concurrent-stream example. | Customize the read loop and review the compact input/output setup for the new tool. Parsed overwrite/append policy and final stream status are retained. |
| [A_SampleD](A_SampleD.java) | The older distributed/MPI-oriented variant of the compact read example. | Wire MPI initialization/finalization at the explicit main-method reminder before using that entry point. Parsed output policy and final stream status are retained. A non-MPI process test does not validate MPI execution. |
| [A_SampleBasic](A_SampleBasic.java) | Basic byte-line text processing. | Replace the sample's first-tab-prefix copy with the intended text operation and header handling. This is not a FASTQ parser. |
| [A_SampleBasic2](A_SampleBasic2.java) | Output-only text generation with a bounded cycle loop. | Supply generation logic in `doSomething`; the default hook is inert and opens no input stream. |
| [A_SampleSummary](A_SampleSummary.java) | Read scanning followed by an aggregate text report. | Fill the per-read loop and replace `outputResults`; this emits a summary rather than forwarding input reads. |
| [A_Sample_Generator](A_Sample_Generator.java) | Parallel synthetic-read generation, optionally using preloaded template reads. | Implement validation, paired-output choices and `generateRead`. The assertion/TODO stops require a concrete generator's decisions. |
| [BBTool_ST](BBTool_ST.java) | An abstract base for a single-threaded read-processing tool. | Supply defaults, argument parsing, startup/shutdown/statistics hooks, read-pair processing, shared-header policy and a concrete entry point. The private example parser is documentation, not the live parser. |

Legacy concurrent-stream examples remain useful when adapting existing tools.
Their presence does not make their APIs the preferred starting point for every
new tool, and similar examples may preserve different processing shapes.

## Shared helpers

| Helper | Contract |
| --- | --- |
| [Accumulator](Accumulator.java) | Supplies result accumulation, an optional coordination lock and a success report. The interface itself does not synchronize work or discover worker failures. |
| [ThreadWaiter](ThreadWaiter.java) | Starts or joins caller-supplied workers and optionally accumulates results. The accumulator's read lock covers execution/waiting and its write lock covers accumulation; workers must not need that write lock to finish. This is not a general cancellation mechanism. |
| [DoWorker](DoWorker.java) | The single `doWork()` method contract; implementations supply the operation. |
| [ThreadPoolJob](ThreadPoolJob.java) | A subclass supplies `doWork`. The package-local execution path stores the result and returns the input to a bounded queue only after normal completion. Callers handle poison markers before execution; cleanup is not guaranteed after an exception. |

## When making a concrete tool

Keep standard `flag=value` parsing through `PreParser`/`Parser`; local argument
splitting must preserve everything after the first `=`. Choose and document the
units of limits and counters: input fragment limits, individual mates, SAM records
and nonempty text lines are different quantities.

Retain reader/writer ownership and completion checks when changing a processing
hook. Review ordering, paired-output and SAM-header needs together with the chosen
backends. Several examples modify process-wide parsing or I/O settings; copying a
standalone command-line skeleton into a reusable service requires an explicit
policy for those settings and exceptional cleanup.

Validate a concrete hook with small literal inputs, exact expected outputs and
nearby unchanged controls before measuring speed. Empty input, paired counts,
filtering and output-disabled operation deserve checks when supported. Keep
worker state local and batch transfers coarse; use BBTools byte parsers, builders
and primitive collections where their contracts fit. A template is a starting
implementation, not proof that a customized tool is correct or fast.
