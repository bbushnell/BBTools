# Parsing utilities

This package handles command-line settings and reusable delimited-field parsing.
These are separate layers: recognizing an option does not guarantee that a tool
uses the resulting setting.

## Command-line dispatch

`PreParser` handles the early argument pass and configuration expansion.
The program then splits each argument and normally lowercases its name before
calling `Parser.parse(arg, name, value)`. Program-specific handlers may run first;
for example, ReformatStreamer gives its processor a chance before the shared parser.

`Parser.parse` tries static handlers first: JVM options, multi-writer policy,
quality encoding, compression, SAM, FASTA, common global options, histograms and
quality recalibration. It then tries instance handlers for files, common settings,
trimming, interleaving, mapping and cardinality. `parseK` is a separate entry point.
Unrecognized arguments return false so the program can reject or handle them.

Configure before starting worker threads. Both static and instance handlers can
change process-wide settings in other packages; constructing another Parser does
not reset those settings or isolate simultaneous tools. Call `processQuality()`
after argument parsing to apply pending quality offsets.

Configuration files contain one argument per nonempty, non-comment line. They
are expanded at their position in the argument list; later options can override
earlier ones according to the relevant handler. A comma separates config files.
Expansion is one level, not recursive. Equals signs inside configuration filenames
and unparsed argument values survive this expansion. Each tool still owns subsequent
argument splitting; legacy callers such as ReformatStreamer currently split values
again on every equals sign, so this does not repair all input/output filename paths.

## Contracts that are easy to confuse

- `ftr` stores an inclusive rightmost retained coordinate. `ftr2` stores a count
  removed from the right end of each read, including both mates in BBDuk.
- A `trimq` array supplies a sequence of thresholds; BBMerge uses successive
  thresholds for retry attempts, rather than assigning thresholds to base positions.
- `trimclip` concerns alignment-clipped ends. It does not detect adapters.
- `untrim` requests restoration from saved trimming data in supporting callers.
- `testsize` requests throughput reporting based on file sizes; it does not skip work.
- `interleaved=auto` enables detection without marking a forced interleaving choice.
  A caller can still adjust pairing based on separate input/output filenames.
- `Parse.equalsIgnoreCase` with a byte array and a String compares ASCII letters
  without regard to case. Other ASCII characters remain literal; non-ASCII values
  and either null input return false. It does not decode a byte-array encoding.
  The String/String overload retains Java's Unicode case-insensitive comparison.

LineParser variants deliberately retain different scanning and traversal strategies.
Use the individual class contracts for delimiter ordering, empty fields, setters,
reuse and forward-only access. Reuse a parser per worker; do not share its mutable
current-line state between threads. Regression and timing harnesses live in
`test.parse`, outside the production package.
