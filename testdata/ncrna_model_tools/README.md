# Small ncRNA model-tool examples

These files are synthetic teaching data derived from the first consensus in the
public BBTools `resources/callgenes_6S/RF01685/consensus.fa` library. The complete
library, paired models and seeds are already distributed with BBTools; no private
corpus or external download is needed.

`genes.fa` contains twenty differently named copies of that sequence.
`contexts.fa` adds twenty `ACGT` flank bases at each end. The roster/locus files
declare twenty **synthetic** source groups and positive TaxIDs. Those identifiers,
the constant GC metadata and the teaching panel are not biological annotations.
`graph.tsv` is the native ORIGINAL support graph at identity 1 against all nineteen
installed models. `panel.fa` contains nineteen synthetic contigs: one shipped
consensus per contig with 200 independently generated bases on each side, using
Java Random seed 42. `panel.truth.tsv` gives their one-based inclusive positions
(start 201, stop 200 plus consensus length, plus strand). Coordinates come from
construction, never from predicted calls. This is a synthetic genome stand-in;
the consensus FASTA itself is not passed as a genome in the current examples.

The files were generated with the native BBTools fixture producer and Reformat;
the shared reducer produced the graph. The example recipes use two training
iterations and intentionally prohibit endpoint movement, making OFF the correct
boundary-builder result. Do not use them as production training settings or as
evidence of precision, recall or generalization on biological genomes.

From the BBTools installation root, follow the complete commands in:

- `docs/guides/ReduceNcrnaModelsGuide.txt`: 19 models become one while retaining
  support for all 20 synthetic genes. The further protection example restores
  15 models to retain all 19 panel loci; serial regeneration retains those calls
  and the 20 supported genes. `protected.config` supplies every shared input,
  including HBM and truth; `consensus.config` supplies the regeneration recipe.
- `docs/guides/BuildNcrnaBoundaryNetsGuide.txt`: learn and compare networks, retain
  19 matched panel loci with no extra calls, export the correct OFF choice, then
  call the supplied panel using that exported family.

The panel producer is `dev/testbuildncrnaboundarynets.sh mode=examples out=newdir`.
It writes only the two small panel inputs from installed resources; the guides
use the already supplied files, so readers need not run this developer fixture.

Each command writes to a fresh directory. The builder guide creates an absolute
genome-list path at runtime so these example files remain portable.
