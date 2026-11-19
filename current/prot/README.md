# Protein tools

## ProkCC 1.8.8

Install `prokcc_v1.8.8.tar` by extracting it into the BBTools `resources/`
directory, then run `prokcc.sh in=bin.fa out=quality.tsv`. The archive creates
one `prokcc/` directory containing its configuration, models and supporting
tables. Keep these resources together; their versioned hashes bind the model
to the inputs used during evaluation.

This release uses a single ten-head composite with 18-bit weights and its
evaluated 32-bit subnet bundle. The bundle contains 132 networks and is gzip
compressed even though its filename ends in `.bbnets`. Model normalization
and prediction interpretation are embedded in the network. Error estimates
remain uncalibrated and are not confidence intervals.

The supplied configuration selects the eight-feature subnet layout and the
tracked-to-called gene ratio. These are model-input contracts, not options to
mix with an older model or bundle. Relative resource paths resolve against the
configuration directory. To use an earlier release, retain its complete
resource directory and pass its configuration explicitly with `config=`.

## Protein substitution matrix interchange

`convertproteinmatrix.sh` reads native version1 TSV or NCBI substitution-matrix
text and writes either representation. Input detection uses file content, not
the extension. `.mat` is recommended for NCBI output. Native parsing retains its
exact version/header rules and rejects unknown metadata; an invalid native file
is never retried as a permissive text matrix.

Export an accepted POSSUM55 native matrix with its actual background source:

```sh
convertproteinmatrix.sh in=possum55.tsv out=possum55.mat outformat=ncbi \
  background=probabilities.tsv band=25to55
```

The background input is the versioned normalized monomer probability table that
produced the matrix, not frequencies inferred from rounded scores. For the
accepted R5 matrix, the retained source is
`Yelan/q46_averaged_r5_derivation_20261115_v1/run/avg10/damped/monomer/probabilities.tsv`
under the project's Dori scratch root. Its matching original matrix is
`run/derive/avg10_damped/matrices/25to55.tsv` in that same derivation root.
The reader checks the directed cell sums and saved marginals before computing
`p[i]=(queryMarginal[i]+referenceMarginal[i])/(2*total)`. Copy or regenerate those
data through the documented project recipe; the converter has no private-path
fallback. The output comments retain the supplied source and all20 probabilities.

POSSUM55 expands to **Prokaryotic Ortholog Sequence SUbstitution Matrix**. It uses
the25–55% identity band and `round(8*log2(observed/expected))`, in eighth-bit score
units, with a linear gap penalty of16 score units. Its native `scale=4` is the
multiplier relative to BLOSUM62's half-bit score units; it is not the8 in that
formula. Conversion never recalculates the400 canonical integer scores.
Before attaching that named formula, the converter verifies that the supplied
`25to55` probability state reproduces all400 scores with those units. This checks
the matrix/background pairing; it does not authenticate an arbitrary file's
biological origin. Use the retained accepted pair for the release resource.
The generic Java `ncbiText` API supplies no derivation claim from a name alone.

NCBI output includes B, Z, X and star. B uses conditional D/N background weights;
Z uses E/Q; X uses the complete background. A pair involving ambiguous symbols
is the expectation over both distributions, rounded once, half away from zero.
This also defines B/B, B/Z, B/X and X/X without intermediate rounding. Star
scores the canonical matrix minimum against every residue symbol, and +1 against
itself. Nonfinite probabilities, nonnormalized total mass and undefined D/N or
E/Q distributions fail loudly.

Our NCBI exports contain a versioned `#bbgap_matrix` comment with name, scale and
gap, so the ordinary matrix loader can read them directly. Plain third-party
NCBI files lack that binding and need explicit metadata:

```sh
convertproteinmatrix.sh in=external.mat out=external.tsv outformat=native \
  name=external scale=1 gap=4
```

Choose those values from the actual scoring contract; the example does not
establish suitable units or gap costs for an arbitrary matrix. Conflicting
metadata is rejected. Native output contains the usual20×20 scores and needs
no background file. Neither command overwrites existing output by default;
`ow=t` permits replacement, but an input alias is always rejected.

The input table must contain all20 canonical row and column labels exactly once.
Reordered axes, spaces/tabs, blank lines and comment lines are supported in NCBI
text. Optional B/Z/X/star/J rows and columns must be complete and contain valid
integers, but their scores are discarded from the runtime binding. Unknown or
duplicate labels, missing cells, overflow and invalid canonical scores fail.

**Runtime scoring is unchanged.** BBTools folds sequence B/Z/J to X; X retains
its scaled BLOSUM62 scoring policy. Internal stops remain invalid, and the
protein-sequence reader strips one conventional terminal FASTA stop. Importing
an NCBI matrix does not activate its exported ambiguity or stop cells. Custom
consensus scoring uses `ProteinSubstitutionMatrix` through `HbmPositionModel`;
the separate `AAAligner` retains its built-in BLOSUM62 behavior.

**Statistical E-value parameters are not provided.** No lambda/K calibration is
inferred from this file, and text interchange does not imply that a particular
BLAST executable accepts a custom matrix or supplies valid statistics for it.

For Java callers, `ProteinSubstitutionMatrix.read(path)` accepts native or
self-described NCBI input. `ProteinSubstitutionMatrixIO.read(path, name, scale,
gap)` supplies explicit metadata for plain NCBI input and checks it against any
embedded metadata. `nativeText(matrix)` and `ncbiText(matrix, background, source)`
return builders for final output; backgrounds use BBTools encoded residue order.
These methods do not mutate the matrix or the supplied background.
