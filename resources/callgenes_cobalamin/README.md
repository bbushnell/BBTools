# Cobalamin riboswitch resources

This package provides cobalamin discovery models and boundary refinement.
`cobalamin=t` enables the opt-in family; calling remains disabled by default.

The fixed discovery libraries contain 354 RF00174, 47 RF01482 and 11 RF01689
consensus/HBM models and share forward 17-mer seeds. Their identity cutoffs are
0.64 for RF00174 and 0.66 for RF01482/RF01689. This package does not change other
families or choose a new shortlist policy. The manifest binds the shared seeds,
consensus sequences, compressed HBMs, cutoffs and positional configurations.

All three libraries and the shared seeds were built from all 95,128 admitted
genes in the 40,686-genome working roster, with no genomes held out.
Model reduction used those full gene sets; caller-protection checks used
1,000-genome and 10,000-genome panels, followed by measured model add-backs.
These panels were not excluded from resource construction.

RF00174 and RF01482 use supported per-model endpoint offsets followed by separate
first-stage networks for each end. The offset support floor is ten matches;
unsupported models retain zero applied offsets. Each network has 69 inputs,
64/48/48 hidden layers and 22 outputs: 21 positions at radius ten plus an
out-of-range head. Three K9 tables per end, span geometry and selected-model HBM features
retain the same training and serving definitions.

The packaged policy is A at both ends, with zero positional and head margins.
The best in-range move must strictly beat the current position, other positions
and the out-of-range head; ties hold. A valid model-centered span is retained
when the head wins or neural features are unavailable. Invalid centering retains
the original span. Every call is offered to refinement, subject to the existing
physical-end, feature and span guards. Refined ends enter ordinary path selection.
There is no deferred reporting or second-stage gate. RF01689 has no positional
configuration and receives no direct boundary refinement from this package.

Offsets use all 79,863 recorded matched witnesses from the 40,686-genome working
roster. Tables and networks include all 91,072 main-member and 2,270 second-member
original teaching loci; none of the former evaluation genomes is withheld.
Initialization 58 was fixed before training. Same-locus monitor draws are not
independent validation. The 1,000-genome and historical 300-genome runs are
in-sample packaging checks with different reference definitions, not accuracy or
generalization estimates. Their results belong in the accompanying project
record, including gains, losses, exact/within-three and endpoint-distance changes.

Weighted tables and HBMs use BGZF compression. Each compressed payload is checked
against its original bytes after decompression. `cobalaminroot` selects a group
directory; missing or inconsistent declared resources fail loading. The integrated
release passed the full build, 33 finite checks and unchanged-default whole-GFF
regressions on 203 and 223 genomes. Caller timing and reference comparisons on
10,000 genomes are recorded with their measured gains, losses and limitations;
the candidate's endpoint results are in-sample, not generalization estimates.
