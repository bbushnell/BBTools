# Cobalamin riboswitch resources

`cobalamin=t` enables this opt-in family. Calling remains disabled by default.
The selected configuration uses 1,928 shared forward 17-mer keys, at least two
distinct seed keys per window, K9 model ranking and a shortlist of up to 100
models. The three member libraries contain 508 RF00174, 39 RF01482 and 16 RF01689
consensus/HBM models. Compressed HBM resources are distributed; HBM verification
is disabled with its independent threshold set to 1.01.

Each member's `cutoffs.tsv` sets equal pass and borderline identity thresholds:

| Member | Identity |
|---|---:|
| RF00174 | 0.64 |
| RF01482 | 0.66 |
| RF01689 | 0.66 |

The group `manifest.tsv` binds the shared seeds and each member's consensus,
compressed models and cutoff table by path, byte count and sha80. These are
scalar member cutoffs. The selected package contains no per-model length or
identity profile.

On the reviewed 300-genome engineering panel, the selected production
configuration matched 366 of 422 supported reference loci with 18 unmatched
calls: 95.31% precision and 86.73% recall. A further 23 calls overlapped unresolved
reference evidence and remained unscored. This all-corpus construction benchmark
is not an independent holdout. Both endpoints were exact for 19 matched loci.
All other RNA calls were unchanged. Five annotated CDS starts improved, and five
predictions without matching source CDS annotations were omitted. Absence of an
annotation does not establish absence of biological function.

Measured process CPU increased by about 24% at both one and sixteen threads,
which exceeded the selected 5% default-on limit. The 50-model alternative lost
five true positives (1.18 percentage points of recall), so the shortlist remained
100. Resource construction identities and caller identities are different
measurements; these cutoffs came from caller validation.

`cobalaminroot` selects another group directory. Explicit `ncrnaidpass` and
`ncrnaidborderline` overrides take precedence for the selected `ncrnafamily`.
Omit those identity overrides to use the packaged member values. A declared
missing or invalid resource fails loading. The canonical construction, benchmark,
cost and decision records are in the ncRNA project's `rrnafallback` repository.
