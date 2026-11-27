# Cobalamin riboswitch resources

Enable this optional family with `callgenes.sh in=assembly.fa out=calls.gff cobalamin=t`.
Cobalamin calling remains disabled by default.

The three libraries contain 354 RF00174, 47 RF01482 and 11 RF01689 consensus/HBM
models. They share 2,710 forward 17-mer discovery keys and require two distinct
keys. A merged DDL ranking walk, protected per-key windows and the packaged path
coefficients select candidates for ordinary gene-path competition. The native
identity floors are 0.63 for RF00174 and 0.66 for the other two members.

RF00174 uses anchored covariance-model endpoint refinement with learned corridors,
input-padding overrides and a gathering-threshold hold. RF01482 keeps its compact
K9 positional tables, supported model offsets and separate networks for each end.
RF01689 keeps its native endpoints. Refined ends enter ordinary path selection.
The second-member networks have 69 inputs, 64/48/48 hidden layers and 22 outputs:
21 positions within radius ten and an out-of-range head. Unsupported offsets remain
zero; ties and unavailable neural features retain the applicable current span.

Acceptance is gated by alignment identity. At identity 0.75 or above, the
acceptance check bypasses further covariance scoring; its lower admission floor
is 0.60. Below a member's
native floor only its own covariance model can accept it. At or above that floor,
other-member covariance evidence requires a passing alignment to that member.
An acceptance witness does not relabel the emitting model or replace its endpoints.
Member labels describe the emitting library, not independent biological subtype truth.

The libraries and discovery seeds were built from all 95,128 admitted original
genes in the 40,686-assembly working collection. Ordering profiles and the compact
second-member endpoint resources include the formerly excluded evaluation inputs.
Measurements on that collection after refitting are in-corpus checks, not independent
generalization estimates. Constants were selected on development benchmarks.
The retained covariance corridors have teachers from 96.32% [39,188/40,686] of
assemblies; 1,498 historical panel/evaluation/alias-group exclusions remain.
The padding-word blacklist also
retains its development fit. These are documented best-effort restrictions.

The manifest and caller policy bind one compatible group of libraries, seeds,
profiles, windows, coefficients, endpoint resources and acceptance settings.
`cobalaminroot` selects a replacement group directory. Missing or inconsistent
resources and incompatible individual overrides fail loading. Payload encoding
is declared by the actual manifest paths; metadata and small policies remain raw.
