# Legacy neural-configuration test model

`selected.bbnet` is the frozen model required by `test.assemble.TadpoleNeuralConfigTest`:
8,456 bytes, sha80 `05d0812fb6b4cbf1c5fe`. It is a compatibility fixture for model
binding and configuration checks, not a new model recommendation or accuracy claim.

Fischl authorized copying the model and its two original provenance tables from
`ecc-nn/results/plan19_architecture_v1/` on November 8, 2026. The model's current
source hash was checked and all three copied files were byte-compared. The original
producer is recorded in that project's scripts 90–94 and decision D106.

`artifact_sha80.tsv` lists the original experiment's artifacts; most of those
artifacts are not included here. `selection.tsv` records the original selection
rule and explicitly says that its external set had become development data and
fresh external validation was required. These tables do not establish a new
independent evaluation.

The test takes a BBTools fixture root and this model path as separate positional
arguments. It also requires the bundled PacBio models, creates temporary test
networks/configuration, and writes `pacbio_smoke.fa` beneath its supplied root.
Use an isolated test invocation and preserve any existing file at that output path.
