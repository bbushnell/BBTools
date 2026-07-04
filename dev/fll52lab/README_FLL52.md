# FLL52 — production integration notes
*Amber, 2026-07-04. Research log: /mnt/c/playground/Amber/plans/tail_mle_research.plan*

## What shipped (staged to C:\workspace\BBTools, awaiting your compile)
- `cardinality/Fll52.java` — the species: 16-bit words = 4-bit localExp + 5x2-bit
  antennae + 2-bit saturating counter of future-bit re-hits (reset on promotion).
  512 bytes = 256 words = 1280 tails, every bit accounted for.
  `cardinality()` = composite-likelihood estimate blended with a complexity
  correction from the modular window net: w=clamp((cHat-0.85)/0.10),
  nHat = w*MLE + (1-w)*added*cHat.  Graceful fallbacks if resources are absent
  (warns once, degrades to moment estimate).
- `cardinality/Fll52MLE.java` — (globalExp, phase)-conditioned composite
  likelihood; auto-loads the nearest-size table from resources on first use.
- `cardinality/Fll52Calibrate.java` — regeneration tool for the tables.
- `cardinality/CardinalityTracker.java` — +3 lines each at the three dispatch
  sites: type "FLL52"/"Fll52".  (Diff vs your copy was clean before staging.)
- `ml/RegressionTrainer.java` — continuous-output trainer emitting real BBNet
  files (CellNet's own serializer; self-verifying round-trip on every save).
  Usage: java ml.RegressionTrainer in=data.tsv out=net.bbnet dims=16,32,1
  [final=rslog|sig] [epochs] [batch] [lr] [wd] [seed].  balance is not needed
  (regression-native).  Measured: 1.75x lower error than tuned ml.Trainer on
  the complexity regression task.

## Resources (staged to C:\releases\bbmap\resources)
- `fll52net.tsv.gz` — the modular window net (16 features -> 32 tanh -> linear).
  One net serves ALL estimator sizes: k = numWords/256 predictions averaged.
- `fll52mle_512w.tsv.gz` — likelihood table for 512-word (1KB) sketches.
- PENDING: fll52mle_{256,1024,2048}w.tsv.gz from cluster job 23600754 (Amber
  will stage on arrival).  Fll52MLE.forWords() picks the nearest size.

## Measured (honest equal bytes, WidthWtAbsErr, 128 instances, 1KB)
```
              HC        LC(minrand iter4)
AVLL_HLDLC   1.584%    1.696%
FLL52_MLE    0.931%   22.9%     (likelihood alone, duplicate-blind)
FLL52_blend  1.098%    1.868%   (the production cardinality() path)
```
Complexity readout accuracy (held-out shapes): ~2-4% at 512B, improving ~1/sqrt(k).

## Regeneration (Dori, 128 cores, ~15 min each)
```
java -Xmx64g -cp . cardinality.Fll52Calibrate <words*5> 8192 4096 128 fll52mle_<words>w.tsv
java -Xmx256g -cp . cardinality.Fll52FarmW 800000 64 128 32 netwin52.txt   # window net
```

## In this lab folder (experimental, NOT staged to the workspace)
- Fll53.java / Fll53FarmW.java — 21-bit organism: 5x3-bit antennae + 1 counter.
- Fll53b.java / Fll53bFarmW.java — 20-bit: 4x3-bit antennae + 2 counters with
  Brian's flow-on-promotion vs reset A/B (Fll53b.FLOW toggle).
- ComplexityFarm[W].java, FllComplexity.java, FllSkewTest.java, FarmEval.java,
  WidthWt52.java — the training farms and referee harnesses.
Gates pending (cluster jobs 23600754/23600843/23600950); winners get the same
production treatment.

## Known limits
- Sketches under 256 words (512B) skip the complexity blend (net window
  doesn't fit); they run pure likelihood = duplicate-blind.
- `add(CardinalityTracker)` (merge) is unimplemented — counters merge is
  ill-defined across windows; needs design if union semantics are wanted.
- Tables/net are trained at k=31 hashing defaults; recalibrate if that changes.
