# pacbio package

`current/pacbio` contains older PacBio assembly, adapter-removal,
site-stacking and synthetic-data utilities. Several classes are legacy pipeline
entry points or numbered variants retained for compatibility; keep behavior
changes separate from documentation and formatting cleanup.

## Package map

| Class | Role and current limits |
| --- | --- |
| `CalcCoverageFromSites` | Reads stacked site records, trims site ends, reports whole-chromosome per-base coverage and correctness summaries, and optionally writes 32-bit total-depth coverage arrays. |
| `GenerateMultiChrom` | Builds synthetic chromosome copies from an input genome build. The command line requires source/copies/build, optionally followed by all three N-gap settings; mutation rates are fixed constants in source. |
| `MakePacBioScript` | Generates legacy PacBio correction shell scripts from parsed input lists and resource settings. |
| `MergeFastaContigs` | Concatenates FASTA/FASTQ contigs into artificial chromosomes with N padding and contig bookkeeping. |
| `MergeReadsAndGenome` | Appends generated genome reads and optional input reads to one output stream. It renumbers retained reads and trims terminal reference-N sequence for `chrom>=1`; it does not perform generic quality correction, and it now propagates input and output close errors. |
| `PartitionFastaFile` | Splits a FASTA stream into length-limited output files. Its positional arguments leave an unused gap before the optional total-output cutoff, and the static total-output cutoff can stop mid-record. |
| `PartitionReads` | Distributes single reads or paired roots across numbered output streams for partitioned downstream work. It poisons writers without joining or checking their error state. |
| `ProcessStackedSitesNormalized` | Normalizes and filters stacked alignment sites across intervals. |
| `RemoveAdapters2` | Legacy threaded PacBio adapter finder/remover. Read splitting is unpaired-only; paired input must use `split=f`. |
| `RemoveAdapters3` | Legacy alternate adapter remover with flanking reverse-complement validation of weak hits. Read splitting is unpaired-only; paired input must use `split=f`. |
| `SiteR` | Compact site record for the stacking/sorting pipeline. It stores chromosome, strand, numeric read id and mate bit; mate 0 is `numericID`, positive mate-1 IDs keep their historical `-numericID` keys, and mate-1 ID 0 is `Long.MIN_VALUE`. It has typed equality helpers for `SiteScore` and `SiteScoreR`, not ordering/comparison methods. |
| `SortSites` | Sorts temporary stacked-site records by genomic position or read id. |
| `SplitOffPerfectContigs` | Uses coverage-based classification to emit high-coverage ranges and write retained lower-coverage ranges into chromosome-array files. If `buildout` is omitted it defaults to the input build and may replace that build's chromosome files and delete its index files. |
| `StackSites` | Legacy site-stacking driver. |
| `StackSites2` | Second legacy site-stacking driver; kept alongside `StackSites`. |

## Shared maintenance notes

- Most package tools predate current BBTools entry-point conventions and are
  used through generated or historical pipelines rather than direct shell
  launchers.
- `SiteR` key encoding is shared with the `var.GenerateVarlets*` consumers;
  preserve the nonnegative numeric id and low-bit mate convention when touching
  that path.
- The numbered classes are intentional historical variants. Do not merge or
  delete them as part of ordinary cleanup.
- Several tools write BBTools genome-build chromosome arrays directly. Review
  build arguments before running them against a real build, and keep generated
  controls and output paths outside Brian's Eclipse tree.
