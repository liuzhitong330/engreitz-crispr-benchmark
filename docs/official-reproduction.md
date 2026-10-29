# Reproducing the saved official example

## Scope and commands

This project ran the **original R merge script and numerical report functions** of
[`EngreitzLab/CRISPR_comparison`](https://github.com/EngreitzLab/CRISPR_comparison)
at commit `50587422e6b11259ead6fbc6f867681c788f39b7`, using the saved inputs in its
`config/config.yml`. This is not merely a Python recreation of the expected result.

From the project root:

```sh
Rscript scripts/install_official.R
Rscript scripts/reproduce_official.R
```

The first command installs missing packages into ignored `.cache/r-library/`,
never into a global library. The second fetches the public repository into
`.cache/CRISPR_comparison/`, checks out the pinned commit, checks for source
modifications, and writes to `results/official/`. If dependencies are missing it
stops with the local install command. Git and a recent R installation are needed;
the completed run used R 4.6.1. The scripts infer the root from their own path, so
they can also be run from another working directory.

The completed run used current local binary dependencies rather than the older
versions in the upstream Conda YAML. `sessionInfo.txt` and `loaded_versions.tsv`
record exact versions. Package installation is not a frozen lockfile: future
installs can select newer versions. Cohort counts and point AUPRC results are
asserted. For strict environment freezing use the recorded versions, or the
upstream Conda workflow. Never commit `.cache/`.

## What is and is not reproduced

The wrapper executes unchanged:

- `mergePredictionsWithExperiment.R`, with a minimal R S4 object supplying the
  same inputs and parameters as Snakemake's example rule;
- `processMergedData`, `calcPRCurves`, `applyCellTypes`, and
  `makePRSummaryTableBS`, as used in the official R Markdown report;
- original numerical bootstrap routines, with the report's 1,000-replicate
  default and a recorded seed of `20261029`.

No upstream source or scientific function body was changed. Plotting helper
function definitions are loaded without unused graphical package imports. The
example's optional-annotation rule only copies its input, so the wrapper reads the
same merged table directly. **The full Snakemake HTML report and browser tracks
were not rendered.** This is an original-code numerical reproduction, not an
identical full workflow/environment execution.

## Reproduced output

| Quantity | Result |
|---|---:|
| Experimental input rows | 5,091 |
| Saved prediction rows | 55,039 |
| TSS-overlapping predicted elements removed | 2,902 |
| TSS-overlapping prediction rows removed | 10,606 |
| Experimental rows absent from gene universe | 4 |
| Valid tested K562 pairs | 3,852 |
| Positive / negative `Regulated` labels | 101 / 3,751 |
| Tested pairs missing ABC predictions (filled 0) | 521 |
| ABC official AUPRC | 0.588235447008955 |
| Distance official AUPRC | 0.401794369883952 |
| ABC precision at configured 0.02 cutoff | 0.489510489510490 |
| ABC recall at configured 0.02 cutoff | 0.693069306930693 |

The original pair-level percentile bootstrap 95% intervals are 0.482064738–
0.679251968 for ABC and 0.293215791–0.493965355 for distance. These are not
gene-cluster intervals. Cluster-aware uncertainty in the extended reanalysis is a
separate analysis.

## Conventions to preserve

- The configured predictor is ABC `powerlaw.Score`, not a newly trained model.
  Overlapping predicted element scores are summed; missing predictions are filled
  with zero **without changing experimental labels**.
- The example filters `ValidConnection == "TRUE"` and keeps missing predictions.
  Only experimentally tested pairs enter the comparison; untested pairs never
  become negatives.
- Official AUPRC uses the upstream yardstick precision–recall endpoint removal
  and `caTools::trapz` convention. It is **not average precision**; report those
  metrics separately.
- Coordinate columns, not the human-readable `name` field, carry the GRCh38
  locations. Some names retain old coordinate strings.
- Upstream merging uses closed IRanges while its promoter filter explicitly
  converts BED-like 0-based intervals. Its distance calculation also follows
  GenomicRanges conventions. The original code is retained for this reproduction;
  a new implementation should audit any convention differences explicitly.
- Gene-universe exclusions are HBBP1 and three PVT1-TSS1 rows; all four were already
  invalid connections and do not change the valid evaluation cohort.
- This is a reanalysis of provided predictions and evidence, **not independent
  validation, new experimental evidence, or proof of generalization**.

The original code emitted a dplyr `case_match()` deprecation warning under dplyr
1.2.1; no workaround or source patch was necessary.

## Saved artifacts

`results/official/` contains original merged rows, the valid included-predictor
table, ROCR PR points, performance summary, counts, filter diagnostics, and
session/provenance records. These numerical outputs let the extended analysis
cross-check its pair selection and scoring against original-code results.
