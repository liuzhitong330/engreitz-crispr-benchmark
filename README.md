# Enhancer–gene predictions, checked against perturbations

A focused K562 reanalysis of the saved **ABC power-law** predictions in the
[Engreitz Lab CRISPR benchmark](https://github.com/EngreitzLab/CRISPR_comparison).
Prepared by Cathy Liu with AI-assisted implementation and independent code/result
checks. Not an official lab product, a new predictor, or independent biological validation.

**Question:** on exactly the same experimentally tested enhancer–gene pairs, how
does saved ABC compare with ranking by distance to the target gene, and which
disagreements deserve a closer look?

## Result in brief

On **3,852 valid tested pairs (101 down-regulation calls; 54 genes)**, average
precision is **0.602 for ABC versus 0.414 for distance**. The original R numerical
report is also reproduced: its different, trapezoidal AUPRC convention gives
**0.588 versus 0.402**. These are two differently defined metrics, not contradictory runs.

ABC ranks better within all three supplied distance ranges in this example.
However, coverage is incomplete: **521 pairs have no matched prediction, including
seven positive experiments**. The main comparison preserves the official zero-fill
policy; a clearly conditional covered-only comparison uses the same 3,331 pairs
for both methods. No untested pairs become negatives.

- [Interactive report](https://liuzhitong330.github.io/engreitz-crispr-benchmark/)
- [Results, limitations and follow-up](RESULTS.md)
- [Plain-language walkthrough](docs/walkthrough.md)
- [Exact original-code numerical reproduction](docs/official-reproduction.md)

## Run

Python 3.12 was used. From the repository root:

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/reproduce.py
python -m unittest discover -s tests -v
python -m http.server 8000
```

Open `http://localhost:8000`. The analysis fetches nine small public files from
commit `50587422e6b11259ead6fbc6f867681c788f39b7` and verifies SHA-256 checksums.
It rebuilds the cohort, metrics, two cohort-policy comparisons, 2,000-draw paired
gene-cluster bootstrap intervals, three figures and browser data. It asserts
agreement with the checked-in original-R output for **every pair's label, score,
distance and prediction coverage**. It does not call a remote AI/model service.

To rerun the original R routines as well (R, Git and network access required):

```sh
Rscript scripts/install_official.R
Rscript scripts/reproduce_official.R
python scripts/reproduce.py
```

The R wrapper executes unmodified upstream merge/numerical report functions;
it does not render the complete Snakemake HTML report. R package versions and
the original report's 1,000-draw **pair-level** bootstrap are recorded separately.
The installer resolves current compatible packages; it is not a frozen R environment.
Python dependencies are pinned; `requirements-lock.txt` records the full tested set.

## What is included

| File | Purpose |
| --- | --- |
| `scripts/reproduce.py` | Auditable cohort construction, ranking comparison, uncertainty, figures |
| `data/input_manifest.json` | Source commit, file paths and input SHA-256 checksums |
| `results/cohort.csv` | Every evaluated pair, original label/effect/p-value/power and coverage |
| `results/excluded_rows.csv` | Original excluded rows and reasons |
| `results/metrics.csv` | AP, prevalence, sample counts and gene-bootstrap intervals |
| `results/examples.csv` | Five reproducibly selected illustrative cases, with selection rules |
| `results/qc.json` | Coordinate, identifier, cell-type and missing-prediction audit |
| `results/official/` | Original R numerical outputs and environment/provenance |
| `figures/01-pr.*`, `02-distance.*`, `03-disagreement.*` | Three useful figures (SVG and PNG) |

## Boundaries

The source filename mentions Fulco2019, but contains multiple studies. We retain
the supplied example—not the much larger K562 training file. Coordinate columns
are GRCh38; embedded coordinates in row names are historical and must not be
parsed. `Regulated=FALSE` means no significant **negative** effect under the
source definition: 31 such rows are significant increases. Ensembl IDs are absent,
so matching uses checked exact gene symbols, not inferred ID conversions.

This benchmark is historically related to ABC development. A new bootstrap or
distance stratification does not make it independent validation. Current lab
resources also include newer rE2G/scE2G models; this project does not compare them.
See [third-party attribution](THIRD_PARTY_NOTICES.md) for sources and licensing.
