# Verification record

Verified 29 October 2026. This checks computation and presentation, not independent
biological validity.

## Numerical checks

- Original upstream R scientific functions ran without edits against pinned commit
  `50587422e6b11259ead6fbc6f867681c788f39b7`; clean exit and asserted report values.
  Full Snakemake HTML/browser-track jobs were not run.
- Separate Python implementation agrees on **every one of 3,852 pair IDs, labels,
  scores, distance values and prediction-presence flags**.
- All nine input/config/license SHA-256 checksums pass. Fresh network retrieval
  and cached/local-checkout routes both successfully rebuild the analysis.
- Nine committed unit/evidence tests pass, including tied scores, weighted AP,
  interval boundaries, original labels, missing-positive retention and R parity.
- Independent reviewer calculations agree with all reported AP values. Sixty
  weighted/tied-score cases were also checked against literal resampling plus
  sklearn AP.
- Figures use the exported data; all three PNGs were visually inspected. No
  synthetic data or untested negatives were introduced.

## Interface checks

- Ramani-shell validator and JavaScript syntax check pass.
- All eight cohort/distance filter combinations were exercised in the browser;
  displayed counts and scores agree with `metrics.csv`.
- All five curated cases pass a DOM smoke test. Browser inspection verifies that
  the NFE2 coverage gap remains experimentally positive and zero-filled, and
  COPZ1 is a significant increase with a nonpositive directional label.
- Actual browser CSV export contains 3,852 data rows, not only the selected gene.
- Desktop 1280px and mobile 390px layouts were inspected. At 390px the page width
  equals the viewport, with no horizontal document overflow. Native labeled
  controls, accessible SVG descriptions and live text results are included.
- Invalid-gene handling and clearing the gene filter were exercised in a separate
  DOM harness. Public URLs are recorded separately after deployment verification.

## Deliberate limits

The result is a saved-example reanalysis, not a claim of independent model
validation, a new regulatory mechanism or personal wet-lab experience. In the
near and far bins the paired bootstrap interval includes zero. Far-bin intervals
omit explicitly counted no-positive resamples. Original pair-bootstrap and added
gene-bootstrap intervals are labeled separately. R dependencies are recorded,
not fully frozen; exact point-result assertions guard against unnoticed drift.
