# A plain-language walkthrough

## 1. What biological question are we asking?

An enhancer is a piece of DNA that can help control a gene. The nearest gene is
not necessarily the gene it regulates. Engreitz Lab studies these connections
using computational predictions and perturbation experiments.

Here we ask a narrow question: **among pairs that were actually tested in K562
cells, does the saved ABC score put the experimentally supported connections near
the top of the list better than simply choosing the closest ones?** We also ask
which cases disagree and whether a missing prediction explains some apparent misses.

We are not testing whether every possible enhancer–gene pair is connected. Most
possible pairs were never tested and must not be called negatives.

## 2. What does each input represent?

**The experiment table:** one row is a tested DNA region plus a target gene.
It records genomic coordinates, the measured expression effect, significance,
adjusted p-value, estimated power and study source. The source's `Regulated` label
is true for a significant negative expression effect. Silencing a genuine
activating element can reduce its target gene's expression, which motivates this
label. A false label does not always mean no effect: some genes increased instead.

**The prediction table:** candidate regions, genes and already-computed ABC
`powerlaw.Score` values. ABC combines regulatory activity with an estimate of
element–gene contact and compares candidate elements for a gene. This saved example
uses a contact estimate based on genomic distance—not newly measured Hi-C contacts.
We did not train ABC or calculate new activity measurements.

**The gene annotation:** the supplied locations of transcription start sites
(TSSs), where transcription begins. We use these to compute the distance baseline.

**The configuration and mapping files:** the original rules for matching K562
predictions to K562 experiments, summing overlapping predicted elements, filling
missing scores and marking the fixed ABC threshold. These make the comparison
reproducible rather than a collection of improvised choices.

## 3. What did we actually do?

1. Pinned the lab repository to one exact version and recorded file checksums.
2. Ran its original R merge and numerical evaluation functions unchanged.
3. Independently reconstructed the cohort in Python and checked every pair's
   ID, label, ABC score, distance and coverage against the R output.
4. Ranked the same pairs twice: high ABC first; short genomic distance first.
5. Compared the rankings overall, within three distance ranges and after excluding
   unmatched predictions from both methods. No score was trained or threshold tuned.
6. Inspected five reproducibly chosen agreements/disagreements in their original
   assay context and made the pipeline, tables and figures available to rerun.

## 4. Why do these analysis choices matter?

**Coordinates:** the readable row names contain old coordinates. Matching on
those names as though they were current locations would be wrong. We instead use
the explicit GRCh38 columns and the supplied gene annotation.

**Same tested pairs:** comparing ABC only where it has predictions with distance
on every experiment would give the methods different tests. The main result uses
every eligible tested pair for both methods. The matched-only check again gives
both methods exactly the same, smaller test.

**Missing is not negative:** a missing ABC entry receives score zero under the
official evaluation rule, but its experimental label stays unchanged. Seven such
missing pairs have a positive experimental label. This distinction is visible.

**Direction of effect:** `Regulated=FALSE` means the source did not call a
significant decrease. There are 31 significant increases in this group, so “ABC
predicted an effect where none exists” would be too strong.

**Distance ranges:** close pairs have a much higher positive rate than distant
pairs. Compare ABC with distance inside each range; do not interpret lower AP
in a difficult, low-prevalence range as though everything else were equal.

**Uncertainty:** many pairs share the same gene. We resample whole genes, using
the same resample for both methods, instead of pretending every pair is independent.
This is still only a description of this selected dataset. In the far bin some
resamples contain no positives; their AP is undefined and they are counted explicitly.

## 5. How do I read the main numbers?

Imagine reading down a ranked list of pairs. **Precision** asks, “Of the pairs
seen so far, what fraction had a down-regulation call?” **Recall** asks, “What
fraction of all 101 down-regulation calls have I recovered?”

**Average precision** summarizes precision as additional true responses enter
the ranking. Higher is better on the same test set. Here it is **0.602 for ABC
and 0.414 for distance**. This is not 60.2% versus 41.4% accuracy, and neither
number is a probability that one specific pair regulates a gene.

The official R report gives a different pair of values, **0.588 and 0.402**, because
it uses a different area calculation and endpoint convention. Both are reproduced
and clearly labeled. Do not mix them in a sentence as if they were the same metric.

ABC's point estimate is higher in all three distance bins. However, the near
and far-bin uncertainty intervals include no advantage, and the far bin has only
13 positives. “ABC definitely wins at every distance” is not justified.

## 6. What do the example cases teach me?

- **MYC:** a supported connection about 163 kb away still has an ABC score above
  the supplied threshold. Distance is useful but not sufficient.
- **FTL:** an experimentally supported nearby pair falls just below the threshold.
  Threshold decisions can hide useful ranking information.
- **NFE2:** the experiment is positive but no matching saved prediction exists.
  Check coverage before trying to explain a model's biological reasoning.
- **DNASE2:** a high score does not give a significant down-call in this assay.
  This is a useful follow-up candidate, not proof that either method is wrong.
- **COPZ1:** a significant increase is scored as nonpositive because this benchmark
  specifically asks about decreases. Always check what a label means.

## 7. What can I accurately say to the lab?

> I reproduced the saved K562 example from your CRISPR comparison repository and
> built a small explorer comparing ABC with a distance baseline on the same tested
> pairs. ABC ranked down-regulating pairs better overall in this example; the tool
> makes the distance-specific uncertainty, missing predictions and individual
> prediction–experiment disagreements easy to inspect.

If asked about ownership, be straightforward: **“I developed this with AI-assisted
coding support. I can walk through the cohort, matching rules, checks and limits.”**
Only say you can independently implement a method after you have actually learned
and practiced it. The source notebooks/code and this walkthrough are a starting
point for that preparation, not a substitute for it.

Avoid: “I independently validated ABC,” “I discovered a new enhancer,” “ABC works
better in every cell type,” “the negative experiments prove no regulation,” or
“I trained a new AI model.” None of those was done here.

**A good conversation question:** “Would the same distance-matched comparison be
useful on a genuinely separate perturbation dataset, and which assay or coverage
differences would you want controlled before interpreting the disagreements?”
