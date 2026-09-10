# How the projection works

The algorithm puts literature expectations and observed abundance on a common directional scale. It retains both the participant-level contributions and the study-level contrasts so that an overall disease pattern can be examined relation by relation.

## 1. Resolve a literature expectation

For each disease–genus relation, sum the signed extraction records within each publication and take the sign. This produces one publication vote: −1, 0 or +1. Sum the publication votes with equal weight and take the sign again to obtain the relation's expected direction, **K**.

- **K = +1:** higher abundance is expected in the disease.
- **K = −1:** lower abundance is expected in the disease.
- **K = 0:** no net signed expectation from the available votes.

Apply user-supplied publication exclusions before these operations. A publication vote represents a publication; it need not represent an independent biological cohort.

## 2. Estimate an external healthy-reference interval

For a target study and genus, collect healthy controls from **other studies in the same profiling resource**. Restrict that genus's reference distribution to detectable, positive abundance values. Its 5th and 95th percentiles define the lower and upper bounds, **L** and **U**.

The same bounds are applied to cases and controls in the target study. The controls defining the reference interval are therefore distinct from the controls used in its case–control comparison.

The default minimum of 100 controls applies to the complete external pool. Genus-specific positive reference support is reported separately. If a genus has no positive reference observation, the fallback bounds are L = 0 and U = 10⁻⁶ percentage points.

## 3. Classify the observed abundance

For a positive observed abundance **x**:

| Rule | State S |
| --- | --- |
| x ≤ L | Low, **−1** |
| x ≥ U, after the Low rule | High, **+1** |
| Otherwise | Neutral, **0** |

An undetected genus (x = 0) is neutral. Low takes precedence when the two reference bounds coincide. Consequently, boundary values and ties have deterministic assignments.

The neutral state means that an observation makes no directional contribution under this rule. It does not imply absence of biological variation or that every value in the interval is clinically healthy.

## 4. Combine the two directions

The contribution of a genus is the product **S × K**:

- +1 if the observed state and literature direction agree;
- −1 if they oppose one another;
- 0 for a neutral observation or no signed literature expectation.

In particular, **Low × lower expected abundance = (−1) × (−1) = +1**.

Summing contributions across genera gives a participant's raw literature-concordance score for the disease. For disease-level association analysis, standardize scores across cases and controls from the same study within each study–disease task, using the sample standard deviation (ddof = 1).

## 5. Examine disease-level score associations

Study-stratified conditional logistic regression relates the standardized score to case status. The odds ratio describes the association for a one-within-study-standard-deviation increase in score. It estimates a common within-study association across the contributing studies.

The exported results include 95% confidence intervals, two-sided Wald P values and Benjamini–Hochberg adjustment across evaluated diseases. A non-estimable disease contributes P = 1 to the adjustment family while its raw P value remains missing. Constant-score tasks do not provide information about a score association. Inspect the analysis status and contributing-task information before interpreting an estimate.

A fit that reports non-convergence is marked `not_estimable_nonconverged` and its estimate and raw P value are left missing, even if the optimizer returned a finite coefficient.

## 6. Examine individual relations

For each study, disease and genus, subtract the controls' mean abundance state from the cases' mean abundance state. Average these shifts with **equal weight for each study**. A larger study therefore does not automatically determine a disease's relation profile. By default, relation-landscape summaries include diseases supported by at least two profiling resources; this threshold is configurable.

Compare the sign of that average shift with K to classify a linked relation as agreement, opposition or exactly zero shift. Keep unlinked relations separate and exclude them from linked-relation denominators.

This is a comparison of reference-defined states. Its direction can differ from a comparison of raw mean or median abundance. Many participants may have neutral states while differing proportions of Low and High observations produce a non-zero group shift.

### Graded alignment and categorical summaries

For an optional display scale, rank genera by the sum of their percentile ranks for total absolute cohort shift and number of signed disease links. Select up to 28 genera by default. Let **c** be the 95th percentile of absolute shifts in those selected cells, with a floor of 10⁻⁸. The graded alignment is **clip(A / c, −1, +1) × K**, where A is the equal-study mean shift. The full supplied genus axis remains in the output; `display_selected` identifies the selected subset.

Categorical agreement uses the sign of the raw shift, so changing this display cap cannot turn agreement into opposition. A linked relation with A = 0 is `exact_zero`. A relation with K = 0 is `unlinked` and has missing graded alignment. Disease and genus summaries count aligned, opposed and exactly zero relations among linked cells, and report coverage separately.

### Study-level directional agreement

A second summary first computes **sign(A_study) × K** for each linked relation within each study. It averages those signs within a study, including exactly zero shifts, then averages studies equally for the disease. This measures directional agreement across studies; it can differ from classifying the sign of the shift after studies have been averaged. An explicitly named `nonzero_only` summary also excludes zero shifts. These exports do not supply bootstrap confidence intervals.

## Two summaries answer different questions

**Disease-level score association** asks whether the combined literature-defined pattern is associated with case status within studies.

**Relation-level agreement** asks which individual genus contrasts follow or oppose the literature direction. A positive overall association can coexist with substantial disagreement among individual relations.

Use the [interpretation guide](interpretation.md) to connect these summaries to possible follow-up studies.

[Back to the README](../README.md)
