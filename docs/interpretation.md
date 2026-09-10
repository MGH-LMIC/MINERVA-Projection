# Reading the results

The framework asks two related questions: does a published disease pattern align with case status, and which individual microbial relations support or oppose that pattern?

## A participant score summarizes a pattern

Each genus contributes **+1** when its observed abundance state matches the literature expectation, **−1** when it opposes that expectation, and **0** when the state is neutral or there is no signed expectation. Summing these contributions gives a participant's raw literature-concordance score for a disease.

A positive contribution can come from **low abundance**. If the literature expects a genus to be less abundant in the disease and the participant is classified Low, the two directions align: the contribution is +1.

Case–control association is evaluated within studies using standardized scores. An association with case status describes a group-level pattern. Participants can overlap substantially in their scores.

## A relation shift describes the group contrast

For one disease and genus:

1. Assign Low (−1), neutral (0) or High (+1) to each participant using the external healthy reference.
2. In each study, subtract the control group's mean state from the case group's mean state.
3. Average these differences with equal weight for each contributing study.
4. Compare the resulting sign with the signed literature expectation.

For example, suppose the mean state is −0.20 in cases and +0.10 in controls. The case-minus-control shift is −0.30. It agrees with a literature expectation of lower abundance and opposes an expectation of higher abundance.

Many observations may remain inside the healthy range. The group shift depends on the **balance of Low and High states**, so a difference can still arise from the observations outside that range. A positive shift can reflect more High states, fewer Low states, or both.

## Agreement and opposition guide different follow-up questions

| Finding | A useful next question |
| --- | --- |
| Agreement supported across studies | Which relation merits a mechanistic experiment or biomarker evaluation? |
| Opposition to the literature expectation | Does the source evidence describe the same comparison, and does the association vary by population or context? |
| Exactly zero group shift | Is there no net state difference under this reference rule, despite variation among individuals? |
| No signed literature direction | What evidence would be needed to evaluate this relation? |

Agreement describes correspondence under this projection. It does not establish causality, independent replication of every source publication, or the efficacy and safety of an intervention. Relation-level direction summaries are descriptive; their signs alone are not statistical significance tests.

## Keep the denominators visible

When summarizing agreement, report the number of **linked relations** used in the denominator. An unlinked relation and a linked relation with exactly zero shift both make no contribution to an alignment summary, but they represent different evidence states. Keep them distinguishable in exported tables.

Compare the full set of measured genera with any selected display subset separately. Selecting genera with larger shifts can change the proportions of agreement, opposition and zero shifts.

[Back to the README](../README.md)
