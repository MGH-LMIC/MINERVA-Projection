# Input and output guide

The command-line workflow accepts abundance profiles, participant metadata and signed literature extraction records. Prepare one consistent set of genus and disease labels across the files.

## 1. Participant metadata

`participants.csv` contains these columns:

| Column | Meaning |
| --- | --- |
| `sample_id` | Unique identifier matching an abundance row |
| `source` | Profiling resource; references are constructed within this resource |
| `study_id` | Study identifier; controls and cases are compared within this study |
| `group` | `case` or `control` |
| `disease` | Disease label for a case; empty for a healthy control |

Include the target case–control studies and other healthy-control studies needed for the external reference. Use a common `source` only for profiles intended to share a reference scale and processing resource.

Prepare one selected profile per participant before running the analysis. For example, resolve repeat sampling and technical replicates in the input cohort. The identifiers should describe that prepared cohort.

Healthy controls within a study can contribute to more than one disease comparison. A comparison-record count can therefore differ from the number of unique participants.

## 2. Abundance matrix

`abundance.csv` has a `sample_id` column followed by one numeric column per genus. Use the same genus labels as the literature file wherever evidence is available. Genera without a signed literature direction can remain in the matrix; they receive zero literature contribution and remain distinguishable as unlinked relations.

- Supply **percentage points**, from 0 to 100: an abundance of 2% is entered as `2`, not `0.02`.
- Each row must sum to at most 100. A subset of genera may sum to less than 100; retain their original abundance scale rather than renormalizing that subset.
- Values must be finite and non-negative. Missing abundance values are not accepted.
- A measured zero means undetected and receives a neutral state under the primary rule.

If your matrix uses fractions summing to one, convert it to percentage points before running the workflow.

## 3. Literature extraction records

`votes.csv` contains one row per extraction record:

| Column | Meaning |
| --- | --- |
| `record_id` | Unique extraction identifier; duplicated identifiers are rejected |
| `disease` | Disease label matching the participant metadata |
| `genus` | Genus label matching an abundance column |
| `publication_id` | Identifier that groups records from the same publication |
| `direction` | `+1` for higher expected abundance, `−1` for lower expected abundance, or `0` for a neutral record; use numeric CSV values `1`, `-1`, `0` |

Records from one publication are summed within a disease–genus relation and reduced to one sign. Those publication votes are then summed and reduced to the relation's consensus sign. A publication therefore contributes one vote per relation, even if several extraction records describe it.

Check that each direction is appropriate for the intended disease–abundance comparison. Keep the original supporting evidence in your own source records so that an unexpected result can be examined in context.

### Optional publication exclusions

To exclude publications associated with your evaluation or reference cohorts, provide a CSV with a `publication_id` column. Exclusions are applied before consensus is calculated.

```bash
python -m minerva_projection run --participants inputs/participants.csv --abundance inputs/abundance.csv --votes inputs/votes.csv --exclude-publications inputs/exclusions.csv --output results
```

Publication identifiers must match those in `votes.csv`. The exclusion list is supplied by the user; it should reflect the cohort and evidence provenance available for the intended analysis.

## Run the analysis

```bash
python -m minerva_projection run --participants inputs/participants.csv --abundance inputs/abundance.csv --votes inputs/votes.csv --output results
```

Optional thresholds:

| Option | Default | Meaning |
| --- | --- | --- |
| `--min-group-size` | `10` | Minimum cases and controls from the same study in each study–disease task |
| `--min-reference-size` | `100` | Minimum external healthy-reference pool from the same source, excluding the target study |
| `--display-genera` | `28` | Maximum number of selected genera used to determine the graded-alignment display scale |
| `--min-landscape-sources` | `2` | Minimum profiling resources supporting a disease for relation-landscape summaries |

The reference-pool minimum applies before restricting the reference to positive observations of each genus. Inspect genus-specific support in `references.csv`; a large overall pool does not imply strong support for every genus.

Choose a new output directory, or an empty existing directory. A non-empty output directory is rejected to protect previous results.

## Inspect the outputs

Start with `run_summary.json` and `task_eligibility.csv`, then follow the workflow through these tables:

| Output | What to inspect |
| --- | --- |
| `publication_votes.csv` | One vote for each publication and disease–genus relation |
| `consensus.csv` | Combined literature expectation for each relation |
| `task_eligibility.csv` | Which study–disease comparisons satisfy the analysis criteria |
| `references.csv` | Healthy-reference bounds and support for each target study and genus |
| `states.csv` | Observed Low, neutral and High states |
| `contributions.csv` | Signed literature–observation contributions |
| `scores.csv` | Raw and within-task standardized participant scores |
| `disease_associations.csv` | Study-stratified score associations with case status |
| `study_relation_shifts.csv` | Case-minus-control mean state differences within studies |
| `relation_alignment.csv` | Relation-level comparisons of cohort shifts and literature expectations |
| `disease_alignment_summary.csv` | Disease-wise summaries over linked relations |
| `genus_alignment_summary.csv` | Genus-wise summaries over linked disease relations |
| `study_directional_alignment.csv` | Directional agreement summarized within each study–disease task |
| `disease_directional_alignment.csv` | Equal-study mean of the study-level directional-agreement summaries |
| `run_summary.json` | Execution and analysis summary |

The three relation-landscape tables (`relation_alignment.csv`, `disease_alignment_summary.csv` and `genus_alignment_summary.csv`) are written when at least one disease satisfies `--min-landscape-sources`. With a single profiling resource, the default run still exports eligible score associations and study-level shifts. Set `--min-landscape-sources 1` explicitly if you also want an exploratory relation landscape for that resource.

If no task is eligible, the command exits with status 2 and writes the eligibility report and run summary; it does not create scientific estimates. The `reasons` column explains exclusions, including insufficient group size, insufficient external controls or no signed relation on the supplied genus axis.

[Back to the README](../README.md)
