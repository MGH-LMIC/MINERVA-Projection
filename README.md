# MINERVA-Projection

**Literature-to-cohort projection for microbiome discovery.**

What happens when published microbiome knowledge is projected onto human cohort data? MINERVA-Projection turns literature directions into explicit expectations, compares them with observed abundance states, and shows where the two agree or point in opposite directions.

Use it to examine a disease pattern as a whole and the individual microbial relations within it. A strong overall pattern can contain both agreement and opposition. Keeping those relations visible helps identify what to investigate next.

[Input and output guide](docs/data-format.md)

## From a published direction to a cohort observation

The literature may report a genus as **more abundant** or **less abundant** in a disease. A participant's observed abundance is classified relative to a healthy-reference interval estimated from **other studies in the same profiling resource**.

| Literature expectation | Observed state | Contribution |
| --- | --- | --- |
| More abundant in disease | High | **+1** — agreement |
| Less abundant in disease | Low | **+1** — agreement |
| More abundant in disease | Low | **−1** — opposition |
| Less abundant in disease | High | **−1** — opposition |
| Either direction | Within range or undetected | **0** |
| No signed direction | Any state | **0** |

**Low abundance is not automatically a negative contribution.** If both the literature and the observation point lower, they agree.

<p align="center">
  <img src="docs/assets/projection.svg" width="520" alt="Six projection cases: higher with High and lower with Low give plus one; opposite directions give minus one; neutral observations or no signed direction give zero.">
</p>

## What you can analyze

- **Publication consensus:** combine extraction records into one vote per publication, then combine publication votes into a signed disease–genus expectation.
- **Participant scores:** estimate external healthy-reference bounds, assign abundance states, and sum the literature-aligned contributions for each disease.
- **Disease patterns:** standardize scores within a study–disease comparison and estimate their association with case status using study-stratified conditional logistic regression.
- **Individual relations:** compare case and control mean states within each study, average these contrasts equally across studies, and summarize agreement, opposition and exactly zero shifts.

The implementation provides inspectable tables at each step, from publication votes and reference bounds to scores and relation-level summaries.

## Try the synthetic example

From this repository's root, install the package in your Python 3.10 or newer environment:

```bash
python -m pip install -e .
```

Run the included example:

```bash
python -m minerva_projection run --participants examples/synthetic/participants.csv --abundance examples/synthetic/abundance.csv --votes examples/synthetic/votes.csv --output demo_output
```

**Every example record is synthetic.** The condition labels, genus labels, publication identifiers and abundance values are invented for demonstration. They are not participant records or extracted literature evidence.

Start with `demo_output/run_summary.json`, then inspect `scores.csv`, `disease_associations.csv` and `relation_alignment.csv`. The [input guide](docs/data-format.md) explains all output files and how to run your own data.

To generate the deterministic example in a different directory:

```bash
python -m minerva_projection synthetic --output my_synthetic_inputs
```

## Analyze your own cohort

Provide three CSV files:

| File | Contents |
| --- | --- |
| Participants | Sample, source and study identifiers; case/control group; disease label |
| Abundance | One row per sample and one column per genus, in **percentage points** |
| Literature votes | Disease–genus extraction records with publication identifiers and directions −1, 0 or +1 |

Include healthy-control studies from the same profiling resource to build an external reference for each target study. Harmonize genus and disease labels across the files before analysis. The default evaluation thresholds are 10 cases, 10 controls from the same study and 100 external reference controls. Relation-landscape summaries require disease support from at least two profiling resources by default; disease-level score associations can use one.

The [input and output guide](docs/data-format.md) includes the schemas, optional publication exclusions, eligibility rules, execution command and a description of every exported table.

## Interpreting agreement

Agreement supplies a reason to investigate a relation further; opposition identifies a question to examine in the source evidence or across cohorts. Neither category alone establishes causality or validates a treatment. Relation-level signs describe the observed direction under this reference rule; they are not significance tests.

This repository implements the core projection and association algorithms for user-supplied data. Its synthetic example demonstrates that workflow; it is not an exact rerun of the manuscript's full analyses.

## Project

Developed in connection with **A literature-to-cohort framework for microbiome discovery**.

Organization: [MGH-LMIC](https://github.com/MGH-LMIC/).

Repository: [MINERVA-Projection](https://github.com/MGH-LMIC/MINERVA-Projection).

Interactively explore the results of this study in [MINERVA Landscape](http://100.30.173.66:3000/).

## License

The code and documentation in this repository are available under the [MIT License](LICENSE).
