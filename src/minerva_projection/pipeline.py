"""CSV workflow around the public core, using caller-supplied inputs only."""
from __future__ import annotations
import csv
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .core import (
    require_columns, publication_consensus, healthy_reference, abundance_states,
    literature_scores, standardize_within_study, conditional_logistic_associations,
    equal_study_relation_shifts, directional_alignment, category_summary,
    study_directional_alignment, display_cap,
)


def read_csv(path: Path) -> pd.DataFrame:
    """Reject duplicate headers rather than allowing automatic name mangling."""
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        header = next(csv.reader(handle), [])
    if len(header) != len(set(header)):
        raise ValueError("Duplicate CSV column names")
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def validate_inputs(participants: pd.DataFrame, abundance: pd.DataFrame,
                    votes: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray, list[str], pd.DataFrame]:
    require_columns(participants, ["sample_id", "source", "study_id", "group", "disease"])
    require_columns(abundance, ["sample_id"])
    require_columns(votes, ["record_id", "disease", "genus", "publication_id", "direction"])
    m = participants.copy().reset_index(drop=True)
    for col in ["sample_id", "source", "study_id", "group"]:
        if m[col].isna().any() or m[col].astype(str).str.strip().eq("").any():
            raise ValueError(f"Missing participant {col}")
        m[col] = m[col].astype(str)
    m["disease"] = m.disease.fillna("").astype(str)
    if not m.group.isin(["case", "control"]).all():
        raise ValueError("group must be case or control")
    if m.loc[m.group.eq("case"), "disease"].str.strip().eq("").any():
        raise ValueError("Each case requires a disease label")
    if m.sample_id.duplicated().any() or abundance.sample_id.duplicated().any():
        raise ValueError("sample_id must be unique; resolve biological/technical repeats before use")
    if set(m.sample_id) != set(abundance.sample_id.astype(str)):
        raise ValueError("Participant and abundance sample_id sets must match exactly")
    genera = [str(c) for c in abundance.columns if c != "sample_id"]
    if not genera:
        raise ValueError("At least one genus abundance column is required")
    numeric = abundance.set_index("sample_id").reindex(m.sample_id)[genera].apply(pd.to_numeric, errors="raise")
    values64 = numeric.to_numpy(float)
    if not np.isfinite(values64).all() or (values64 < 0).any() or (values64 > 100).any():
        raise ValueError("Abundance must be finite percentage points between 0 and 100")
    if (values64.sum(axis=1) > 100 + 1e-5).any():
        raise ValueError("Abundance row sums exceed 100 percentage points")
    # Match the manuscript pipeline's percentage-point float32 state inputs.
    values = values64.astype(np.float32)
    v = votes.copy()
    v["direction"] = pd.to_numeric(v.direction, errors="raise")
    # A JSON pair avoids collisions between repeated study names across sources.
    m["study_key"] = [json.dumps([s, t], ensure_ascii=False, separators=(",", ":"))
                      for s, t in zip(m.source, m.study_id)]
    return m, values, genera, v


def run_analysis(participants: pd.DataFrame, abundance: pd.DataFrame,
                 votes: pd.DataFrame, *, excluded_publications: set[str] | None = None,
                 min_group_size: int = 10, min_reference_size: int = 100,
                 display_genera: int = 28, min_landscape_sources: int = 2) -> tuple[dict[str, pd.DataFrame], dict]:
    """Run primary estimands. Ineligible tasks are reported, never silently kept.

    Input deduplication, biological-unit validation, cohort-publication linkage
    and taxonomic/disease harmonization remain caller responsibilities.
    """
    if min(min_group_size, min_reference_size, display_genera, min_landscape_sources) < 1:
        raise ValueError("Thresholds and display_genera must be positive integers")
    m, values, genera, v = validate_inputs(participants, abundance, votes)
    consensus, papers = publication_consensus(v, excluded_publications)
    matrix = consensus.pivot(index="genus", columns="disease", values="direction").reindex(genera).fillna(0)
    tasks = m[m.group.eq("case")][["source", "study_id", "study_key", "disease"]].drop_duplicates()
    tasks = tasks.sort_values(["disease", "source", "study_id"])
    eligibility, refs, scores, states, contributions, shifts = [], [], [], [], [], []
    cached = {}
    for task in tasks.itertuples(index=False):
        same = m.study_key.eq(task.study_key).to_numpy()
        case = np.flatnonzero(same & m.group.eq("case").to_numpy() & m.disease.eq(task.disease).to_numpy())
        control = np.flatnonzero(same & m.group.eq("control").to_numpy())
        refmask = m.source.eq(task.source) & m.group.eq("control") & ~m.study_id.eq(task.study_id)
        k = matrix[task.disease].to_numpy(np.int8) if task.disease in matrix else np.zeros(len(genera), np.int8)
        failures = []
        if len(case) < min_group_size: failures.append("too_few_cases")
        if len(control) < min_group_size: failures.append("too_few_controls")
        if int(refmask.sum()) < min_reference_size: failures.append("too_few_external_controls")
        if not k.any(): failures.append("no_signed_relation_on_supplied_genus_axis")
        eligibility.append(dict(disease=task.disease, source=task.source, study_id=task.study_id,
                                study_key=task.study_key, n_case=len(case), n_control=len(control),
                                n_external_controls=int(refmask.sum()), n_signed_genera=int((k != 0).sum()),
                                eligible=not failures, reasons=";".join(failures)))
        if failures:
            continue
        if task.study_key not in cached:
            reference = healthy_reference(m, values, task.source, task.study_id, min_reference_size)
            cached[task.study_key] = reference
            refs.append(pd.DataFrame(dict(source=task.source, study_id=task.study_id, study_key=task.study_key,
                genus=genera, lower=reference.lower, upper=reference.upper,
                n_external_controls=len(reference.indices), n_detected_controls=reference.detected_counts,
                n_reference_studies=reference.n_reference_studies,
                fallback_no_detections=reference.detected_counts == 0)))
        ref = cached[task.study_key]
        indices = np.r_[case, control]
        labels = np.r_[np.ones(len(case), int), np.zeros(len(control), int)]
        signal = abundance_states(values[indices], ref.lower, ref.upper)
        contribution, raw = literature_scores(signal, k)
        block = m.iloc[indices][["sample_id", "source", "study_id", "study_key"]].copy()
        block["disease"], block["label"], block["target_score"] = task.disease, labels, raw
        scores.append(block)
        keys = pd.DataFrame(dict(sample_id=np.repeat(m.iloc[indices].sample_id.to_numpy(), len(genera)),
            source=task.source, study_id=task.study_id, study_key=task.study_key, disease=task.disease,
            label=np.repeat(labels, len(genera)), genus=np.tile(genera, len(indices))))
        state_table = keys.copy()
        state_table["abundance_percent"] = values[indices].ravel()
        state_table["state"] = signal.ravel()
        states.append(state_table)
        contribution_table = keys.copy()
        contribution_table["direction"] = np.tile(k, len(indices))
        contribution_table["contribution"] = contribution.ravel()
        contributions.append(contribution_table)
        delta = signal[labels == 1].mean(axis=0) - signal[labels == 0].mean(axis=0)
        shifts.append(pd.DataFrame(dict(disease=task.disease, source=task.source, study_key=task.study_key,
            study_id=task.study_id, genus=genera, observed_shift=delta, direction=k)))
    outputs = {"consensus": consensus, "publication_votes": papers,
               "task_eligibility": pd.DataFrame(eligibility),
               "references": pd.concat(refs, ignore_index=True) if refs else pd.DataFrame()}
    summary = dict(data_origin="user_supplied", abundance_unit="percentage_points",
        input_samples=len(m), input_genera=len(genera), input_extraction_records=len(v),
        excluded_publications_requested=len(excluded_publications or set()),
        retained_publications=int(papers.publication_id.nunique()),
        retained_relations=int(consensus.direction.ne(0).sum()),
        projection_relations_on_genus_axis=int(np.count_nonzero(matrix.to_numpy())),
        candidate_tasks=len(tasks), eligible_tasks=len(scores),
        min_group_size=min_group_size, min_reference_size=min_reference_size,
        min_landscape_sources=min_landscape_sources,
        primary_thresholds_used=min_group_size == 10 and min_reference_size == 100,
        scope="primary projection and association; no automatic evidence extraction or full manuscript reproduction")
    if not scores:
        summary["status"] = "no_eligible_tasks"
        return outputs, summary
    pred = standardize_within_study(pd.concat(scores, ignore_index=True))
    fp = pd.concat(shifts, ignore_index=True)
    study_alignment, disease_alignment = study_directional_alignment(fp)
    outputs.update(scores=pred, states=pd.concat(states, ignore_index=True),
        contributions=pd.concat(contributions, ignore_index=True),
        disease_associations=conditional_logistic_associations(pred), study_relation_shifts=fp,
        study_directional_alignment=study_alignment,
        disease_directional_alignment=disease_alignment)
    support = pred.groupby("disease").source.nunique()
    landscape_diseases = support[support >= min_landscape_sources].index
    landscape_fp = fp[fp.disease.isin(landscape_diseases)]
    summary.update(status="complete", evaluated_unique_samples=int(pred.sample_id.nunique()),
        evaluated_comparison_rows=len(pred), evaluated_diseases=int(pred.disease.nunique()),
        landscape_diseases=len(landscape_diseases))
    if not landscape_fp.empty:
        relations = equal_study_relation_shifts(landscape_fp)
        cap, selected = display_cap(relations, display_genera)
        cells = directional_alignment(relations, cap)
        cells["display_selected"] = cells.genus.isin(selected)
        outputs.update(relation_alignment=cells,
                       disease_alignment_summary=category_summary(cells, "disease"),
                       genus_alignment_summary=category_summary(cells, "genus"))
        summary.update(display_cap=cap, display_genera=selected,
                       display_cap_quantile=.95, display_cap_floor=1e-8)
    return outputs, summary


def write_outputs(output: Path, tables: dict[str, pd.DataFrame], summary: dict) -> None:
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output directory is not empty; choose a new directory")
    output.mkdir(parents=True, exist_ok=True)
    for name, table in tables.items():
        table.to_csv(output / f"{name}.csv", index=False)
    (output / "run_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False,
                                                       allow_nan=False) + "\n", encoding="utf-8")
