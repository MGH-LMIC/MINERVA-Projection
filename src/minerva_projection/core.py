"""Literature projection primitives with explicit input and denominator checks.

Abundances are percentage points. Literature directions are supplied by the
caller; this module does not extract, retrieve or validate scientific evidence.
"""
from __future__ import annotations

from dataclasses import dataclass
import warnings
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm
from statsmodels.discrete.conditional_models import ConditionalLogit
from statsmodels.stats.multitest import multipletests
from statsmodels.tools.sm_exceptions import ConvergenceWarning


def require_columns(frame: pd.DataFrame, columns: list[str]) -> None:
    missing = set(columns) - set(frame.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    if frame.columns.duplicated().any():
        raise ValueError("Duplicate column names are not allowed")


def require_signs(values: np.ndarray, name: str) -> np.ndarray:
    a = np.asarray(values)
    if not np.isin(a, [-1, 0, 1]).all():
        raise ValueError(f"{name} must contain only -1, 0, +1")
    return a.astype(np.int8)


def publication_consensus(
    records: pd.DataFrame, excluded_publications: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One signed vote per publication per relation, then equal-publication sum.

    Required record columns: record_id, disease, genus, publication_id, direction.
    record_id identifies an extraction record and must be globally unique.
    Zero-sum publications contribute zero; zero-sum relations have direction zero.
    Exclusion is performed before either aggregation and applies globally.
    Returns (relation_consensus, publication_votes).
    """
    keys = ["disease", "genus", "publication_id"]
    require_columns(records, ["record_id", *keys, "direction"])
    r = records.copy()
    for c in ["record_id", *keys]:
        if r[c].isna().any() or r[c].astype(str).str.strip().eq("").any():
            raise ValueError(f"{c} must be present and non-empty")
        r[c] = r[c].astype(str)
    if r.record_id.duplicated().any():
        raise ValueError("Duplicate record_id; resolve extraction identity before analysis")
    r["direction"] = require_signs(r.direction.to_numpy(), "direction")
    if excluded_publications:
        r = r[~r.publication_id.isin({str(x) for x in excluded_publications})]
    papers = r.groupby(keys, as_index=False, sort=True).agg(
        record_direction_sum=("direction", "sum"), n_records=("direction", "size"))
    papers["publication_vote"] = np.sign(papers.record_direction_sum).astype(np.int8)
    relations = papers.groupby(["disease", "genus"], as_index=False, sort=True).agg(
        vote_sum=("publication_vote", "sum"), n_publications=("publication_id", "size"),
        n_positive_votes=("publication_vote", lambda x: int((x > 0).sum())),
        n_negative_votes=("publication_vote", lambda x: int((x < 0).sum())),
        n_zero_votes=("publication_vote", lambda x: int((x == 0).sum())))
    relations["direction"] = np.sign(relations.vote_sum).astype(np.int8)
    return relations, papers


def present_conditional_bounds(reference: np.ndarray, q_low: float = 5,
                               q_high: float = 95) -> tuple[np.ndarray, np.ndarray]:
    """Linear quantiles among strictly positive values; float32 stored bounds.

    An undetected genus uses lower=0 and upper=1e-6 percentage points. This
    fallback does not establish a biological normal interval.
    """
    a = np.asarray(reference)
    if a.ndim != 2 or not np.isfinite(a).all() or (a < 0).any():
        raise ValueError("reference must be a finite, nonnegative 2D matrix")
    if not 0 <= q_low < q_high <= 100:
        raise ValueError("Require 0 <= q_low < q_high <= 100")
    lower = np.zeros(a.shape[1], dtype=np.float32)
    upper = np.full(a.shape[1], 1e-6, dtype=np.float32)
    for j in range(a.shape[1]):
        positive = a[:, j][a[:, j] > 0]
        if positive.size:
            lower[j], upper[j] = np.quantile(positive, [q_low / 100, q_high / 100])
    return lower, upper


@dataclass(frozen=True)
class HealthyReference:
    lower: np.ndarray
    upper: np.ndarray
    indices: np.ndarray
    detected_counts: np.ndarray
    n_reference_studies: int


def healthy_reference(metadata: pd.DataFrame, abundance: np.ndarray,
                      target_source: str, target_study: str,
                      min_reference_size: int = 100) -> HealthyReference:
    """Select controls in the same source, excluding the entire target study.

    The minimum applies to the total pool, not each genus's detected count.
    metadata and abundance must already be aligned row-for-row.
    """
    require_columns(metadata, ["source", "study_id", "group"])
    a = np.asarray(abundance)
    if a.ndim != 2 or len(metadata) != a.shape[0]:
        raise ValueError("metadata and abundance rows must align")
    if not metadata.group.isin(["case", "control"]).all():
        raise ValueError("group must be case or control")
    if min_reference_size < 1:
        raise ValueError("min_reference_size must be positive")
    mask = (metadata.source.eq(target_source) & metadata.group.eq("control")
            & ~metadata.study_id.eq(target_study)).to_numpy()
    indices = np.flatnonzero(mask)
    if len(indices) < min_reference_size:
        raise ValueError("Insufficient external same-source controls")
    lower, upper = present_conditional_bounds(a[indices])
    return HealthyReference(lower, upper, indices, (a[indices] > 0).sum(axis=0),
                            int(metadata.iloc[indices].study_id.nunique()))


def abundance_states(abundance: np.ndarray, lower: np.ndarray,
                     upper: np.ndarray) -> np.ndarray:
    """Low=-1; neutral=0; high=+1. Absence is neutral; low wins tied bounds."""
    a = np.asarray(abundance)
    lo, hi = np.asarray(lower), np.asarray(upper)
    if a.ndim != 2 or lo.shape != (a.shape[1],) or hi.shape != lo.shape:
        raise ValueError("Bounds must be 1D and match the genus dimension")
    if not all(np.isfinite(x).all() for x in (a, lo, hi)) or (a < 0).any():
        raise ValueError("Abundances and bounds must be finite; abundances nonnegative")
    if (lo < 0).any() or (lo > hi).any():
        raise ValueError("Require 0 <= lower <= upper")
    low = (a > 0) & (a <= lo)
    high = (a > 0) & ~low & (a >= hi)
    return high.astype(np.int8) - low.astype(np.int8)


def literature_scores(states: np.ndarray, directions: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return signed contributions and their raw sum; no abundance weighting."""
    s, k = require_signs(states, "states"), require_signs(directions, "directions")
    if s.ndim != 2 or k.shape != (s.shape[1],):
        raise ValueError("directions must match the genus dimension of states")
    contributions = s * k
    # Wide accumulation avoids integer overflow for arbitrary user-supplied axes.
    return contributions, contributions.sum(axis=1, dtype=np.int64)


def standardize_within_study(scores: pd.DataFrame) -> pd.DataFrame:
    """Standardize pooled cases and controls inside each disease/study task.

    Required columns: disease, study_key, target_score. Constant or single-row
    tasks remain z=0 and are marked uninformative; ddof=1, threshold 1e-10.
    """
    require_columns(scores, ["disease", "study_key", "target_score"])
    out = scores.reset_index(drop=True).copy()
    if not np.isfinite(out.target_score.to_numpy(float)).all():
        raise ValueError("Scores must be finite")
    out["within_study_z"] = 0.0
    out["informative_study"] = False
    for index in out.groupby(["disease", "study_key"], sort=True).groups.values():
        values = out.loc[index, "target_score"].to_numpy(float)
        sd = float(values.std(ddof=1)) if len(values) > 1 else np.nan
        if np.isfinite(sd) and sd > 1e-10:
            out.loc[index, "within_study_z"] = (values - values.mean()) / sd
            out.loc[index, "informative_study"] = True
    return out


def _conditional_fit(group: pd.DataFrame) -> dict:
    fields = ["pooled_log_or", "pooled_se", "pooled_or", "ci_low", "ci_high",
              "wald_p_two_sided", "likelihood_ratio_p"]
    base = {x: np.nan for x in fields}
    data = group[group.informative_study].copy()
    if data.empty:
        return {**base, "status": "not_estimable_no_informative_study"}
    if data.label.nunique() < 2:
        return {**base, "status": "not_estimable_no_case_control_contrast"}
    # Primary tasks require cases and controls inside every stratum.
    if data.groupby("study_key").label.nunique().lt(2).any():
        raise ValueError("Every disease/study task must include both groups")
    try:
        model = ConditionalLogit(data.label.to_numpy(int),
            data[["within_study_z"]].to_numpy(float),
            groups=pd.factorize(data.study_key, sort=True)[0])
        with warnings.catch_warnings(record=True) as fit_warnings:
            warnings.simplefilter("always")
            result = model.fit(method="bfgs", maxiter=1000, disp=False)
        # ConditionalResults does not retain mle_retvals in every statsmodels
        # version, so also inspect the optimizer's explicit warning category.
        reported = getattr(result, "mle_retvals", {}) or {}
        failed = reported.get("converged") is False or any(
            issubclass(w.category, ConvergenceWarning) for w in fit_warnings)
        if failed:
            return {**base, "status": "not_estimable_nonconverged"}
        beta, se = float(result.params[0]), float(result.bse[0])
        if not np.isfinite([beta, se]).all() or se <= 0:
            return {**base, "status": "not_estimable_nonfinite_estimate"}
        lr = max(0., 2 * (float(result.llf) - float(model.loglike(np.zeros(1)))))
        return dict(pooled_log_or=beta, pooled_se=se, pooled_or=float(np.exp(beta)),
                    ci_low=float(np.exp(beta - 1.96 * se)), ci_high=float(np.exp(beta + 1.96 * se)),
                    wald_p_two_sided=float(2 * norm.sf(abs(beta / se))),
                    likelihood_ratio_p=float(chi2.sf(lr, 1)), status="estimable")
    except (ValueError, np.linalg.LinAlgError, RuntimeError) as error:
        return {**base, "status": f"not_estimable_{type(error).__name__}"}


def conditional_logistic_associations(scores: pd.DataFrame) -> pd.DataFrame:
    """Primary unadjusted conditional OR, two-sided Wald P, 95% CI and BH.

    All supplied diseases form one BH family; non-estimable P values enter
    that adjustment as 1, while their unadjusted P stays missing.
    Explicit optimizer nonconvergence is returned as non-estimable, even if
    the solver supplied finite coefficients.
    """
    require_columns(scores, ["disease", "study_key", "target_score", "label"])
    if not scores.label.isin([0, 1]).all():
        raise ValueError("label must be 0 or 1")
    if scores.empty:
        raise ValueError("No tasks supplied")
    data = standardize_within_study(scores)
    rows = []
    for disease, group in data.groupby("disease", sort=True):
        rows.append(dict(disease=disease, n_studies=int(group.study_key.nunique()),
                         n_case=int(group.label.eq(1).sum()), n_control=int(group.label.eq(0).sum()),
                         **_conditional_fit(group)))
    out = pd.DataFrame(rows)
    out["fdr_bh"] = multipletests(out.wald_p_two_sided.fillna(1), method="fdr_bh")[1]
    return out


def equal_study_relation_shifts(study_shifts: pd.DataFrame) -> pd.DataFrame:
    """Mean case-minus-control abundance-state shift; each study has weight one."""
    cols = ["disease", "study_key", "genus", "observed_shift", "direction"]
    require_columns(study_shifts, cols)
    if study_shifts.duplicated(["disease", "study_key", "genus"]).any():
        raise ValueError("One shift per disease/study/genus is required")
    require_signs(study_shifts.direction.to_numpy(), "direction")
    if not np.isfinite(study_shifts.observed_shift.to_numpy(float)).all():
        raise ValueError("Shifts must be finite")
    grouped = study_shifts.groupby(["disease", "genus"], sort=True)
    if grouped.direction.nunique().gt(1).any():
        raise ValueError("Literature direction must be fixed across studies")
    return grouped.agg(observed_shift=("observed_shift", "mean"),
                       direction=("direction", "first"),
                       n_studies=("study_key", "nunique")).reset_index()


def display_cap(relations: pd.DataFrame, n_genera: int = 28) -> tuple[float, list[str]]:
    """Select genera by sum of cohort-magnitude and literature-degree ranks.

    Compute the 95th percentile of absolute selected-cell shifts, with a 1e-8
    floor. Only display scaling uses this cap; sign categories do not.
    """
    require_columns(relations, ["disease", "genus", "observed_shift", "direction"])
    if relations.empty or n_genera < 1:
        raise ValueError("At least one relation and one display genus are required")
    obs = relations.pivot(index="disease", columns="genus", values="observed_shift")
    kg = relations.pivot(index="disease", columns="genus", values="direction")
    if obs.isna().any().any() or kg.isna().any().any():
        raise ValueError("A complete disease-by-genus coordinate is required")
    rank = obs.abs().sum().rank(pct=True) + kg.ne(0).sum().rank(pct=True)
    selected = rank.nlargest(n_genera).index.tolist()
    cap = max(float(np.quantile(np.abs(obs[selected].to_numpy()), .95)), 1e-8)
    return cap, selected


def directional_alignment(relations: pd.DataFrame, cap: float) -> pd.DataFrame:
    """Separate unlinked from exact zero; graded alignment is clip(A/c)*K."""
    require_columns(relations, ["disease", "genus", "observed_shift", "direction"])
    if not np.isfinite(cap) or cap <= 0:
        raise ValueError("cap must be finite and positive")
    out = relations.copy()
    k = require_signs(out.direction.to_numpy(), "direction")
    if not np.isfinite(out.observed_shift.to_numpy(float)).all():
        raise ValueError("Shifts must be finite")
    out["linked"] = k != 0
    out["scaled_shift"] = np.clip(out.observed_shift / cap, -1, 1)
    out["graded_alignment"] = np.where(out.linked, out.scaled_shift * k, np.nan)
    signed = np.sign(out.observed_shift.to_numpy()) * k
    out["category"] = np.select([~out.linked, signed > 0, signed < 0],
                                 ["unlinked", "aligned", "opposed"], default="exact_zero")
    return out


def category_summary(cells: pd.DataFrame, entity: str) -> pd.DataFrame:
    """Count sign categories among linked cells; report coverage separately."""
    if entity not in {"disease", "genus"}:
        raise ValueError("entity must be disease or genus")
    require_columns(cells, [entity, "linked", "category", "graded_alignment"])
    rows = []
    for name, group in cells.groupby(entity, sort=True):
        linked = group[group.linked]
        row = {entity: name, "n_coordinate_cells": len(group), "n_linked": len(linked),
               "coverage_fraction": len(linked) / len(group),
               "mean_graded_alignment": float(linked.graded_alignment.mean())}
        for label in ["aligned", "opposed", "exact_zero"]:
            count = int(linked.category.eq(label).sum())
            row[f"n_{label}"] = count
            row[f"{label}_fraction_among_linked"] = count / len(linked) if len(linked) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def study_directional_alignment(study_shifts: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Mean sign(A_study)*K per study, then equal-study disease mean.

    This differs from classifying the sign of the study-averaged shift. Zero
    shifts stay in the primary denominator; an explicitly labelled nonzero
    sensitivity summary is returned alongside it. No bootstrap CI is fitted.
    """
    # Reuse validation but retain individual studies for this estimand.
    equal_study_relation_shifts(study_shifts)
    x = study_shifts[study_shifts.direction.ne(0)].copy()
    x["signed"] = np.sign(x.observed_shift) * x.direction
    rows = []
    for (disease, study), group in x.groupby(["disease", "study_key"], sort=True):
        nonzero = group[group.signed.ne(0)]
        rows.append(dict(disease=disease, study_key=study, n_relations=len(group),
                         directional_alignment=float(group.signed.mean()),
                         nonzero_only_alignment=float(nonzero.signed.mean())))
    study = pd.DataFrame(rows, columns=["disease", "study_key", "n_relations", "directional_alignment", "nonzero_only_alignment"])
    disease = study.groupby("disease", as_index=False, sort=True).agg(
        equal_study_directional_alignment=("directional_alignment", "mean"),
        equal_study_nonzero_only_alignment=("nonzero_only_alignment", "mean"),
        n_studies=("study_key", "size"))
    return study, disease
