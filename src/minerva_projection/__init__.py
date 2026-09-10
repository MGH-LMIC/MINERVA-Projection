"""Reusable algorithms; no built-in literature graph or participant data."""

from .core import (
    abundance_states, healthy_reference, literature_scores,
    publication_consensus, present_conditional_bounds,
    standardize_within_study, conditional_logistic_associations,
    equal_study_relation_shifts, directional_alignment,
    category_summary, study_directional_alignment, display_cap,
)

__version__ = "0.1.0"
__all__ = [
    "abundance_states", "healthy_reference", "literature_scores",
    "publication_consensus", "present_conditional_bounds",
    "standardize_within_study", "conditional_logistic_associations",
    "equal_study_relation_shifts", "directional_alignment",
    "category_summary", "study_directional_alignment", "display_cap",
]
