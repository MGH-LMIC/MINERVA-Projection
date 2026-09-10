"""Boundary and analytic checks using invented values only.

Run from the repository root: python -m unittest discover -s tests -v
"""
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
from types import SimpleNamespace
import warnings
import numpy as np
import pandas as pd
from statsmodels.tools.sm_exceptions import ConvergenceWarning

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from minerva_projection.core import (
    publication_consensus, present_conditional_bounds, healthy_reference,
    abundance_states, literature_scores, standardize_within_study,
    conditional_logistic_associations, equal_study_relation_shifts,
    directional_alignment, category_summary, study_directional_alignment, display_cap,
)


class CoreTests(unittest.TestCase):
    def test_publication_balance_and_global_exclusion(self):
        records = pd.DataFrame([
            ("r1", "D", "G", "synthetic-paper-a", 1),
            ("r2", "D", "G", "synthetic-paper-a", 1),
            ("r3", "D", "G", "synthetic-paper-a", 1),
            ("r4", "D", "G", "synthetic-paper-b", -1),
            ("r5", "D", "G", "synthetic-paper-c", -1),
            ("r6", "D", "G", "synthetic-paper-tie", 1),
            ("r7", "D", "G", "synthetic-paper-tie", -1),
        ], columns=["record_id", "disease", "genus", "publication_id", "direction"])
        relations, papers = publication_consensus(records)
        self.assertEqual(relations.direction.iloc[0], -1)
        self.assertEqual(relations.vote_sum.iloc[0], -1)
        self.assertEqual(papers.publication_vote.tolist(), [1, -1, -1, 0])
        changed, _ = publication_consensus(records, {"synthetic-paper-b", "synthetic-paper-c"})
        self.assertEqual(changed.direction.iloc[0], 1)

    def test_repeated_record_is_not_silently_double_counted(self):
        r = pd.DataFrame([["r", "D", "G", "P", 1]] * 2,
                         columns=["record_id", "disease", "genus", "publication_id", "direction"])
        with self.assertRaisesRegex(ValueError, "Duplicate record_id"):
            publication_consensus(r)

    def test_detected_quantiles_and_empty_genus_fallback(self):
        reference = np.array([[0, 0], [1, 0], [2, 0], [3, 0], [4, 0], [5, 0]], np.float32)
        low, high = present_conditional_bounds(reference)
        np.testing.assert_array_equal(low, np.array([1.2, 0], np.float32))
        np.testing.assert_array_equal(high, np.array([4.8, 1e-6], np.float32))

    def test_external_reference_excludes_target_and_other_source(self):
        m = pd.DataFrame({"source": ["A"] * 101 + ["B"],
                          "study_id": ["target"] + ["external"] * 100 + ["external"],
                          "group": ["control"] * 102})
        a = np.array([99] + list(range(1, 101)) + [99], np.float32)[:, None]
        r = healthy_reference(m, a, "A", "target")
        np.testing.assert_array_equal(r.indices, np.arange(1, 101))
        np.testing.assert_array_equal(r.lower, np.array([5.95], np.float32))
        np.testing.assert_array_equal(r.upper, np.array([95.05], np.float32))
        with self.assertRaisesRegex(ValueError, "Insufficient"):
            healthy_reference(m.iloc[:-2], a[:-2], "A", "target")

    def test_all_projection_cases(self):
        # lower+low, higher+high, lower+high, higher+low, within, absent, unlinked.
        a = np.array([[1, 9, 9, 1, 5, 0, 9]], np.float32)
        state = abundance_states(a, np.full(7, 2), np.full(7, 8))
        c, total = literature_scores(state, np.array([-1, 1, -1, 1, 1, -1, 0]))
        np.testing.assert_array_equal(state, [[-1, 1, 1, -1, 0, 0, 1]])
        np.testing.assert_array_equal(c, [[1, 1, -1, -1, 0, 0, 0]])
        np.testing.assert_array_equal(total, [0])

    def test_equal_bounds_low_priority_and_inclusive_edges(self):
        a = np.array([[0], [1], [2]], np.float32)
        np.testing.assert_array_equal(abundance_states(a, np.array([1]), np.array([1])), [[0], [-1], [1]])
        np.testing.assert_array_equal(abundance_states(np.array([[2], [8]]), np.array([2]), np.array([8])), [[-1], [1]])

    def test_invalid_numeric_inputs_fail(self):
        for bad in [np.array([[np.nan]]), np.array([[-1]])]:
            with self.assertRaises(ValueError):
                abundance_states(bad, np.array([0]), np.array([1]))
        with self.assertRaises(ValueError):
            literature_scores(np.array([[2]]), np.array([1]))

    def test_score_accumulation_cannot_wrap_int16(self):
        _, score = literature_scores(np.ones((1, 40000), np.int8), np.ones(40000, np.int8))
        self.assertEqual(score[0], 40000)

    def test_ddof_one_and_uninformative_strata(self):
        x = pd.DataFrame(dict(disease=["D"] * 6, study_key=["A"] * 3 + ["B"] * 2 + ["C"],
                              target_score=[0, 2, 4, 7, 7, 9]))
        out = standardize_within_study(x)
        np.testing.assert_array_equal(out.within_study_z, [-1, 0, 1, 0, 0, 0])
        self.assertEqual(out.informative_study.tolist(), [True] * 3 + [False] * 3)

    def test_conditional_logistic_against_analytic_matched_pairs(self):
        rows = []
        for i in range(4):
            labels = [0, 1] if i < 3 else [1, 0]
            for score, label in zip([0, 1], labels):
                rows.append(dict(disease="D", study_key=f"pair-{i}", target_score=score, label=label))
        # Three favorable and one unfavorable matched pair: beta=log(3)/sqrt(2).
        out = conditional_logistic_associations(pd.DataFrame(rows)).iloc[0]
        self.assertEqual(out.status, "estimable")
        self.assertAlmostEqual(out.pooled_log_or, np.log(3) / np.sqrt(2), places=4)
        self.assertAlmostEqual(out.pooled_se, np.sqrt(2 / 3), places=4)
        self.assertAlmostEqual(out.fdr_bh, out.wald_p_two_sided, places=12)

    def test_bh_keeps_uninformative_disease_in_family(self):
        rows = []
        for disease, scores in [("variable", [0, 1, 0, 1]), ("constant", [1, 1, 1, 1])]:
            for i, score in enumerate(scores):
                rows.append(dict(disease=disease, study_key="study", target_score=score, label=i // 2))
        out = conditional_logistic_associations(pd.DataFrame(rows)).set_index("disease")
        self.assertEqual(out.loc["constant", "status"], "not_estimable_no_informative_study")
        self.assertTrue(np.isnan(out.loc["constant", "wald_p_two_sided"]))
        self.assertEqual(out.loc["constant", "fdr_bh"], 1)

    def test_finite_but_nonconverged_fit_is_not_reported_as_estimable(self):
        x = pd.DataFrame(dict(disease=["D"] * 4, study_key=["A"] * 4,
                              target_score=[0, 1, 0, 1], label=[0, 0, 1, 1]))
        def failed_fit(*args, **kwargs):
            warnings.warn("Synthetic optimizer failure", ConvergenceWarning)
            return SimpleNamespace(params=np.array([.5]), bse=np.array([.1]))
        with patch("minerva_projection.core.ConditionalLogit.fit", side_effect=failed_fit):
            result = conditional_logistic_associations(x).iloc[0]
        self.assertEqual(result.status, "not_estimable_nonconverged")
        self.assertTrue(np.isnan(result.pooled_log_or))
        self.assertEqual(result.fdr_bh, 1)

    def test_equal_study_weight_and_noncommuting_direction(self):
        fp = pd.DataFrame(dict(disease=["D", "D"], study_key=["small", "large"], genus=["G", "G"],
                               observed_shift=[-1., .2], direction=[1, 1], n_samples=[10, 1000]))
        agg = equal_study_relation_shifts(fp)
        self.assertAlmostEqual(agg.observed_shift.iloc[0], -.4)
        self.assertEqual(directional_alignment(agg, .5).category.iloc[0], "opposed")
        _, d = study_directional_alignment(fp)
        self.assertEqual(d.equal_study_directional_alignment.iloc[0], 0)

    def test_unlinked_and_exact_zero_denominators(self):
        r = pd.DataFrame(dict(disease=["D"] * 4, genus=list("ABCD"),
                              observed_shift=[.2, 0, .8, -.1], direction=[1, 1, 0, 1]))
        cells = directional_alignment(r, .1)
        self.assertEqual(cells.category.tolist(), ["aligned", "exact_zero", "unlinked", "opposed"])
        self.assertTrue(np.isnan(cells.graded_alignment.iloc[2]))
        summary = category_summary(cells, "disease").iloc[0]
        self.assertEqual(summary.n_linked, 3)
        self.assertEqual(summary.coverage_fraction, .75)
        self.assertAlmostEqual(summary.aligned_fraction_among_linked, 1 / 3)
        self.assertAlmostEqual(summary.exact_zero_fraction_among_linked, 1 / 3)
        self.assertEqual(cells.category.tolist(), directional_alignment(r, 100).category.tolist())

    def test_display_cap_quantile_and_zero_floor(self):
        r = pd.DataFrame(dict(disease=["A", "A", "B", "B"], genus=["G1", "G2"] * 2,
                              observed_shift=[0, .2, -.1, .4], direction=[1, 1, -1, 1]))
        cap, selected = display_cap(r, 2)
        self.assertAlmostEqual(cap, .37)
        self.assertEqual(set(selected), {"G1", "G2"})
        r.observed_shift = 0
        self.assertEqual(display_cap(r, 2)[0], 1e-8)

    def test_inconsistent_directions_and_duplicate_shifts_fail(self):
        x = pd.DataFrame(dict(disease=["D"] * 2, genus=["G"] * 2, study_key=["A", "B"],
                              observed_shift=[0, 0], direction=[-1, 1]))
        with self.assertRaisesRegex(ValueError, "fixed across studies"):
            equal_study_relation_shifts(x)
        x.direction = 1
        x.study_key = "A"
        with self.assertRaisesRegex(ValueError, "One shift"):
            equal_study_relation_shifts(x)


if __name__ == "__main__":
    unittest.main()
