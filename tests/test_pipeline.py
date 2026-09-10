from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from minerva_projection.synthetic import generate_synthetic
from minerva_projection.pipeline import read_csv, run_analysis, validate_inputs, write_outputs
from minerva_projection.cli import main


class PipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.path = Path(cls.tmp.name) / "synthetic"
        generate_synthetic(cls.path)
        cls.m = read_csv(cls.path / "participants.csv")
        cls.a = read_csv(cls.path / "abundance.csv")
        cls.v = read_csv(cls.path / "votes.csv")

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_synthetic_primary_workflow(self):
        tables, summary = run_analysis(self.m, self.a, self.v)
        self.assertEqual(summary["status"], "complete")
        self.assertEqual(summary["eligible_tasks"], 8)
        self.assertEqual(summary["landscape_diseases"], 2)
        self.assertTrue(tables["task_eligibility"].n_external_controls.ge(100).all())
        self.assertTrue(tables["task_eligibility"].n_case.ge(10).all())
        self.assertEqual(len(tables["relation_alignment"]), 16)
        self.assertTrue(tables["relation_alignment"].query("genus == 'GenusF'").category.eq("unlinked").all())
        self.assertTrue(tables["states"].query("genus == 'GenusG'").state.eq(0).all())
        self.assertTrue(tables["references"].query("genus == 'GenusG'").fallback_no_detections.all())
        self.assertTrue(tables["disease_associations"].status.eq("estimable").all())
        grouping = ["sample_id", "disease", "study_key"]
        raw = tables["contributions"].groupby(grouping).contribution.sum().sort_index()
        scores = tables["scores"].set_index(grouping).target_score.sort_index()
        np.testing.assert_array_equal(raw, scores)
        output = Path(self.tmp.name) / "result"
        write_outputs(output, tables, summary)
        self.assertTrue((output / "run_summary.json").exists())
        with self.assertRaisesRegex(ValueError, "not empty"):
            write_outputs(output, tables, summary)

    def test_reference_ineligibility_is_reported(self):
        tables, summary = run_analysis(self.m, self.a, self.v, min_reference_size=1000)
        self.assertEqual(summary["status"], "no_eligible_tasks")
        self.assertTrue(tables["task_eligibility"].reasons.str.contains("too_few_external_controls").all())
        self.assertNotIn("scores", tables)

    def test_group_ineligibility_is_reported(self):
        tables, summary = run_analysis(self.m, self.a, self.v, min_group_size=20)
        self.assertEqual(summary["eligible_tasks"], 0)
        self.assertTrue(tables["task_eligibility"].reasons.str.contains("too_few_cases").all())

    def test_zero_direction_is_not_a_fitted_task(self):
        v = self.v.copy()
        v.direction = 0
        tables, summary = run_analysis(self.m, self.a, v)
        self.assertEqual(summary["eligible_tasks"], 0)
        self.assertTrue(tables["task_eligibility"].reasons.str.contains("no_signed_relation").all())

    def test_missing_abundance_and_duplicate_samples_fail(self):
        with self.assertRaisesRegex(ValueError, "sets must match"):
            validate_inputs(self.m, self.a.iloc[1:], self.v)
        with self.assertRaisesRegex(ValueError, "must be unique"):
            validate_inputs(pd.concat([self.m, self.m.iloc[:1]]), self.a, self.v)
        bad = self.a.copy()
        bad.loc[0, "GenusA"] = "nan"
        with self.assertRaises(ValueError):
            validate_inputs(self.m, bad, self.v)

    def test_abundance_row_order_is_aligned_by_id(self):
        _, values, _, _ = validate_inputs(self.m, self.a, self.v)
        _, shuffled, _, _ = validate_inputs(self.m, self.a.iloc[::-1], self.v)
        np.testing.assert_array_equal(values, shuffled)

    def test_synthetic_generation_is_repeatable(self):
        other = Path(self.tmp.name) / "again"
        generate_synthetic(other)
        for name in ["participants.csv", "abundance.csv", "votes.csv"]:
            self.assertEqual((other / name).read_bytes(), (self.path / name).read_bytes())

    def test_duplicate_csv_header_is_rejected(self):
        path = Path(self.tmp.name) / "duplicate.csv"
        path.write_text("sample_id,G,G\nx,1,2\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Duplicate CSV"):
            read_csv(path)

    def test_cli_synthetic_generation(self):
        output = Path(self.tmp.name) / "cli"
        self.assertEqual(main(["synthetic", "--output", str(output)]), 0)
        self.assertTrue((output / "SYNTHETIC_DATA_NOTICE.json").exists())


if __name__ == "__main__":
    unittest.main()
