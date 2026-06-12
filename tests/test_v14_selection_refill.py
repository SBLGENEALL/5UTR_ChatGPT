import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "01_pipeline/scripts/10_select_2000_cluster_diverse_library.py"


def encoded_sequence(number):
    alphabet = "ACG"
    suffix = ""
    for _ in range(10):
        suffix = alphabet[number % 3] + suffix
        number //= 3
    return "AC" * 25 + suffix


class V14LogicalScoringSelectionTests(unittest.TestCase):
    def test_logical_scoring_reaches_2000_with_floors_and_caps(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            input_csv = workdir / "candidates.csv"
            validation_csv = workdir / "classifier_validation.csv"
            pd.DataFrame(
                [{"selection_metric": 1.0, "split_mode": "gene_seq_cluster_split"}]
            ).to_csv(validation_csv, index=False)
            fields = [
                "utr_id",
                "gene_name",
                "utr5_sequence_tss_corrected",
                "seq_cluster_id",
                "robust_public_te_rank",
                "heavy_ensemble_score",
                "protein_abundance_rank",
                "protein_residual_rank",
                "multi_omics_utr_rank",
                "is_expressed_public",
            ]
            with input_csv.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                for i in range(2400):
                    writer.writerow(
                        {
                            "utr_id": f"utr_{i}",
                            "gene_name": f"gene_{i // 3}",
                            "utr5_sequence_tss_corrected": encoded_sequence(i),
                            "seq_cluster_id": f"cluster_{i // 2}",
                            "robust_public_te_rank": 0.5 + i / 2400 if i < 900 else "",
                            "heavy_ensemble_score": 0.5 + i / 2800 if i < 1100 else "",
                            "protein_abundance_rank": 0.7 if 1600 <= i < 2200 else "",
                            "protein_residual_rank": 0.7 if 1700 <= i < 2200 else "",
                            "multi_omics_utr_rank": 0.7 if i < 1600 else "",
                            "is_expressed_public": True,
                        }
                    )

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--input",
                    str(input_csv),
                    "--n",
                    "2000",
                    "--max-per-cluster",
                    "1",
                    "--allow-cluster-fill",
                    "2",
                    "--max-per-gene",
                    "3",
                    "--classifier-validation",
                    str(validation_csv),
                ],
                cwd=workdir,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, msg=result.stdout + "\n" + result.stderr)

            qc_path = workdir / "07_library_design/tables/v1.4_selection_policy_qc.csv"
            qc = pd.read_csv(qc_path).iloc[0]
            self.assertEqual(int(qc["selected_n"]), 2000)
            self.assertEqual(int(qc["shortage_n"]), 0)
            self.assertEqual(int(qc["J_fill_selected_n"]), 0)
            self.assertGreater(int(qc["multiomics_without_robust_candidate_n"]), 0)
            self.assertGreater(int(qc["multiomics_without_robust_selected_count"]), 0)
            self.assertGreaterEqual(int(qc["total_protein_supported_selected_count"]), 300)
            self.assertGreaterEqual(int(qc["total_classifier_supported_selected_count"]), 800)
            self.assertEqual(int(qc["diversity_exploratory_selected_n"]), 100)
            self.assertEqual(int(qc["group_count_H_negative_controls"]), 50)
            self.assertAlmostEqual(float(qc["classifier_weight"]), 0.30, places=8)
            weight_columns = [c for c in qc.index if c.startswith("weight_")]
            self.assertAlmostEqual(sum(float(qc[c]) for c in weight_columns), 1.0, places=8)
            self.assertTrue(any(c.startswith("length_bin_count_") for c in qc.index))
            self.assertTrue(any(c.startswith("gc_bin_count_") for c in qc.index))
            self.assertLessEqual(int(qc["max_per_gene"]), 3)
            self.assertLessEqual(int(qc["max_per_seq_cluster"]), 2)

            selected_path = workdir / (
                "07_library_design/tables/"
                "selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv"
            )
            selected = pd.read_csv(selected_path)
            selected_multiomics_without_robust = selected[
                selected["multi_omics_utr_rank"].notna() &
                selected["robust_public_te_rank"].isna()
            ]
            self.assertGreater(len(selected_multiomics_without_robust), 0)
            self.assertTrue(selected["composite_evidence_score"].between(0, 1).all())
            self.assertFalse(selected["library_group"].astype(str).str.contains("J_fill").any())

            bin_counts = pd.read_csv(
                workdir / "07_library_design/tables/v1.4_selection_length_gc_bins.csv"
            )
            self.assertTrue(
                {"length_bin", "gc_bin", "library_group", "selection_source",
                 "evidence_source_signature"}.issubset(set(bin_counts["dimension"]))
            )
            comparison = pd.read_csv(
                workdir / "07_library_design/tables/v1.4_pr3_vs_pr2_comparison_ready.csv"
            )
            self.assertTrue(
                {"selected_pr2", "selected_pr3", "status_vs_pr2"}.issubset(comparison.columns)
            )


if __name__ == "__main__":
    unittest.main()
