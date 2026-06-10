import importlib.util
from pathlib import Path
import unittest

import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "01_pipeline/scripts/11_v14_qc_audit.py"
SPEC = importlib.util.spec_from_file_location("v14_qc_audit", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class V14QcAuditTests(unittest.TestCase):
    def test_same_chromosome_distinct_loci_are_multi_mapping(self):
        selected = pd.DataFrame([
            {
                "_query_id": "lib1",
                "library_index": 1,
                "gene_name": "GENE1",
                "library_group": "A_publicTE_high_confidence",
            }
        ])
        hits = pd.DataFrame([
            {
                "query_id": "lib1",
                "target_id": "chr1",
                "target_start": 100,
                "target_end": 160,
                "identity": 0.99,
                "coverage": 1.0,
                "bitscore": 100,
            },
            {
                "query_id": "lib1",
                "target_id": "chr1",
                "target_start": 1000,
                "target_end": 1060,
                "identity": 0.98,
                "coverage": 1.0,
                "bitscore": 95,
            },
        ])
        summary = MODULE.summarize_mapping(selected, hits, 0.95, 0.90)
        self.assertTrue(bool(summary.loc[0, "mapped"]))
        self.assertTrue(bool(summary.loc[0, "multi_mapping"]))
        self.assertEqual(int(summary.loc[0, "qualifying_locus_count"]), 2)

    def test_unmapped_sequence_is_flagged(self):
        selected = pd.DataFrame([
            {
                "_query_id": "lib1",
                "library_index": 1,
                "gene_name": "GENE1",
                "library_group": "A_publicTE_high_confidence",
            }
        ])
        summary = MODULE.summarize_mapping(selected, pd.DataFrame(), 0.95, 0.90)
        self.assertFalse(bool(summary.loc[0, "mapped"]))
        self.assertTrue(bool(summary.loc[0, "suspected_non_CHO_or_hallucinated"]))


if __name__ == "__main__":
    unittest.main()
