# CHANGELOG

## 2026-06-04 - PR2 validation checkpoint: TSS extension and expressed-only TE labels

Branch:

* `improved-v1.1-pr2-tss-expression`

Latest commits:

* `d52210d` Add PR2 TSS extension and expressed-only TE labels
* `92b3279` Make PR2 selection expression gates group-specific

Status:

* Code pushed.
* `py_compile` passed by Codex for 02/03/10.
* Synthetic smoke test passed by Codex for group-specific 10 selection.
* Full raw-data validation pending on company workstation.

Changed files:

```text
01_pipeline/scripts/02_tss_correction.py
01_pipeline/scripts/03_map_rna_ribo_public_te.py
01_pipeline/scripts/10_select_2000_cluster_diverse_library.py
docs/VALIDATION_PLAN.md
```

Summary:

* PR2 improves biological reliability of TSS correction and public TE labels.
* PR2 keeps PR1 08 -> 07 -> 10 heavy-score connection intact.
* PR2 keeps final selected CSV/FASTA filenames unchanged.
* PR2 changes final selection from global expression gate to group-specific expression gate.

Details:

1. `02_tss_correction.py`

* Added `trim_to_tss` / `extend_to_tss` behavior.
* Upstream TSS within `max_extend` can extend annotated 5'UTR.
* Added TSS correction mode and QC summary.

2. `03_map_rna_ribo_public_te.py`

* Added `is_expressed_public`.
* Added `expression_qc_reason`.
* TE/residual/`robust_public_te_rank` are calculated only for expressed rows.
* Non-expressed rows are treated as unreliable labels rather than forced low scores.

3. `10_select_2000_cluster_diverse_library.py`

* `base_cand` uses sequence QC only.
* `evidence_cand` uses `base_cand` + expression pass + `robust_public_te_rank`.
* A/B/C/D/E use `evidence_cand`.
* F/G/H can use clean `base_cand` where appropriate.
* J fill uses `evidence_cand` first, then `base_cand`.
* Added `selection_source` reporting.
* Added `evidence_candidate_pool_after_expression_TE_QC` reporting.

Validation required:

* Full pipeline rerun with real raw data.
* Check `selected_n` remains 2000.
* Check final library length remains 50-100.
* Check `max_per_seq_cluster` remains <= 2.
* Check `heavy_ensemble_score` present in final output.
* Check TSS correction summary.
* Check expressed-only label counts.
* Check `selection_source` counts.
* Rerun PR1.1 split disjointness check after PR2.

## 2026-06-04 - PR1.1 validated: split disjointness and final diversity reports

Branch:

* `improved-v1.1-pr1.1`
* merged into `improved-v1.1`

Commit:

* `b061ee7` Add PR1.1 split disjointness and final diversity reports

Changed files:

```text
01_pipeline/scripts/09_cluster_aware_classification_benchmark.py
01_pipeline/scripts/10_select_2000_cluster_diverse_library.py
```

Added outputs:

```text
06_modeling/tables/cluster_split_disjointness_check.csv
06_modeling/tables/final_library_gene_cluster_diversity_summary.txt
```

Validation:

* `gene_seq_cluster_split` passed for all four targets.
* Gene overlap = 0.
* Sequence cluster overlap = 0.
* `pass_required_for_split = True`.
* Final library `selected_n = 2000`.
* `n_unique_seq_clusters = 1943`.
* `max_per_seq_cluster = 2`.
* `n_unique_genes = 1850`.
* `max_per_gene = 4`.
* `heavy_ensemble_score` non-null = 2000 / 2000.
* length range = 50-100.

Interpretation:

* Final library is sequence-cluster-diverse.
* `seq_cluster_id` is a 5'UTR sequence-similarity cluster, not a gene cluster.
* Gene-level diversity is high but not hard-capped.

## 2026-06-03 - PR1 validated: cluster-aware heavy modeling connected to final selection

Branch:

* `improved-v1.1`

Commit:

* `4a29c08` Connect cluster-aware heavy modeling to final selection

Changed files:

```text
01_pipeline/scripts/run_00_full_final_pipeline.py
01_pipeline/scripts/07_heavy_rnafold_kmer6_automl.py
01_pipeline/scripts/10_select_2000_cluster_diverse_library.py
```

Changes:

* Run order changed so 08 Jaccard clustering runs before 07 heavy modeling.
* 07 now reads `tss_corrected_5utr_with_seq_clusters.csv` as primary input.
* 07 adds `seq_cluster_split` and `gene_seq_cluster_split`.
* 07 outputs `tss_corrected_5utr_with_seq_clusters_and_heavy_scores.csv`.
* 10 prioritizes heavy-score integrated input table.
* Final CSV/FASTA filenames remain unchanged.

Validation:

* `py_compile` passed for modified scripts.
* 07 created heavy-score integrated table.
* 09 benchmark completed.
* 10 selected final 2000-member library.
* `heavy_ensemble_score` non-null = 2000 / 2000.
* length range = 50-100.
* max per `seq_cluster_id` = 2.

Key benchmark result:

* `robust_public_te_rank` gene_seq_cluster_split RandomForest: AUC 0.671 / AP 0.582
* `protein_residual_rank` gene_seq_cluster_split ExtraTrees: AUC 0.696 / AP 0.630

Interpretation:

* Sequence features retain predictive signal even under strict gene + sequence-cluster split.
* `protein_residual_rank` is useful auxiliary evidence.
* `protein_abundance_rank` is weaker under strict split and should remain auxiliary.

## 2026-06-03 - Raw Data And LFS Handling

Summary:

* Raw public datasets were added using Git LFS where needed.
* Large NCBI mapping resources were reduced or filtered where possible.
* Code and raw data can be retrieved for reproducible local/company workstation runs.

Notes:

* Generated outputs should generally not be committed.
* Raw data may be kept via Git LFS or local/company storage depending on size/security constraints.
