# 5UTR Engineering MASTER

## Current Source Of Truth

Repository:

* SBLGENEALL/5UTR_Engineering

Main stable branch:

* `main`: initial stable numbered pipeline release

Current working baseline:

* `improved-v1.1`
* PR1 and PR1.1 have been validated and reflected.
* The 08 Jaccard sequence clustering -> 07 heavy RNAfold/k-mer modeling -> 09 cluster-aware benchmark -> 10 final library selection connection has been validated.

Current experimental branch:

* `improved-v1.1-pr2-tss-expression`
* PR2: TSS extension + expressed-only TE/residual labels + group-specific expression gate
* Status: code pushed, `py_compile` and synthetic smoke test passed by Codex, but full raw-data validation is still pending on the company workstation.

## Pipeline Order

```text
00_check_inputs.py
01_build_utr_database.py
02_tss_correction.py
03_map_rna_ribo_public_te.py
04_preprocess_heffner_proteomics.py
05_integrate_proteomics_multiomics.py
06_plot_multiomics_distributions.py
08_jaccard_sequence_cluster_qc.py
07_heavy_rnafold_kmer6_automl.py
09_cluster_aware_classification_benchmark.py
10_select_2000_cluster_diverse_library.py
```

Important:

* 08 must run before 07 because 07 uses `seq_cluster_id`.
* 07 generates `heavy_ensemble_score` and the integrated heavy-score table.
* 10 should preferentially read:

```text
04_te_labeling/tables/tss_corrected_5utr_with_seq_clusters_and_heavy_scores.csv
```

## PR1 Validated Result

Goal:

* Restore the correct connection among 08 clustering, 07 heavy model scoring, and 10 final selection.

Validated outputs:

```text
04_te_labeling/tables/tss_corrected_5utr_with_seq_clusters.csv
04_te_labeling/tables/tss_corrected_5utr_with_seq_clusters_and_heavy_scores.csv
07_library_design/tables/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv
07_library_design/fasta/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.fasta
```

Validation result:

* `selected_n`: 2000
* `heavy_ensemble_score` non-null: 2000 / 2000
* length range: 50-100
* max per `seq_cluster_id`: 2

## PR1.1 Validated Result

Goal:

* Add verification reports without changing modeling or selection logic.

Added outputs:

```text
06_modeling/tables/cluster_split_disjointness_check.csv
06_modeling/tables/final_library_gene_cluster_diversity_summary.txt
```

Key validation:

* `gene_seq_cluster_split` passed for all four targets.
* `gene_overlap_count = 0` and `seq_cluster_overlap_count = 0` in `gene_seq_cluster_split`.
* `pass_required_for_split = True` for all required split checks.

Final library diversity:

* `selected_n`: 2000
* `n_unique_seq_clusters`: 1943
* `max_per_seq_cluster`: 2
* `n_unique_genes`: 1850
* `max_per_gene`: 4

Interpretation:

* `seq_cluster_id` is a 5'UTR sequence-similarity cluster, not a gene cluster.
* Diversity is enforced by sequence cluster.
* Gene-level diversity is reported but not hard-capped.
* The same gene can appear multiple times if different transcript or UTR isoforms belong to different sequence clusters.

Recommended wording:

> The final 2,000-member library spans 1,943 5'UTR sequence-similarity clusters, with no more than 2 candidates per sequence cluster. It is drawn from 1,850 unique source genes, with a maximum of 4 candidates per gene. Diversity is enforced on sequence clusters; gene-level spread is reported but not capped.

## PR2 Current Design

Branch:

* `improved-v1.1-pr2-tss-expression`

Goal:

* Improve TSS correction and public TE label reliability.
* Keep the PR1 08 -> 07 -> 10 connection intact.
* Keep final selected CSV/FASTA filenames unchanged.

Changed files:

```text
01_pipeline/scripts/02_tss_correction.py
01_pipeline/scripts/03_map_rna_ribo_public_te.py
01_pipeline/scripts/10_select_2000_cluster_diverse_library.py
docs/VALIDATION_PLAN.md
```

PR2 changes:

1. `02_tss_correction.py`

* TSS inside annotated UTR: `trim_to_tss`
* Upstream TSS within `max_extend`: `extend_to_tss`
* Adds TSS correction mode and QC summary.
* Uses a conservative default `max_extend`, with CLI override for exploration.

2. `03_map_rna_ribo_public_te.py`

* Adds `is_expressed_public`.
* Adds `expression_qc_reason`.
* Computes TE/residual/`robust_public_te_rank` only for expressed rows.
* Non-expressed or unreliable labels remain `NaN`, not zero-scored.

3. `10_select_2000_cluster_diverse_library.py`

* Uses group-specific expression gating.
* `base_cand`: sequence QC only.
* `evidence_cand`: `base_cand` + `is_expressed_public` + non-null `robust_public_te_rank`.
* A/B/C/D/E are selected from `evidence_cand`.
* F/G are selected from `base_cand` to preserve exploratory/diversity candidates.
* H negative controls do not require the expression gate.
* J fill uses `evidence_cand` first, then `base_cand` if needed.
* Summary reports `candidate_pool_after_QC`, `evidence_candidate_pool_after_expression_TE_QC`, and `selection_source` counts.

## PR2 Full Validation Plan

PR2 must be validated by a full raw-data pipeline rerun because 02 can change UTR sequences, which can affect 08 clustering, 07 heavy model, 09 benchmark, and 10 final library selection.

Run on the company workstation:

```bash
python 01_pipeline/scripts/00_check_inputs.py
python 01_pipeline/scripts/run_00_full_final_pipeline.py > run_pr2_full.log 2>&1
```

After the full run, check:

```text
03_tss_correction/qc/tss_correction_summary.txt
04_te_labeling/qc/robust_public_te_mapping_summary.txt
07_library_design/qc/selected_2000_50_100bp_cluster_diverse_evidence_balanced_summary.txt
06_modeling/tables/final_library_gene_cluster_diversity_summary.txt
06_modeling/tables/cluster_split_disjointness_check.csv
```

PR2 acceptance criteria:

* Input check passes.
* Full pipeline completes without error.
* `selected_n` remains 2000.
* Final length range remains 50-100.
* `max_per_seq_cluster` remains <= 2.
* `heavy_ensemble_score` remains present in the final library.
* `evidence_candidate_pool_after_expression_TE_QC` is not too small.
* `selection_source` counts show A-E from `evidence_cand`, while F/G/H and fill can use `base_cand` as intended.
* `gene_seq_cluster_split` still has `gene_overlap_count = 0` and `seq_cluster_overlap_count = 0`.
* PR2 generated numbers may differ from PR1/PR1.1 because TSS correction and label reliability changed.

## Next Planned Work

PR3:

* Construct-level screening
* uAUG/uORF strict screening
* cryptic splice donor/acceptor screening
* polyA-like motif screening
* restriction enzyme / assembly site screening
* mAb CDS junction context screening

PR4:

* Team release documentation
* README update
* final validated tag
* possible promotion of `improved-v1.1` to `main` after PR2/PR3 validation
