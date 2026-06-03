# 5UTR Engineering MASTER

## 1. Current Project State

Repository:

* `SBLGENEALL/5UTR_Engineering`

Stable branch:

* `main`
* Initial stable numbered pipeline release.

Current working baseline:

* `improved-v1.1`
* PR1 + PR1.1 validated.
* The 08 Jaccard sequence clustering -> 07 heavy RNAfold/k-mer modeling -> 09 cluster-aware benchmark -> 10 final library selection connection has been validated.

Current experimental branch:

* `improved-v1.1-pr2-tss-expression`
* Purpose: PR2 TSS correction + expressed-only public TE labels + group-specific final selection expression gate.
* Status: code pushed and lightweight checks passed, but full raw-data validation is still pending on the company workstation.

Current operating rule:

* Other chats should treat this `MASTER.md` as the source of truth.
* Historical chat context is secondary.
* If a future chat becomes slow, start a new chat and ask it to read `MASTER.md` and `CHANGELOG.md` first.

## 2. Current Validated Baseline: improved-v1.1

PR1 goal:

* Restore correct pipeline connection among 08 sequence clustering, 07 heavy model scoring, and 10 final 2000 library selection.

PR1 validated outputs:

```text
04_te_labeling/tables/tss_corrected_5utr_with_seq_clusters.csv
04_te_labeling/tables/tss_corrected_5utr_with_seq_clusters_and_heavy_scores.csv
07_library_design/tables/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv
07_library_design/fasta/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.fasta
```

PR1 validation result:

* `selected_n`: 2000
* `heavy_ensemble_score` non-null: 2000 / 2000
* length range: 50-100
* max per `seq_cluster_id`: 2

PR1.1 goal:

* Add validation guardrails and reporting without changing modeling or selection logic.

PR1.1 added outputs:

```text
06_modeling/tables/cluster_split_disjointness_check.csv
06_modeling/tables/final_library_gene_cluster_diversity_summary.txt
```

PR1.1 validation result:

* `gene_seq_cluster_split` passed for all four targets.
* `gene_overlap_count = 0`.
* `seq_cluster_overlap_count = 0`.
* `pass_required_for_split = True`.
* `selected_n`: 2000
* `n_unique_seq_clusters`: 1943
* `max_per_seq_cluster`: 2
* `n_unique_genes`: 1850
* `max_per_gene`: 4

Interpretation:

* `seq_cluster_id` is a 5'UTR sequence-similarity cluster, not a gene cluster.
* Final diversity is enforced on sequence clusters.
* Gene-level diversity is reported but not hard-capped.
* The same gene can appear multiple times if distinct transcript/UTR isoforms fall into different sequence clusters.

Recommended wording:

> The final 2,000-member library spans 1,943 5'UTR sequence-similarity clusters, with no more than 2 candidates per sequence cluster. It is drawn from 1,850 unique source genes, with a maximum of 4 candidates per gene. Diversity is enforced on sequence clusters; gene-level spread is reported but not capped.

## 3. Pipeline Order

Expected full pipeline order:

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

## 4. Current PR2 Branch

Branch:

* `improved-v1.1-pr2-tss-expression`

Latest PR2 commits:

* `d52210d` Add PR2 TSS extension and expressed-only TE labels
* `92b3279` Make PR2 selection expression gates group-specific

PR2 goal:

* Improve TSS correction and public TE label reliability.
* Preserve PR1 08 -> 07 -> 10 heavy-score connection.
* Keep final selected CSV/FASTA filenames unchanged.
* Improve final selection logic so expression gating is applied to evidence-supported groups but does not remove clean exploratory/diversity/negative-control candidates unnecessarily.

Changed files:

```text
01_pipeline/scripts/02_tss_correction.py
01_pipeline/scripts/03_map_rna_ribo_public_te.py
01_pipeline/scripts/10_select_2000_cluster_diverse_library.py
docs/VALIDATION_PLAN.md
```

PR2 design summary:

1. `02_tss_correction.py`

* TSS inside annotated UTR: `trim_to_tss`.
* Confident upstream TSS within `max_extend`: `extend_to_tss`.
* Adds TSS correction mode and QC summary.
* TSS extension should be conservative because the final library uses 50-100 bp windows.

2. `03_map_rna_ribo_public_te.py`

* Adds `is_expressed_public`.
* Adds `expression_qc_reason`.
* Computes TE/residual/`robust_public_te_rank` only for expressed rows.
* Non-expressed or unreliable labels should become `NaN`, not zero-scored.

3. `10_select_2000_cluster_diverse_library.py`

* Group-specific expression gate.
* `base_cand`: sequence QC only.
* `evidence_cand`: `base_cand` + `is_expressed_public` + non-null `robust_public_te_rank`.
* A/B/C/D/E selected from `evidence_cand`.
* F/G selected from `base_cand` to preserve exploratory/diversity candidates.
* H negative controls do not require the expression gate.
* J fill uses `evidence_cand` first, then `base_cand` if needed.
* Summary reports `candidate_pool_after_QC`, `evidence_candidate_pool_after_expression_TE_QC`, `selection_source` counts, and expression gate coverage.

## 5. PR2 Company Workstation Validation Plan

PR2 must be validated by full raw-data rerun because 02 can alter UTR boundaries and therefore affects:

* 08 sequence clustering
* 07 heavy model score
* 09 benchmark
* 10 final library

Company workstation procedure:

1. Download ZIP from GitHub branch:

```text
improved-v1.1-pr2-tss-expression
```

2. Extract into a fresh folder, not over an old PR1 output folder.

3. Place real raw data under:

```text
00_raw_data/
```

4. Run:

```bash
python 01_pipeline/scripts/00_check_inputs.py
```

5. If input check passes, run:

```bash
python 01_pipeline/scripts/run_00_full_final_pipeline.py > run_pr2_full.log 2>&1
```

6. After completion, inspect:

```text
run_pr2_full.log
03_tss_correction/qc/tss_correction_summary.txt
04_te_labeling/qc/robust_public_te_mapping_summary.txt
07_library_design/qc/selected_2000_50_100bp_cluster_diverse_evidence_balanced_summary.txt
06_modeling/tables/final_library_gene_cluster_diversity_summary.txt
06_modeling/tables/cluster_split_disjointness_check.csv
```

PR2 first-pass acceptance criteria:

* Input check passes.
* Full pipeline completes without error.
* `selected_n` remains 2000.
* Final library length range remains 50-100.
* `max_per_seq_cluster` remains <= 2.
* `heavy_ensemble_score` exists in final selected library.
* `candidate_pool_after_QC` is not too small.
* `evidence_candidate_pool_after_expression_TE_QC` is sufficient for A-E evidence-supported groups.
* `selection_source` counts show intended mixture of `evidence_cand` and `base_cand`-derived groups.
* H negative controls are not unintentionally removed by expression gate.
* `gene_seq_cluster_split` still has `gene_overlap_count = 0` and `seq_cluster_overlap_count = 0`.
* PR2 output numbers may differ from PR1.1; this is expected because TSS correction and public TE label reliability changed.

Minimum result values to report back from company workstation:

* `selected_n`
* `candidate_pool_after_QC`
* `evidence_candidate_pool_after_expression_TE_QC`
* `selection_source` counts
* `heavy_ensemble_score` non-null count
* length min/max
* `max_per_seq_cluster`
* `n_unique_seq_clusters`
* `n_unique_genes`
* `max_per_gene`
* `is_expressed_public` count
* `expression_qc_reason` counts
* TSS correction mode counts
* any error/warning from `run_pr2_full.log`

## 6. Next Planned Work

PR3:

* Construct-level screening
* uAUG/uORF strict screening
* cryptic splice donor/acceptor screening
* polyA-like motif screening
* restriction enzyme / assembly site screening
* mAb CDS junction context
* classify checks into hard fail, soft warning, and report-only

PR4:

* Team release documentation
* README update
* final validated tag
* possible promotion of `improved-v1.1` or later validated branch to `main`

## 7. Git Strategy

`main`:

* Stable release only.

`improved-v1.1`:

* Next stable candidate.
* PR1 + PR1.1 validated.

`improved-v1.1-pr2-tss-expression`:

* PR2 experimental validation branch.
* Full raw-data validation pending.

Merge policy:

* Do not merge PR2 into `improved-v1.1` until company workstation validation passes.
* Do not merge `improved-v1.1` into `main` until PR2/PR3 validation and documentation are complete.
* Generated outputs should not be committed unless explicitly chosen as small release artifacts.
