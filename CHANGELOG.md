# 5UTR Engineering CHANGELOG

This file records only project-level decisions and validated changes. Long exploratory chat history should not be copied here.

---

## 2026-06-08

### Corrected

- Clarified that the current official baseline remains `v1.2`.
- Clarified that PR3 / PR3-2 validation should be folded into the `v1.2` baseline when accepted.
- Removed `improved-v1.3 candidate` as a current baseline concept.
- Clarified that PR4 was not executed and should not be treated as an active or abandoned development branch.
- Clarified that the next improvement should be named as a `v1.2` follow-up branch/work item, not as standalone `PR5`.

### PR3 / PR3-2 Validation Record

PR3 introduced the uAUG0 production-selection validation gate. PR3-2 is the workstation-validated PR3 follow-up.

Validation results:

```text
selected_n=2000
uaug_positive_n=0
n_unique_seq_clusters=1937
max_per_seq_cluster=2
max_per_gene=3
mean_heavy_ensemble_score=0.5811515
mean_robust_public_te_rank=0.642380
```

Decision:

```text
Treat PR3-2 as accepted PR3 evidence for v1.2 maintenance.
Do not call the current state improved-v1.3.
Do not advance to PR4/PR5 naming until the accepted PR3/PR3-2 state is folded into v1.2.
```

### Branch Cleanup Note

The remote branch `pr5-diversity-qc` was created from a documentation-only commit and is not the official next baseline. It should be ignored or replaced by a properly named `v1.2` follow-up branch after the next work item is defined.

---

## 2026-06-04

### Added

- Added `MASTER.md` as the lightweight project consensus document.
- Established GitHub `main` + `MASTER.md` as the source of truth for future ChatGPT workspaces.
- Added workspace rotation rule: use `5UTR_WORKSPACE_01`, `5UTR_WORKSPACE_02`, etc. and update `MASTER.md` before moving to the next workspace.

### Confirmed

- Repository is the final, team-shareable CHO 5' UTR candidate discovery pipeline.
- Final pipeline starts from raw/reference inputs only.
- Official run entry points:

```bash
bash RUN_FINAL_MAIN.sh
python 01_pipeline/scripts/run_00_full_final_pipeline.py
```

### Pipeline Consensus

```text
CHO genome/GFF annotation
→ TSS-atlas correction
→ RNA/Ribo public TE labeling
→ Heffner proteomics mapping
→ multi-omics label generation
→ RNAfold/k-mer6/tree2000 heavy modeling
→ Jaccard sequence cluster QC
→ cluster-aware classification benchmark
→ final 2,000 cluster-diverse 50-100 bp 5' UTR library
```

### Active Development Note

- Step 03 public TE mapping / robust TE mapping should preserve or reconcile the recent auto sample assignment rule:

```text
rna_day3  -> s01-s03
ribo_day3 -> s04-s06
rna_day6  -> s07-s09
ribo_day6 -> s10-s12
```

### Final Outputs

```text
07_library_design/tables/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv
07_library_design/fasta/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.fasta
```
