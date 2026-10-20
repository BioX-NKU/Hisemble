# Computational analysis entry points

| Reported analysis | Entry point(s) | Principal output |
|---|---|---|
| scATAC Hisemble benchmark | `workflows/scatac/hisemble/run_Hisemble.py` | clustered AnnData and benchmark metrics |
| LSI ablation | `workflows/scatac/hisemble/Hisemble_ablation_lsi.py` | full/no-LSI metrics |
| Neighborhood sensitivity | `workflows/scatac/hisemble/Hisemble_K_Sensitive.py` | metrics across k |
| Convergence (fixed initial W) | `analyses/convergence/run_core_convergence.py`, `export_original_metrics.py`; compatible workflow `workflows/scatac/hisemble/Hisemble_convergence.py` | per-iteration convergence metrics and iteration-20 consensus graph |
| SC3 and CHAI baselines | `workflows/scatac/baselines/` | baseline metrics |
| Scaling/downsampling | `workflows/scatac/downsampling/` | runtime, memory, and clustering metrics |
| scBasset, SnapATAC2, and LSI representations | `workflows/scatac/representations/` | precomputed cell embeddings |
| scATAC preprocessing | `workflows/scatac/preprocess/` | processed AnnData inputs |
| scRNA benchmark (Liver and Pancreas; Fig. 5A--B) | `workflows/scrna/` (`Hisemble_RNA.py`, `SC3_RNA.py`) | clustering metrics and labels |
| Spatial RNA fusion (Fig. 5C) | `workflows/spatial/run_Hisemble_spatial.py --mode RNA` | Hisemble, GraphST and PCA ARI/NMI; fused graph and labels |
| Spatial ATAC fusion (Fig. 5D--E) | `workflows/spatial/run_Hisemble_spatial.py --mode ATAC` | Hisemble and LSI ARI/NMI; fused graph and labels |
| GraphST/PCA comparison | `workflows/spatial/run_GraphST.py` | GraphST and PCA comparison |
| Peak-to-gene downstream analysis | `workflows/downstream/peakbed2gene/peakbed2gene.py` | peak-to-gene result table |
| LSI complementarity | `workflows/scatac/hisemble/Hisemble_ablation_lsi.py` and `analyses/lsi_neighbor_analysis.py` | full/no-LSI metrics and neighborhood summaries |
| Component ablation/scaling | `analyses/component_ablation.py`, `analyses/run_component_ablation.sh`, `analyses/run_component_scaling.sh` | CSV/JSON measurements |
| Seed stability | `analyses/seed_stability.py` | — |
| Degraded-view robustness | `analyses/degraded_view_analysis.py` | — |
| Graph-based lineage analysis | `analyses/lineage_trajectory_analysis.py` | MCA Bone Marrow B-cell pseudotime and ten-root sensitivity |
