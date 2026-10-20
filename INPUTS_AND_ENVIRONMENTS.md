# Inputs and environments

## Hisemble core environment

Follow `Reproduction_Tutorial.md` to install the core environment. Tests use synthetic inputs; biological analyses require the data and environments below.

## scATAC preprocessing and benchmarks

`DATA_ROOT` follows the original layout:

```
ADATA/<dataset>_scBasset.pkl
ADATA/<dataset>_snap.pkl
lsi_results/ADATA/<dataset>_lsi.pkl
0.01_h5ad/<dataset>_processed.h5ad
```

The pickles hold AnnData objects with `obsm['projection']`, `obsm['X_spectral']`, and `obsm['X_lsi']`, respectively. See each script for the expected annotation field and dataset list. Some benchmark scripts instead expect LSI under `ADATA/`; use the filenames and locations expected by the selected script. Benchmark scripts write outputs under `RESULTS_ROOT` in `Hisemble_final/`, baseline-specific folders, and `downsampling/`; downstream scripts read these generated outputs. Preprocessing scripts take h5ad/FASTA manifests; genome assemblies must match each dataset. Place the reference FASTA files under `GENOME_ROOT`.

Install the core dependencies as described in `Reproduction_Tutorial.md` and any additional packages imported by the selected workflow. Install upstream encoders in separate environments; their versions and hardware-specific environments are not pinned here.

## scRNA

`workflows/scrna/Hisemble_RNA.py` and `workflows/scrna/SC3_RNA.py` read `Tabula_Sapiens/TS_<tissue>_30k_processed.h5ad` relative to the working directory and write outputs under `results/Hisemble_RNA/` and `results/SC3_RNA/`, respectively. Their `--tissues` defaults are Liver and Pancreas.

## Spatial omics

`workflows/spatial/run_Hisemble_spatial.py` runs Hisemble fusion for both MISAR-seq modalities. Place `<section>_expr_preprocessed.h5ad` and `<section>_atac_preprocessed.h5ad` under `SPATIAL_ROOT`, for sections E11_0-S1, E13_5-S1, E15_5-S1 and E18_5-S1. RNA uses PCA, GraphST and physical-coordinate graphs; ATAC uses LSI and physical-coordinate graphs. Both produce Hisemble ARI/NMI, cluster assignments and the fused graph. Inputs need `Annotation_for_Combined` and spatial coordinates. Scanpy/scikit-misc and GraphST/PyTorch are needed to compute the upstream RNA representations; ATAC and fusion from cached embeddings use the core environment. GraphST accelerator use is separate from CPU-only Hisemble fusion. The original spatial graph settings are representation k=20, physical k=8 and T=20. `run_GraphST.py` remains a separate GraphST/PCA comparison.

## Downstream analyses

`workflows/downstream/peakbed2gene/peakbed2gene.py` maps BED peaks to gene intervals extended 2 kb upstream of the transcription start site.

Set `ANNOTATION_ROOT` in `config.example.json` to your annotation directory (e.g., `/data/annotations`). Save it as `config.json` and run `prepare_workflows.py` following [Reproduction_Tutorial.md](Reproduction_Tutorial.md). Run the generated script, which contains your configured path.

Provide `gencode.vM25.annotation.gtf` in the configured annotation directory and a four-column, tab-separated BED input file (chromosome, start, end, peak identifier). Set `BED_PATH` in the generated script to your input file; its default is `cluster_0_GREAT.bed` relative to the working directory. If using another annotation filename, also update `GTF_PATH` in that script. The annotation and peaks must use the same species, genome assembly, and compatible chromosome names and coordinates. Output is written to `GO/mapped_peaks_to_genes.csv` relative to the working directory.

## Lineage trajectory analysis

`reproducibility/requirements-trajectory.txt` supplies the dependencies for the lineage trajectory analysis, including Scanpy 1.9.8; Docker installs them automatically. Supply `mca_BoneMarrow_62216_Hisemble_result.h5ad` containing the saved Hisemble consensus graph in `.obsp['connectivities']` and annotations in `.obs['cell_type']`. The script calculates B-cell diffusion pseudotime (DPT) and checks ten root choices, as reported in Supplemental Fig. S10 and Supplemental Table S8.
