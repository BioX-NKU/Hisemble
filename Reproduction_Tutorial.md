# Hisemble reproduction tutorial

Use Docker on macOS or Linux, or through WSL2 on Windows. Conda installation is also described below. The Hisemble core runs on CPUs.

## 1. Obtain the code

Obtain the source through GitHub's **Code** menu, then open a terminal at the repository root. Run all commands below from that directory.

## 2. Start Docker and build

Open Docker Desktop or OrbStack and wait for its engine to start. Then check:

```bash
docker version
docker build -f reproducibility/Dockerfile -t hisemble:reproduce .
```

The image includes Hisemble and dependencies for the lineage trajectory analysis in Supplemental Fig. S10. The first build downloads the required packages.

### Install without Docker

From the code directory, create a Conda environment and install the package:

```bash
conda env create -f reproducibility/environment-core.yml
conda activate hisemble-core
python -m pip install -r reproducibility/requirements-trajectory.txt
python -m pip install --no-deps -e .
```

This installs Hisemble and the dependencies for the lineage trajectory analysis in Supplemental Fig. S10. To run the example locally:

```bash
python examples/quickstart_synthetic.py --output-dir example_output
```

The commands in the following sections use Docker.

## 3. Run the checks

```bash
docker run --rm hisemble:reproduce \
  python -m unittest discover -s /opt/hisemble/tests -p 'test_*.py'
```

Expected outcome: `OK`, covering input validation, view-order invariance, the AnnData interface, and spatial RNA/ATAC fusion with synthetic inputs.

For convergence checks, examples, seed consistency, and trajectory-analysis checks, run:

```bash
docker run --rm hisemble:reproduce \
  python /opt/hisemble/reproducibility/verify_installation.py
```

Expected final line begins `PASS: all packaged tests`. These checks use synthetic inputs to exercise the core, examples, and numerical entry points.

## 4. Run Hisemble and save an AnnData result

```bash
mkdir -p example_output
docker run --rm \
  -v "$PWD/example_output:/output" \
  hisemble:reproduce \
  python /opt/hisemble/examples/quickstart_synthetic.py --output-dir /output
```

Expected output:

```text
cells=180
clusters=4
graph_shape=(180, 180)
graph_nnz=7372
```

The example generates `example_output/synthetic_hisemble.h5ad`, which stores the input embeddings in `.obsm`, the consensus graph in `.obsp['connectivities']`, distances in `.obsp['distances']`, Leiden labels in `.obs['hisemble_clusters']`, and neighbor metadata in `.uns['neighbors']`.

## 5. Re-run from a saved input file

This exercises the command used for your own aligned embeddings:

```bash
docker run --rm \
  -v "$PWD/example_output:/data" \
  hisemble:reproduce \
  python /opt/hisemble/examples/run_hisemble_from_h5ad.py \
    --input /data/synthetic_hisemble.h5ad \
    --output /data/synthetic_hisemble_rerun.h5ad \
    --keys view_a,view_b --metrics gaussian,gaussian \
    --k 15 --iterations 20 --resolution 0.5 --seed 42
```

For biological data, substitute your `.h5ad` filename and actual `.obsm` keys. All embeddings must describe the same cells in the same row order. For scATAC inputs, typical keys are `scBasset,SnapATAC2,LSI`, with `gaussian,gaussian,cosine` metrics. Keep each view paired with its intended metric.

NaN/Inf or invalid dimensions stop validation before graph construction. Finite extreme row norms produce a diagnostic warning by default. No cells are automatically deleted or imputed.

## 6. Reproduce supplementary analyses

Provide the required data and precomputed embeddings, then run the analysis commands below.

Use this aligned AnnData-pickle layout:

```text
ADATA/<dataset>_scBasset.pkl        # obsm['projection']
ADATA/<dataset>_snap.pkl           # obsm['X_spectral']; obs['cell_type']
lsi_results/ADATA/<dataset>_lsi.pkl # obsm['X_lsi']
```

All three files for a dataset must have identical `obs_names` order. Set your local path once:

```bash
HISEMBLE_INPUT_DIR=/absolute/path/to/hisemble_inputs
mkdir -p example_output/analyses
```

### Ten-seed clustering stability

```bash
docker run --rm \
  -v "$HISEMBLE_INPUT_DIR:/data:ro" \
  -v "$PWD/example_output/analyses:/results" \
  hisemble:reproduce \
  python /opt/hisemble/analyses/seed_stability.py \
    --data-root /data --output-dir /results/seed_stability
```

By default, the script runs ten seeds on eight datasets (Supplemental Table S6). For a single-dataset run, append `--datasets mca_Cerebellum_62216`. Each seed repeats the resolution search with fixed embeddings and graph.

For parallel seeds inside Docker, substitute `seed_stability_parallel.py` and append `--seed-workers 2` (or more workers if RAM allows). Each dataset builds its consensus graph once.

### Degraded-LSI analysis

```bash
docker run --rm \
  -v "$HISEMBLE_INPUT_DIR:/data:ro" \
  -v "$PWD/example_output/analyses:/results" \
  hisemble:reproduce \
  python /opt/hisemble/analyses/degraded_view_analysis.py \
    --data-root /data --output-dir /results/degraded_view
```

Defaults reproduce the 66-run design on Mouse Cortex, FCA Cerebellum, and MCA Cerebellum. Use `--datasets mca_Cerebellum_62216` for a smaller first run.

### Components and controlled scaling

The component entry point is `component_ablation.py`, with `--data-root`, `--dataset`, `--variant`, and `--output-dir`. Variants are `full`, `no_diffusion`, and `no_dynamic_sparsification`; default k=15, T=20, seed=42. The full component runner is `run_component_ablation.sh`, configured for complete Mouse Cortex, FCA Intestine, and MCA Cerebellum with a four-hour per-run cap.

Retain run logs to distinguish timeouts from failures. Dense variants may require substantial memory.

`run_component_scaling.sh` runs the MCA Cerebellum downsampling analysis. `lsi_neighbor_analysis.py --data-root /data --dataset mca_Cerebellum_62216 --output-dir /results/lsi` evaluates LSI neighborhood complementarity. Run these through the same input/output mounts used above.

### Lineage trajectory analysis

Set a directory containing `mca_BoneMarrow_62216_Hisemble_result.h5ad` with `.obsp['connectivities']` and `.obs['cell_type']`.
```bash
HISEMBLE_GRAPH_DIR=/absolute/path/to/saved_hisemble_graphs
docker run --rm \
  -v "$HISEMBLE_GRAPH_DIR:/graphs:ro" \
  -v "$PWD/example_output/analyses:/results" \
  hisemble:reproduce \
  python /opt/hisemble/analyses/lineage_trajectory_analysis.py \
    --data-root /graphs --output-dir /results/lineage_trajectory
```

The script computes diffusion pseudotime (DPT) and checks ten root choices (Supplemental Fig. S10). Results are saved to the directory specified by `--output-dir`. Compare the summaries with Supplemental Table S8.

## 7. Prepare modality-specific workflows

Choose a workflow using [ANALYSIS_MAP.md](ANALYSIS_MAP.md).

Save your paths as `config.json`, using `config.example.json` as a template, then prepare a new working directory:

```bash
python prepare_workflows.py --config /absolute/path/config.json \
  --output /absolute/path/configured_workflows
export PYTHONPATH="/absolute/path/configured_workflows${PYTHONPATH:+:$PYTHONPATH}"
```

`WORKFLOW_ROOT` in the configuration should point to this same output directory. Launch each configured script in the input directory described in `INPUTS_AND_ENVIRONMENTS.md`. The shared helper is `configured_workflows/utils.py`.

## Troubleshooting

- Docker engine unavailable: start Docker Desktop or OrbStack and retry `docker version`.
- Missing embedding key: inspect your `.obsm` keys and supply the exact names.
- NaN/Inf or incompatible cells: correct or explicitly align the upstream inputs before fusion.
- Memory pressure: begin with a small dataset and fewer workers. Atlas-scale and dense-ablation runs may require more RAM.
- Run all commands from one extracted package; do not mix files from different versions.

## Spatial RNA and ATAC fusion

Run `workflows/spatial/run_Hisemble_spatial.py` with `--data-root`, `--output-dir` and `--mode RNA`, `ATAC` or `both`. This is the Hisemble spatial fusion entry point; `run_GraphST.py` only runs a comparison.

## Convergence analysis

See [the convergence guide](analyses/convergence/README.md) for commands and metric definitions. Use `--iterations 50` for the range shown in the figures.
