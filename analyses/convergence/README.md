# Convergence analysis

The runner uses fixed initial affinity matrices, synchronous updates, top-q sparsification, and lazy row normalization. `vendor_core/` contains the numerical routines used by this runner; `source_manifest.json` records their checksums, which are checked before a run. The alternative entry point `workflows/scatac/hisemble/Hisemble_convergence.py` uses the same fixed-W propagation.

Run the commands below from `analyses/convergence/`.

Run one full dataset in a fresh output directory:

```sh
python run_core_convergence.py --data-root /path/to/ensemble --output-dir /path/to/run/mca_Cerebellum_62216 --dataset mca_Cerebellum_62216 --iterations 50 --chunk-size 2048 --seed 42
python export_original_metrics.py --input-dir /path/to/run --output-dir /path/to/run/original_metrics --iterations 50
```

Run all 24 datasets:

```sh
python launch_all.py --data-root /path/to/ensemble --run-root /path/to/new_run --iterations 50 --workers 4 --threads 8
```

The batch launcher exports the paper metrics after successful numerical runs.

Paper metric definitions:

- `Relative Change = sum(abs(P_next-P_old))/(N*q)`, where k=15 and q=2k=30.
- `Delta_Q` is recorded before each P update. Record 1 is 0; record t>=2 uses diagnostic `Q_delta_fro` at iteration t-1.

The full numerical runner retains additional diagnostics in its raw CSV. `export_original_metrics.py` exports the paper's P/Q metrics, checks trajectory completeness and converts the Q timing. Resource measurements include diagnostics and do not replace benchmark measurements.

All views use identical cells and row order; scBasset and SnapATAC2 use Gaussian affinity and LSI uses nonnegative cosine affinity. The runner checks embeddings, saves consensus_T20.npz, and records partial/complete status without silently resuming incomplete runs. `datasets.json` lists the 24 benchmark dataset IDs. Default diagnostic length is 100; specify `--iterations 50` to reproduce the displayed range. Biological input pickles require AnnData/pandas versions compatible with those used to save them.

Run synthetic-fixture checks with `python test_runner_equivalence.py` to compare the diagnostic runner with the numerical core.
