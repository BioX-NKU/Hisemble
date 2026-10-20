from __future__ import annotations

import contextlib
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    'spatial_workflow', ROOT / 'workflows/spatial/run_Hisemble_spatial.py'
)
spatial = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(spatial)


class SpatialWorkflowTests(unittest.TestCase):
    def fixture(self):
        rng = np.random.default_rng(24)
        data = ad.AnnData(rng.poisson(2, size=(64, 60)).astype(np.float32))
        data.obs_names = [f'cell_{index}' for index in range(data.n_obs)]
        data.obs['x'] = np.tile(np.arange(8), 8)
        data.obs['y'] = np.repeat(np.arange(8), 8)
        data.obs['Annotation_for_Combined'] = np.repeat(['a', 'b', 'c', 'd'], 16)
        data.obsm['X_pca'] = rng.normal(size=(64, 30))
        data.obsm['emb'] = rng.normal(size=(64, 12))
        return data

    def test_rna_has_three_independent_view_graphs(self):
        data = self.fixture()
        embeddings = spatial.prepare_embeddings(data, 'RNA')
        self.assertEqual(list(embeddings), ['PCA', 'GraphST', 'Spatial'])
        with contextlib.redirect_stdout(io.StringIO()):
            graphs, fused = spatial.fuse_views(embeddings, data.obs_names, 5, 3, 2)
        self.assertEqual(list(graphs), ['PCA', 'GraphST', 'Spatial'])
        self.assertEqual(fused.shape, (64, 64))
        self.assertGreater(np.abs((graphs['Spatial'] - graphs['PCA']).data).sum(), 0)
        self.assertTrue(np.isfinite(fused.data).all())
        self.assertLessEqual(abs((fused - fused.T).data).max(initial=0), 1e-12)

    def test_atac_computes_lsi_and_two_view_fusion(self):
        data = self.fixture()
        embeddings = spatial.prepare_embeddings(data, 'ATAC')
        self.assertEqual(list(embeddings), ['Spatial', 'LSI'])
        self.assertEqual(embeddings['LSI'].shape, (64, 29))
        with contextlib.redirect_stdout(io.StringIO()):
            _, fused = spatial.fuse_views(embeddings, data.obs_names, 5, 3, 2)
            _, initial = spatial.fuse_views(embeddings, data.obs_names, 5, 3, 0)
        self.assertGreater(np.abs((fused - initial).data).sum(), 0)

    def test_invalid_coordinates_rejected_before_graph_construction(self):
        data = self.fixture()
        data.obs.loc[data.obs_names[0], 'x'] = np.nan
        embeddings = spatial.prepare_embeddings(data, 'RNA')
        with self.assertRaisesRegex(ValueError, 'Spatial.*NaN/Inf'):
            spatial.fuse_views(embeddings, data.obs_names, 5, 3, 2)

    def test_duplicate_cell_identifiers_rejected(self):
        data = self.fixture()
        data.obs_names = ['duplicate'] * data.n_obs
        with self.assertRaisesRegex(ValueError, 'identifiers must be unique'):
            spatial.prepare_embeddings(data, 'ATAC')

    def test_cli_generates_metrics_and_fused_graph_for_both_modalities(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            self.fixture().write_h5ad(path / 'E11_0-S1_expr_preprocessed.h5ad')
            self.fixture().write_h5ad(path / 'E11_0-S1_atac_preprocessed.h5ad')
            output = path / 'output'
            with contextlib.redirect_stdout(io.StringIO()):
                spatial.main([
                    '--data-root', str(path), '--output-dir', str(output),
                    '--sections', 'E11_0-S1', '--k-neighbors', '5',
                    '--k-spatial', '3', '--iterations', '2',
                ])
            metrics = pd.read_csv(output / 'spatial_metrics.csv')
            self.assertEqual(len(metrics), 5)
            self.assertEqual(set(metrics.modality), {'RNA', 'ATAC'})
            hisemble = metrics[metrics.method == 'Hisemble']
            self.assertEqual(hisemble.set_index('modality').n_views.to_dict(), {'RNA': 3, 'ATAC': 2})
            self.assertTrue(np.isfinite(metrics[['ARI', 'NMI']].to_numpy()).all())
            for mode, views in [('RNA', 3), ('ATAC', 2)]:
                result = ad.read_h5ad(output / f'E11_0-S1_{mode}_hisemble.h5ad')
                self.assertEqual(result.n_obs, 64)
                self.assertIn('Hisemble', result.obs)
                self.assertEqual(len(result.uns['hisemble_spatial']['views']), views)
                self.assertEqual(result.uns['hisemble_spatial']['resolution_fallback'], 'last_tested')
                self.assertIn('spatial_connectivities', result.obsp)
                self.assertIn('hisemble_connectivities', result.obsp)

    def test_missing_sections_fail_before_computation(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(FileNotFoundError, 'Missing input sections'):
                spatial.main(['--data-root', directory, '--output-dir', directory])


if __name__ == '__main__':
    unittest.main()
