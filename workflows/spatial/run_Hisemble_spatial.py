#!/usr/bin/env python3
"""Hisemble fusion of MISAR-seq spatial RNA and spatial ATAC representations."""
from __future__ import annotations

import argparse
from pathlib import Path
import random

import anndata as ad
import igraph as ig
import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from hisemble.cluster import find_best_matrix
from hisemble.core import hisemble_fusion
from hisemble.graph import build_knn_graph_chunked_cpu
from hisemble.io import validate_embeddings

SECTIONS = ['E11_0-S1', 'E13_5-S1', 'E15_5-S1', 'E18_5-S1']


def spatial_coordinates(adata):
    for columns in [('array_col', 'array_row'), ('x', 'y')]:
        if all(name in adata.obs.columns for name in columns):
            return adata.obs[list(columns)].to_numpy(dtype=float)
    if 'spatial' in adata.obsm:
        coordinates = np.asarray(adata.obsm['spatial'], dtype=float)
        if coordinates.shape == (adata.n_obs, 2):
            return coordinates
    raise ValueError('Provide array_col/array_row, x/y, or obsm["spatial"] with two coordinates.')


def prepare_embeddings(adata, mode, seed=42, components=30, graphst_key='emb',
                       device='auto', graphst_seed=41, graphst_epochs=600):
    if not adata.obs_names.is_unique:
        raise ValueError('Cell identifiers must be unique.')
    coordinates = spatial_coordinates(adata)
    adata.obsm['spatial'] = coordinates
    if mode == 'ATAC':
        if 'X_lsi' not in adata.obsm:
            if components < 2 or components > adata.n_vars:
                raise ValueError('ATAC SVD components must be between 2 and the feature count.')
            tfidf = TfidfTransformer(norm='l2', use_idf=True).fit_transform(adata.X)
            adata.obsm['X_lsi'] = TruncatedSVD(
                n_components=components, random_state=seed
            ).fit_transform(tfidf)[:, 1:]
        return {'Spatial': coordinates, 'LSI': np.asarray(adata.obsm['X_lsi'])}
    if mode != 'RNA':
        raise ValueError('mode must be RNA or ATAC')
    if 'X_pca' not in adata.obsm:
        import scanpy as sc
        if components >= min(adata.shape):
            raise ValueError('RNA PCA components must be smaller than both matrix dimensions.')
        sc.pp.highly_variable_genes(
            adata, flavor='seurat_v3', n_top_genes=min(3000, adata.n_vars)
        )
        sc.tl.pca(adata, n_comps=components, random_state=seed)
    if graphst_key not in adata.obsm:
        import torch
        from GraphST.GraphST import GraphST
        if 'highly_variable' not in adata.var:
            raise ValueError('For GraphST training with cached PCA, supply the original highly_variable gene mask.')
        selected_device = ('cuda:0' if torch.cuda.is_available() else 'cpu') if device == 'auto' else device
        model = GraphST(
            adata.copy(), device=torch.device(selected_device),
            random_seed=graphst_seed, epochs=graphst_epochs,
        )
        learned = model.train()
        if not learned.obs_names.equals(adata.obs_names):
            raise ValueError('GraphST changed cell identifiers or their order.')
        adata.obsm[graphst_key] = np.asarray(learned.obsm['emb'])
    return {
        'PCA': np.asarray(adata.obsm['X_pca']),
        'GraphST': np.asarray(adata.obsm[graphst_key]),
        'Spatial': coordinates,
    }


def fuse_views(embeddings, cell_ids, k_neighbors=20, k_spatial=8,
               iterations=20, seed=42, chunk_size=2048):
    validate_embeddings(
        embeddings, cell_ids_by_view={name: cell_ids for name in embeddings}
    )
    n_cells = len(cell_ids)
    if not (1 <= k_neighbors < n_cells and 1 <= k_spatial < n_cells):
        raise ValueError('Both neighbor counts must be positive and smaller than the cell count.')
    if iterations < 0:
        raise ValueError('iterations must be non-negative')
    graphs = {}
    for name, embedding in embeddings.items():
        graphs[name] = build_knn_graph_chunked_cpu(
            embedding, k_neighbors=k_spatial if name == 'Spatial' else k_neighbors,
            metric='cosine' if name == 'LSI' else 'gaussian',
            chunk_size=chunk_size, random_seed=seed,
        )
    fused = hisemble_fusion(
        list(graphs.values()), k_neighbors=k_neighbors,
        n_iterations=iterations, random_seed=seed,
    )
    return graphs, fused


def cluster_graph(graph, target_count, seed):
    rows, columns = graph.nonzero()
    network = ig.Graph(n=graph.shape[0], edges=list(zip(rows, columns)), directed=False)
    network.es['weight'] = graph[rows, columns].A1
    return np.asarray(find_best_matrix(network, target_count, seed=seed)).astype(str)


def process_sample(path, mode, args):
    np.random.seed(args.seed)
    random.seed(args.seed)
    adata = ad.read_h5ad(path)
    if args.label_key not in adata.obs or adata.obs[args.label_key].isna().any():
        raise ValueError(f'{path.name}: supply complete annotations in obs[{args.label_key!r}]')
    labels = adata.obs[args.label_key].astype(str).to_numpy()
    embeddings = prepare_embeddings(
        adata, mode, seed=args.seed, components=args.components,
        graphst_key=args.graphst_key, device=args.device,
        graphst_seed=args.graphst_seed, graphst_epochs=args.graphst_epochs,
    )
    graphs, fused = fuse_views(
        embeddings, adata.obs_names.to_numpy(), k_neighbors=args.k_neighbors,
        k_spatial=args.k_spatial, iterations=args.iterations,
        seed=args.seed, chunk_size=args.chunk_size,
    )
    target_count = len(np.unique(labels))
    methods = {name: graph for name, graph in graphs.items() if name != 'Spatial'}
    methods['Hisemble'] = fused
    section = path.name.split('_expr_preprocessed')[0].split('_atac_preprocessed')[0]
    rows = []
    adata.obs['Ground_Truth'] = pd.Categorical(labels)
    for name, graph in methods.items():
        predicted = cluster_graph(graph, target_count, args.seed)
        adata.obs[name] = pd.Categorical(predicted)
        rows.append({
            'section': section, 'modality': mode, 'method': name,
            'n_cells': adata.n_obs, 'n_views': len(embeddings) if name == 'Hisemble' else 1,
            'target_clusters': target_count, 'n_clusters': len(np.unique(predicted)),
            'ARI': adjusted_rand_score(labels, predicted),
            'NMI': normalized_mutual_info_score(labels, predicted),
            'k_neighbors': args.k_neighbors, 'k_spatial': args.k_spatial,
            'iterations': args.iterations, 'seed': args.seed,
        })
    adata.obsp['spatial_connectivities'] = graphs['Spatial']
    adata.obsp['hisemble_connectivities'] = fused
    adata.obsp['connectivities'] = fused
    adata.uns['hisemble_spatial'] = {
        'mode': mode, 'views': list(embeddings),
        'metrics': ['cosine' if name == 'LSI' else 'gaussian' for name in embeddings],
        'k_neighbors': args.k_neighbors, 'k_spatial': args.k_spatial,
        'iterations': args.iterations, 'seed': args.seed,
        'graphst_seed': args.graphst_seed, 'graphst_epochs': args.graphst_epochs,
        'resolution_fallback': 'last_tested',
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = args.output_dir / f'{section}_{mode}'
    pd.DataFrame(rows).to_csv(str(stem) + '_metrics.csv', index=False)
    if not args.no_save_h5ad:
        adata.write_h5ad(str(stem) + '_hisemble.h5ad', compression='gzip')
    print(pd.DataFrame(rows)[['section', 'modality', 'method', 'ARI', 'NMI']].to_string(index=False))
    return rows


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, default=Path('__HISEMBLE_SPATIAL_ROOT__'))
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--mode', choices=['RNA', 'ATAC', 'both'], default='both')
    parser.add_argument('--sections', nargs='+', default=SECTIONS)
    parser.add_argument('--label-key', default='Annotation_for_Combined')
    parser.add_argument('--k-neighbors', type=int, default=20)
    parser.add_argument('--k-spatial', type=int, default=8)
    parser.add_argument('--iterations', type=int, default=20)
    parser.add_argument('--components', type=int, default=30)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--graphst-key', default='emb')
    parser.add_argument('--graphst-seed', type=int, default=41)
    parser.add_argument('--graphst-epochs', type=int, default=600)
    parser.add_argument('--device', default='auto')
    parser.add_argument('--chunk-size', type=int, default=2048)
    parser.add_argument('--no-save-h5ad', action='store_true')
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    modes = ['RNA', 'ATAC'] if args.mode == 'both' else [args.mode]
    tasks = [(args.data_root / f"{section}_{'expr' if mode == 'RNA' else 'atac'}_preprocessed.h5ad", mode)
             for section in args.sections for mode in modes]
    missing = [str(path) for path, _ in tasks if not path.is_file()]
    if missing:
        raise FileNotFoundError('Missing input sections: ' + ', '.join(missing))
    rows = []
    for path, mode in tasks:
        rows.extend(process_sample(path, mode, args))
        pd.DataFrame(rows).to_csv(args.output_dir / 'spatial_metrics.csv', index=False)


if __name__ == '__main__':
    main()
