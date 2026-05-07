"""Python benchmark helpers for reticulate-driven spatial clustering wrappers."""

from __future__ import annotations

import importlib
import inspect
import random
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd


def _is_sequence(value):
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes))


def _as_sample_ids(sample_ids, n_items):
    if sample_ids is None:
        return [f"sample_{idx + 1}" for idx in range(n_items)]
    if isinstance(sample_ids, (str, bytes)):
        sample_ids = [sample_ids]
    sample_ids = [str(x) for x in list(sample_ids)]
    if len(sample_ids) != n_items:
        raise ValueError(
            f"Expected {n_items} sample IDs but received {len(sample_ids)}."
        )
    return sample_ids


def _expand_sample_arg(value, sample_ids, name):
    n_items = len(sample_ids)
    if isinstance(value, Mapping):
        missing = [sid for sid in sample_ids if sid not in value]
        if missing:
            raise ValueError(
                f"Argument '{name}' is missing entries for: {', '.join(missing)}."
            )
        return [value[sid] for sid in sample_ids]
    if _is_sequence(value):
        if len(value) != n_items:
            raise ValueError(
                f"Argument '{name}' must have length 1 or {n_items}, "
                f"not {len(value)}."
            )
        return list(value)
    return [value] * n_items


def _seed_all(seed):
    import torch

    random.seed(int(seed))
    np.random.seed(int(seed))
    torch.manual_seed(int(seed))
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(int(seed))


def _resolve_device(device):
    import torch

    if device is None:
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if str(device) == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(str(device))


def _load_image(image_source):
    if image_source is None:
        return None
    if isinstance(image_source, str):
        cv2 = importlib.import_module("cv2")
        image = cv2.imread(image_source)
        if image is None:
            raise ValueError(f"Failed to read image from '{image_source}'.")
        return image
    return np.asarray(image_source)


def _build_adata(counts, coords, barcodes, genes, sample_id):
    anndata = importlib.import_module("anndata")
    sparse = importlib.import_module("scipy.sparse")

    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] < 2:
        raise ValueError("Each coordinate matrix must have at least two columns.")

    n_cells = coords.shape[0]
    if sparse.issparse(counts):
        x_mat = counts.T.tocsr()
        n_genes, n_cells_counts = counts.shape
    else:
        counts = np.asarray(counts)
        if counts.ndim != 2:
            raise ValueError("Each count matrix must be two-dimensional.")
        n_genes, n_cells_counts = counts.shape
        x_mat = counts.T

    if n_cells_counts != n_cells:
        raise ValueError(
            "Counts and coordinates disagree on the number of cells/spots: "
            f"{n_cells_counts} vs {n_cells}."
        )

    barcodes = [str(x) for x in list(barcodes)]
    genes = [str(x) for x in list(genes)]
    if len(barcodes) != n_cells:
        raise ValueError(
            f"Expected {n_cells} barcodes for sample '{sample_id}' but received "
            f"{len(barcodes)}."
        )
    if len(genes) != n_genes:
        raise ValueError(
            f"Expected {n_genes} genes but received {len(genes)}."
        )

    obs = pd.DataFrame(index=barcodes)
    obs["sample_id"] = str(sample_id)
    obs["x_array"] = coords[:, 0]
    obs["y_array"] = coords[:, 1]
    obs["x_pixel"] = coords[:, 0]
    obs["y_pixel"] = coords[:, 1]

    var = pd.DataFrame(index=genes)
    adata = anndata.AnnData(X=x_mat, obs=obs, var=var)
    adata.var_names_make_unique()
    adata.obsm["spatial"] = coords[:, :2]
    return adata


def _ensure_sparse_array_property():
    sparse = importlib.import_module("scipy.sparse")
    for cls_name in ("csr_matrix", "csc_matrix", "coo_matrix"):
        cls = getattr(sparse, cls_name, None)
        if cls is not None and not hasattr(cls, "A"):
            cls.A = property(lambda self: self.toarray())


def _obs_values(adata, column):
    if column not in adata.obs:
        return None
    return [str(x) for x in adata.obs[column].astype(str).tolist()]


def _embedding_values(adata, key):
    if key not in adata.obsm:
        return None
    return np.asarray(adata.obsm[key])


def _as_nested_lists(array_like):
    if array_like is None:
        return None
    return np.asarray(array_like).tolist()


def _normalise_log1p(adata, target_sum=1e4, prefer_per_cell=False):
    sc = importlib.import_module("scanpy")
    if prefer_per_cell and hasattr(sc.pp, "normalize_per_cell"):
        sc.pp.normalize_per_cell(adata)
    else:
        sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    return sc


def _accepts_arg(func, arg_name):
    try:
        return arg_name in inspect.signature(func).parameters
    except (TypeError, ValueError):
        return False


def _init_graphst():
    graphst_root = importlib.import_module("GraphST")
    graphst_obj = getattr(graphst_root, "GraphST", None)
    if graphst_obj is not None and hasattr(graphst_obj, "GraphST"):
        graphst_ctor = graphst_obj.GraphST
    elif callable(graphst_obj):
        graphst_ctor = graphst_obj
    else:
        graphst_ctor = importlib.import_module("GraphST.GraphST").GraphST
    clustering_fn = importlib.import_module("GraphST.utils").clustering
    return graphst_ctor, clustering_fn


def run_graphst(
    count_matrices,
    coord_matrices,
    barcodes,
    genes,
    sample_ids=None,
    n_clusters=None,
    seed=567,
    device="cpu",
    cluster_method="mclust",
    refinement=True,
    radius=50,
    datatype="Slide",
    epochs=600,
    start=0.1,
    end=3.0,
    increment=0.01,
):
    sample_ids = _as_sample_ids(sample_ids, len(count_matrices))
    _ensure_sparse_array_property()
    n_clusters = _expand_sample_arg(n_clusters, sample_ids, "n_clusters")
    datatype = _expand_sample_arg(datatype, sample_ids, "datatype")

    graphst_ctor, clustering_fn = _init_graphst()
    torch_device = _resolve_device(device)

    results = {
        "method": "GraphST",
        "sample_ids": sample_ids,
        "barcodes": {},
        "clusters": {},
        "cluster_columns": {},
        "embeddings": {},
        "params": {
            "seed": int(seed),
            "cluster_method": str(cluster_method),
            "refinement": bool(refinement),
            "radius": int(radius),
            "epochs": int(epochs),
        },
    }

    for sid, counts, coords, cells, n_clust, data_type in zip(
        sample_ids, count_matrices, coord_matrices, barcodes, n_clusters, datatype
    ):
        _seed_all(seed)
        adata = _build_adata(counts, coords, cells, genes, sid)

        model = graphst_ctor(
            adata,
            device=torch_device,
            epochs=int(epochs),
            random_seed=int(seed),
            datatype=str(data_type),
        )
        adata = model.train()
        clustering_fn(
            adata,
            int(n_clust),
            radius=int(radius),
            method=str(cluster_method),
            start=float(start),
            end=float(end),
            increment=float(increment),
            refinement=bool(refinement),
        )

        results["barcodes"][sid] = [str(x) for x in adata.obs_names.tolist()]
        results["clusters"][sid] = _obs_values(adata, "domain")
        results["cluster_columns"][sid] = {
            "domain": _obs_values(adata, "domain"),
        }
        if "mclust" in adata.obs:
            results["cluster_columns"][sid]["mclust"] = _obs_values(adata, "mclust")
        if "leiden" in adata.obs:
            results["cluster_columns"][sid]["leiden"] = _obs_values(adata, "leiden")
        if "louvain" in adata.obs:
            results["cluster_columns"][sid]["louvain"] = _obs_values(adata, "louvain")
        results["embeddings"][sid] = _as_nested_lists(_embedding_values(adata, "emb"))

    return results


def run_spagcn(
    count_matrices,
    coord_matrices,
    barcodes,
    genes,
    sample_ids=None,
    n_clusters=None,
    seed=567,
    images=None,
    use_histology=False,
    refine=False,
    refine_shape="square",
    alpha=1.0,
    beta=49.0,
    p=0.5,
    l_value=None,
    res=None,
    init_spa=True,
    init="louvain",
    tol=5e-3,
    lr=0.05,
    max_epochs=200,
    search_res_start=0.7,
    search_res_step=0.1,
    search_res_tol=5e-3,
    search_res_epochs=20,
    l_search_start=0.01,
    l_search_end=1000,
    l_search_tol=0.01,
    l_search_max_run=100,
    min_cells=3,
):
    sc = importlib.import_module("scanpy")
    spg = importlib.import_module("SpaGCN")
    _ensure_sparse_array_property()

    sample_ids = _as_sample_ids(sample_ids, len(count_matrices))
    n_clusters = _expand_sample_arg(n_clusters, sample_ids, "n_clusters")
    images = _expand_sample_arg(images, sample_ids, "images") if images is not None else [None] * len(sample_ids)
    l_value = _expand_sample_arg(l_value, sample_ids, "l_value") if l_value is not None else [None] * len(sample_ids)
    res = _expand_sample_arg(res, sample_ids, "res") if res is not None else [None] * len(sample_ids)
    refine_shape = _expand_sample_arg(refine_shape, sample_ids, "refine_shape")

    results = {
        "method": "SpaGCN",
        "sample_ids": sample_ids,
        "barcodes": {},
        "clusters": {},
        "cluster_columns": {},
        "embeddings": {},
        "params": {
            "seed": int(seed),
            "use_histology": bool(use_histology),
            "refine": bool(refine),
            "alpha": float(alpha),
            "beta": float(beta),
            "p": float(p),
            "max_epochs": int(max_epochs),
        },
    }

    for sid, counts, coords, cells, n_clust, image, l_cur, res_cur, shape in zip(
        sample_ids,
        count_matrices,
        coord_matrices,
        barcodes,
        n_clusters,
        images,
        l_value,
        res,
        refine_shape,
    ):
        adata = _build_adata(counts, coords, cells, genes, sid)
        img = _load_image(image)
        histology_flag = bool(use_histology) and img is not None

        x_pixel = adata.obsm["spatial"][:, 0].tolist()
        y_pixel = adata.obsm["spatial"][:, 1].tolist()
        x_array = adata.obs["x_array"].tolist()
        y_array = adata.obs["y_array"].tolist()

        adj = spg.calculate_adj_matrix(
            x=x_pixel,
            y=y_pixel,
            x_pixel=x_pixel,
            y_pixel=y_pixel,
            image=img,
            beta=float(beta),
            alpha=float(alpha),
            histology=histology_flag,
        )

        spg.prefilter_genes(adata, min_cells=int(min_cells))
        spg.prefilter_specialgenes(adata)
        if hasattr(sc.pp, "normalize_per_cell"):
            sc.pp.normalize_per_cell(adata)
        else:
            sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        if l_cur is None:
            l_cur = spg.search_l(
                float(p),
                adj,
                start=float(l_search_start),
                end=float(l_search_end),
                tol=float(l_search_tol),
                max_run=int(l_search_max_run),
            )
        if res_cur is None:
            res_cur = spg.search_res(
                adata,
                adj,
                l_cur,
                int(n_clust),
                start=float(search_res_start),
                step=float(search_res_step),
                tol=float(search_res_tol),
                lr=float(lr),
                max_epochs=int(search_res_epochs),
                r_seed=int(seed),
                t_seed=int(seed),
                n_seed=int(seed),
            )

        _seed_all(seed)
        clf = spg.SpaGCN()
        clf.set_l(l_cur)
        clf.train(
            adata,
            adj,
            init_spa=bool(init_spa),
            init=str(init),
            n_clusters=int(n_clust),
            res=float(res_cur),
            tol=float(tol),
            lr=float(lr),
            max_epochs=int(max_epochs),
        )
        y_pred, _ = clf.predict()
        pred_labels = [str(x) for x in list(y_pred)]

        cluster_columns = {"pred": pred_labels}
        chosen_labels = pred_labels
        if refine:
            adj_2d = spg.calculate_adj_matrix(x=x_array, y=y_array, histology=False)
            refined = spg.refine(
                sample_id=adata.obs.index.tolist(),
                pred=pred_labels,
                dis=adj_2d,
                shape=str(shape),
            )
            chosen_labels = [str(x) for x in list(refined)]
            cluster_columns["refined_pred"] = chosen_labels

        results["barcodes"][sid] = [str(x) for x in adata.obs_names.tolist()]
        results["clusters"][sid] = chosen_labels
        results["cluster_columns"][sid] = cluster_columns
        results["embeddings"][sid] = None

    return results


def run_stagate(
    count_matrices,
    coord_matrices,
    barcodes,
    genes,
    sample_ids=None,
    n_clusters=None,
    seed=567,
    device="cpu",
    rad_cutoff=150,
    cluster_method="mclust",
    n_top_genes=3000,
    resolution=1.0,
    use_rep="STAGATE",
    target_sum=1e4,
    n_epochs=None,
):
    sc = importlib.import_module("scanpy")
    stagate = importlib.import_module("STAGATE_pyG")
    _ensure_sparse_array_property()

    sample_ids = _as_sample_ids(sample_ids, len(count_matrices))
    n_clusters = _expand_sample_arg(n_clusters, sample_ids, "n_clusters")
    rad_cutoff = _expand_sample_arg(rad_cutoff, sample_ids, "rad_cutoff")

    results = {
        "method": "STAGATE",
        "sample_ids": sample_ids,
        "barcodes": {},
        "clusters": {},
        "cluster_columns": {},
        "embeddings": {},
        "params": {
            "seed": int(seed),
            "cluster_method": str(cluster_method),
            "target_sum": float(target_sum),
            "use_rep": str(use_rep),
        },
    }

    train_fn = stagate.train_STAGATE
    torch_device = _resolve_device(device)

    for sid, counts, coords, cells, n_clust, cutoff in zip(
        sample_ids, count_matrices, coord_matrices, barcodes, n_clusters, rad_cutoff
    ):
        adata = _build_adata(counts, coords, cells, genes, sid)

        n_hvg = min(int(n_top_genes), int(adata.n_vars))
        if n_hvg > 0:
            sc.pp.highly_variable_genes(
                adata, flavor="seurat_v3", n_top_genes=n_hvg, subset=False
            )
        sc.pp.normalize_total(adata, target_sum=float(target_sum))
        sc.pp.log1p(adata)

        stagate.Cal_Spatial_Net(adata, rad_cutoff=float(cutoff))
        if hasattr(stagate, "Stats_Spatial_Net"):
            stagate.Stats_Spatial_Net(adata)

        _seed_all(seed)
        train_kwargs = {}
        if _accepts_arg(train_fn, "device"):
            train_kwargs["device"] = torch_device
        if n_epochs is not None:
            if _accepts_arg(train_fn, "n_epochs"):
                train_kwargs["n_epochs"] = int(n_epochs)
            elif _accepts_arg(train_fn, "num_epoch"):
                train_kwargs["num_epoch"] = int(n_epochs)
            elif _accepts_arg(train_fn, "epochs"):
                train_kwargs["epochs"] = int(n_epochs)
        if _accepts_arg(train_fn, "random_seed"):
            train_kwargs["random_seed"] = int(seed)

        adata = train_fn(adata, **train_kwargs)

        if str(cluster_method) == "mclust":
            adata = stagate.mclust_R(
                adata, used_obsm=str(use_rep), num_cluster=int(n_clust)
            )
            cluster_key = "mclust"
        else:
            sc.pp.neighbors(adata, use_rep=str(use_rep))
            if str(cluster_method) == "leiden":
                sc.tl.leiden(adata, resolution=float(resolution), random_state=int(seed))
                cluster_key = "leiden"
            elif str(cluster_method) == "louvain":
                sc.tl.louvain(
                    adata, resolution=float(resolution), random_state=int(seed)
                )
                cluster_key = "louvain"
            else:
                raise ValueError(
                    "cluster_method must be one of 'mclust', 'leiden', or 'louvain'."
                )

        results["barcodes"][sid] = [str(x) for x in adata.obs_names.tolist()]
        results["clusters"][sid] = _obs_values(adata, cluster_key)
        results["cluster_columns"][sid] = {cluster_key: _obs_values(adata, cluster_key)}
        results["embeddings"][sid] = _as_nested_lists(
            _embedding_values(adata, str(use_rep))
        )

    return results
