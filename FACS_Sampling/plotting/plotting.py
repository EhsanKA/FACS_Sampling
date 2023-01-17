import scanpy as sc
sc.set_figure_params(figsize=(8, 8), fontsize=15, )

def plot_umap(adata, dictionary, celltype_key ='Refined_clustering'):
    sc.pp.neighbors(adata, n_neighbors = 30)
    sc.tl.umap(adata)
    ref_cat = adata.obs[celltype_key].cat.categories
    dic = [dictionary[key] for key in ref_cat]
    adata.uns[celltype_key + '_colors'] = dic
    # sc.pl.umap(adata, color=celltype_key)


