import scanpy as sc
sc.set_figure_params(figsize=(8, 8), fontsize=15, )

def setting_reference_colors(adata, cluster_key = 'Refined_clustering'):

    adatas = [adata[adata.obs[cluster_key].isin([clust])] for clust in adata.obs[cluster_key].cat.categories]

    for dat in adatas:
        if dat.n_obs > 100000:
            sc.pp.subsample(dat, n_obs=500)
        elif dat.n_obs > 10000:
            sc.pp.subsample(dat, n_obs=200)
        else:
            sc.pp.subsample(dat, n_obs=50)

    adata_downsampled = adatas[0].concatenate(*adatas[1:])

    sc.pp.neighbors(adata_downsampled, n_neighbors=30)
    sc.tl.umap(adata_downsampled)
    sc.pl.umap(adata_downsampled, color=cluster_key)

    down_cat = adata_downsampled.obs[cluster_key].cat.categories
    down_colors = adata_downsampled.uns[cluster_key+'_colors']
    dictionary = dict(zip(down_cat, down_colors))

    return dictionary


def plot_umap(ref_adata, dictionary, celltype_key = 'Refined_clustering'):
    sc.pp.neighbors(ref_adata, n_neighbors = 30)
    sc.tl.umap(ref_adata)
    ref_cat = ref_adata.obs[celltype_key].cat.categories
    dic = [dictionary[key] for key in ref_cat]
    ref_adata.uns[celltype_key+'_colors'] = dic
    sc.pl.umap(ref_adata, color=celltype_key)


