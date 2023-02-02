import scanpy as sc
sc.set_figure_params(figsize=(8, 8), fontsize=15, )
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def change_cmap(adata, dictionary, celltype_key ='labels'):
    sc.pp.neighbors(adata, n_neighbors = 30)
    sc.tl.umap(adata)
    ref_cat = adata.obs[celltype_key].cat.categories
    dic = [dictionary[key] for key in ref_cat]
    adata.uns[celltype_key + '_colors'] = dic
    # sc.pl.umap(adata, color=celltype_key)


def kde_plot_genes(adata, n_bins=14, threshold=0.15):
    sc.set_figure_params(figsize=(50, 40), fontsize=15)
    for count, gene in enumerate(adata.var_names.values, start=1):
        plt.subplot(6,7,count)
        mini = adata[:, [gene]]
        std = np.std(np.array(mini.X))
        mean = np.mean(np.array(mini.X))

        measures = []
        var = (n_bins - 2) / 2
        for i in range(n_bins - 1):
            measures.append(mean - var * std)
            var -= 1

        sns.kdeplot(np.array(mini.X.reshape(-1)), fill=True,)

        for i in range(len(measures)):
            plt.axvline(x=measures[i], color='red', linestyle='--')

        plt.axhline(y=threshold, color='b', linestyle='-')

        plt.title("{}".format(gene))
    plt.savefig("kde_plots_genes.png")