import pandas as pd
import scanpy as sc

sc.set_figure_params(figsize=(8, 8), fontsize=15, )

### ArchR color palettes https://rdrr.io/github/GreenleafLab/ArchR/src/R/ColorPalettes.R
# 20-colors
color_paletts = {
    20: ["#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC", "#90D5E4",
         "#89C75F", "#F37B7D", "#9983BD", "#D24B27", "#3BBCA8", "#6E4B9E", "#0C727C","#7E1416", "#D8A767", "#3D3D3D"],
    15:['#371377', '#7700FF','#9E0142', '#FF0080', '#DC494C', "#F88D51", "#FAD510", "#FFFF5F", '#88CFA4', '#238B45',
        "#02401B", "#0AD7D3", "#046C9A", "#A2A475", 'grey35'],}

# todo clean this part
# # 16-colors
# bear = c("1" = "#faa818", "2" = "#41a30d", "3" = "#fbdf72", "4" = "#367d7d", "5" = "#d33502", "6" = "#6ebcbc", "7" = "#37526d", "8" = "#916848", "9" = "#f5b390", "10" = "#342739", "11" = "#bed678", "12" = "#a6d9ee", "13" = "#0d74b6", "14" = "#60824f", "15" = "#725ca5", "16" = "#e0598b"),
# 12-colors
# paired = c("9" = "#A6CDE2", "1" = "#1E78B4", "3" = "#74C476", "12" = "#34A047", "11" = "#F59899", "2" = "#E11E26", "10" = "#FCBF6E", "4" = "#F47E1F", "5" = "#CAB2D6", "8" = "#6A3E98", "6" = "#FAF39B", "7" = "#B15928"),
# # 11-colors
# grove = c("11" = "#1a1334", "9" = "#01545a", "1" = "#017351", "6" = "#03c383", "8" = "#aad962", "2" = "#fbbf45", "10" = "#ef6a32", "3" = "#ed0345", "7" = "#a12a5e", "5" = "#710162", "4" = "#3B9AB2"),
# # 7-colors
# summerNight = c("1" = "#2a7185", "2" = "#a64027", "3" = "#fbdf72", "4" = "#60824f", "5" = "#9cdff0", "6" = "#022336", "7" = "#725ca5")


def read_tcell(
        address="/fast/AG_Haghverdi/Ehsan_Karimiara/subarna_left/data_to_transfer/data_Focal/Final_Focal_Tcell_Selected_Patients_All_Cells.csv"):
    df1 = pd.read_csv(address)

    df_new = df1[
        ['FJComp-APC-A_CCR7_asinh_aligned', 'FJComp-APC-H7-A_CD3_asinh_aligned', 'FJComp-APC-R700-A_C127_asinh_aligned',
         'FJComp-BB515-A_CD45RA_asinh_aligned', 'FJComp-BB700-A_CD69_asinh_aligned',
         'FJComp-BUV396-A_CD4_asinh_aligned',
         'FJComp-BUV496-A_CD8_asinh_aligned', 'FJComp-BUV563-A_lineage_asinh_aligned',
         'FJComp-BUV615-A_CD314_asinh_aligned',
         'FJComp-BUV661-A_CXCR3_asinh_aligned', 'FJComp-BUV737-A_CD38_asinh_aligned',
         'FJComp-BUV805-A_CD45RO_asinh_aligned',
         'FJComp-BV421-A_CD28_asinh_aligned', 'FJComp-BV510-A_efluor506_asinh_aligned',
         'FJComp-BV605-A_CD103_asinh_aligned',
         'FJComp-BV650-A_PD1_asinh_aligned', 'FJComp-BV711-A_CD94_asinh_aligned', 'FJComp-BV750-A_TCRab_asinh_aligned',
         'FJComp-BV786-A_CD95_asinh_aligned', 'FJComp-BYG584-A_CD39_asinh_aligned',
         'FJComp-BYG670-A_ITBG7_asinh_aligned',
         'FJComp-BYG790-A_TCRgd_asinh_aligned', 'FJComp-PE-CF594-A_CD25_asinh_aligned', 'Refined_clustering']]

    del df1
    df_new = df_new.replace(regex=[' T cells', 'memory', 'tissue-resident'], value=['', 'mem', 'tissue-res'])
    adata = create_adata(df_new, obs_features=['Refined_clustering'])

    return adata


def create_adata(df, obs_features=None):
    if obs_features is None:
        obs_features = ["Refined_clustering"]
    obs = pd.DataFrame()
    obs[obs_features] = df[obs_features]
    obs = obs.astype('category')
    var_names = df.drop(columns=obs_features).columns.values
    var = pd.DataFrame(index=var_names)
    X = df.drop(columns=obs_features).values
    adata = sc.AnnData(X, obs=obs, var=var)

    return adata


def set_reference_colors(adata, cluster_key='labels', make_plot=True, palette='ArchR'):
    adata.obs[cluster_key] = adata.obs[cluster_key].astype('category')
    adatas = [adata[adata.obs[cluster_key].isin([clust])]
              for clust in adata.obs[cluster_key].cat.categories]
    for dat in adatas:
        if dat.n_obs > 100000:
            sc.pp.subsample(dat, n_obs=500)
        elif dat.n_obs > 10000:
            sc.pp.subsample(dat, n_obs=200)
        else:
            sc.pp.subsample(dat, n_obs=50)

    adata_stratified = adatas[0].concatenate(*adatas[1:])

    sc.pp.neighbors(adata_stratified, n_neighbors=30)
    sc.tl.umap(adata_stratified)

    num_categories = adata.uns["labels_colors"].__len__()
    if palette is 'ArchR':
        if num_categories in color_paletts.keys():
            colors = color_paletts[num_categories]
    else:
        colors = adata_stratified.uns[cluster_key + '_colors']

    categories = adata_stratified.obs[cluster_key].cat.categories
    dictionary = dict(zip(categories, colors))

    if make_plot:
        sc.pl.umap(adata_stratified, color=cluster_key)

    return dictionary
