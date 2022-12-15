import pandas as pd
import scanpy as sc

sc.set_figure_params(figsize=(8, 8), fontsize=15, )


def read_tcell(address = "/fast/AG_Haghverdi/Ehsan_Karimiara/subarna_left/data_to_transfer/data_Focal/Final_Focal_Tcell_Selected_Patients_All_Cells.csv"):
    df1 = pd.read_csv(address)

    df_new = df1[['FJComp-APC-A_CCR7_asinh_aligned','FJComp-APC-H7-A_CD3_asinh_aligned','FJComp-APC-R700-A_C127_asinh_aligned',
                  'FJComp-BB515-A_CD45RA_asinh_aligned','FJComp-BB700-A_CD69_asinh_aligned','FJComp-BUV396-A_CD4_asinh_aligned',
                  'FJComp-BUV496-A_CD8_asinh_aligned', 'FJComp-BUV563-A_lineage_asinh_aligned','FJComp-BUV615-A_CD314_asinh_aligned',
                  'FJComp-BUV661-A_CXCR3_asinh_aligned','FJComp-BUV737-A_CD38_asinh_aligned', 'FJComp-BUV805-A_CD45RO_asinh_aligned',
                  'FJComp-BV421-A_CD28_asinh_aligned','FJComp-BV510-A_efluor506_asinh_aligned', 'FJComp-BV605-A_CD103_asinh_aligned',
                  'FJComp-BV650-A_PD1_asinh_aligned', 'FJComp-BV711-A_CD94_asinh_aligned', 'FJComp-BV750-A_TCRab_asinh_aligned',
                  'FJComp-BV786-A_CD95_asinh_aligned','FJComp-BYG584-A_CD39_asinh_aligned', 'FJComp-BYG670-A_ITBG7_asinh_aligned',
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
