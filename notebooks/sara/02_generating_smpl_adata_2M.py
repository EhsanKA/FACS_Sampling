import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import os
from FACS_Sampling.preprocessing.utils import create_adata

sc.set_figure_params(figsize=(8,8), fontsize=15, )


# Get the file path from the environment variable
file_path_env = os.getenv('MY_FACS_DATA_PATH')
input_file2 = os.path.join(file_path_env,'datasets', 'AML_data_2M.csv')

data = pd.read_csv(input_file2, low_memory=False)

data = data.drop(columns=['old_index'])
index = data.index
labels = data['population']
label_key = 'population'

obs_features = ['unique_id', 'sex', 'age', 'subtype', 'type', 'blastcount',
                'survival_sorter', 'run','sample_id','alignment_mc_aligned',
                'flowsom_cluster', 'flowsom_metacluster','population']
adata = create_adata(data, obs_features=obs_features)

print(adata.X.max())
x = np.arcsinh(adata.X/500)
adata.X = x

from FACS_Sampling.methods.methods import bin_sample, sample_random, dist_sampling
from FACS_Sampling.preprocessing.utils import set_reference_colors
from FACS_Sampling.plotting.plotting import change_cmap

cluster_key = label_key
ref_adata = adata



def bootstrap(ref_adata, label_key='labels', rep=5, seed=12345):
    np.random.seed(seed)

    groupped_df = pd.DataFrame(data=None, columns=["count", label_key, "method"])
    for i in range(rep):
        new_seed = np.random.randint(100000)

        ps, _ = bin_sample(ref_adata, n_bins=40, s_size=250, seed=new_seed)
        rs = sample_random(ref_adata, s_size=ps.size, seed=new_seed)
        # dist_s = dist_sampling(ref_adata, sample_per_matrix=30, rng=1000, seed=new_seed)
        
        def automation(samples, method_type='Random'):
            
            ref_adata[samples].write(f"../.bootstrap_adata/{method_type}_adata_{rep}__{i}.h5ad")
            a = ref_adata[samples].obs[label_key].value_counts().to_list()
            b = list(ref_adata[samples].obs[label_key].value_counts().keys())
            c = [method_type for _ in a]
            g_df = pd.DataFrame({"count": a, label_key: b, "method": c})
            return g_df
        
        
        g_df1 = automation(rs, method_type='Random')
        g_df2 = automation(ps, method_type='FSBS')
        
        groupped_df = pd.concat([groupped_df, g_df1, g_df2])

    return groupped_df


df_bootstrap = bootstrap(ref_adata, label_key=label_key, rep=20)