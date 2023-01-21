import pandas as pd
import numpy as np
from FACS_Sampling.methods.methods import bin_sample, sample_random, proportional_sampling, corr_sampling


## bootstrap

def bootstrap(ref_adata, label_key='labels', rep=5, seed=12345):
    np.random.seed(seed)

    groupped_df = pd.DataFrame(data=None, columns=["count", "Refined_clustering", "batch"])
    for i in range(rep):
        new_seed = np.random.randint(100000)

        ps, _ = bin_sample(ref_adata, n_bins=14, s_size=20, seed=new_seed)
        rs = sample_random(ref_adata, s_size=ps.size, seed=new_seed)
        prop_s, _ = proportional_sampling(ref_adata, n_bins=14, total_size=150, power=0.01, seed=new_seed)

        def automation(samples, method_type='Uniform'):
            a = ref_adata[samples].obs[label_key].value_counts().to_list()
            b = list(ref_adata[samples].obs[label_key].value_counts().keys())
            c = [method_type for _ in a]
            g_df = pd.DataFrame({"count": a, label_key: b, "method": c})
            return g_df

        g_df1 = automation(rs, method_type='Uniform')
        g_df2 = automation(ps, method_type='FSBS')
        g_df3 = automation(prop_s, method_type='Proportional')

        groupped_df = pd.concat([groupped_df, g_df1, g_df2, g_df3])
    return groupped_df
## label transfer
