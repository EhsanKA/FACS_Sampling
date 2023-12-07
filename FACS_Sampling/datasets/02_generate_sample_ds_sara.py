import numpy as np
import pandas as pd
import scanpy as sc
import os
from FACS_Sampeling.methods.method import bin_sample, sample_random



def generate_adata(rep = 20):
    seed = 12345
    np.random.seed(seed)
    file_path_env = os.getenv('MY_FACS_DATA_PATH')
    input_file = os.path.join(file_path_env,'sara_data', 'adata_ref_sara_2M.h5ad')
    ref_adata = sc.read(input_file)

    for i in range(rep):
        new_seed = np.random.randint(100000)

        ps, _ = bin_sample(ref_adata, n_bins=40, s_size=250, seed=new_seed)
        rs = sample_random(ref_adata, s_size=ps.size, seed=new_seed)
        # dist_s = dist_sampling(ref_adata, sample_per_matrix=30, rng=1000, seed=new_seed)
        output_file1 = os.path.join(file_path_env,'sara_data', 'reps',f"random_adata_{rep}__{i}.h5ad")
        output_file2 = os.path.join(file_path_env,'sara_data', 'reps',f"fsbs_adata_{rep}__{i}.h5ad")
        ref_adata[rs].write(output_file1)
        ref_adata[ps].write(output_file2)

generate_adata(rep=30)