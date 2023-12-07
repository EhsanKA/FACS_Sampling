import os
import numpy as np
import pandas as pd
import scanpy as sc
from FACS_Sampling.methods.methods import bin_sample, sample_random #, dist_sampling
from FACS_Sampling.utils import create_adata


# Get the file path from environment variable
file_path_env = os.getenv('MY_FACS_DATA_PATH')


def creat_ref_adata(size_=2):

    # Construct the input file path
    input_file2 = os.path.join(file_path_env,'sara_data', f"sara_{size_}M.csv")

    # Read the data from the CSV file
    data = pd.read_csv(input_file2, low_memory=False)


    # Define the observation features
    obs_features = ['old_index', 'unique_id', 'sex', 'age', 'subtype', 'type', 'blastcount',
                    'survival_sorter', 'run','sample_id','alignment_mc_aligned',
                    'flowsom_cluster', 'flowsom_metacluster','population']

    # Create an AnnData object from the data
    adata = create_adata(data, obs_features=obs_features)

    # Apply arcsinh transformation to the data
    x = np.arcsinh(adata.X/500)
    adata.X = x

    # Construct the output file path
    output_adata_ref = os.path.join(file_path_env,'sara_data', f"adata_ref_sara_{2}M.h5ad")

    # Write the AnnData object to a file
    adata.write(output_adata_ref)

    return adata

ref_adata_2 = create_ref_adata(size_=2)
ref_adata_5 = create_ref_adata(size_=5)
ref_adata_10 = create_ref_adata(size_=10)



def generate_adata(adata, rep = 20, ref_size=2):
    seed = 12345
    np.random.seed(seed)

    for i in range(rep):
        new_seed = np.random.randint(100000)

        ps, _ = bin_sample(adata, n_bins=40, s_size=250, seed=new_seed)
        rs = sample_random(adata, s_size=ps.size, seed=new_seed)
        # dist_s = dist_sampling(adata, sample_per_matrix=30, rng=1000, seed=new_seed)
        output_file1 = os.path.join(file_path_env,'sara_data', 'reps',f"random_adata_{ref_size}_{rep}__{i}.h5ad")
        output_file2 = os.path.join(file_path_env,'sara_data', 'reps',f"fsbs_adata_{ref_size}_{rep}__{i}.h5ad")
        adata[rs].write(output_file1)
        adata[ps].write(output_file2)


generate_adata(ref_adata_2, rep=30, ref_size=2)
generate_adata(ref_adata_5, rep=30, ref_size=5)
generate_adata(ref_adata_10, rep=30, ref_size=10)