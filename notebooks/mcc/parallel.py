import os
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import time
import scipy.sparse


# Get the file path from environment variable
file_path_env = '/fast/AG_Ohler/ekarimi/projects/FACS_Sampling/data'


REFERENCES = [5, 10, 20, 25, 30]
METHODS = ['random', 'cubic', 'hopper', 'atomic']
methods = ['random', 'cubic', 'hopper']
SIZES = [50000, 100000, 200000, 300000]
REPS = [i for i in range(5)]
label_key = 'celltype'


directory = "mcc/benchmark"
PATH = os.path.join(file_path_env, directory)



import numpy as np
import pandas as pd

def pca_bin_sample_(df, feature_importances, seed=12345):
    np.random.seed(seed)

    # Ensure num_pcs does not exceed the number of columns in df
    num_pcs = min(feature_importances.shape[0], df.shape[1])

    # Function to create bins and digitize
    def create_bins_and_digitize(data, n_bins):
        edges = np.linspace(data.min(), data.max(), n_bins + 1)
        bins = np.digitize(data, edges)
        return bins

    def compute_sample_bins(df, bin_sizes):
        bins = [create_bins_and_digitize(df.iloc[:, i], bin_sizes[i]) for i in range(num_pcs)]

        # Combine bins to form grid cells
        df['grid_cell'] = list(zip(*bins))
        
        return 

    compute_sample_bins(df, feature_importances)
    return

    
def set_min_to_two(pca):
    if len(pca.explained_variance_ratio_)< 40:
        out = np.ceil(pca.explained_variance_ratio_*100).astype(int)

    else:
        mul = 2.0/ pca.explained_variance_ratio_[19]
        out = np.ceil(pca.explained_variance_ratio_*mul).astype(int)
        
    return out[out>2]


import random

def find_threshold_index(sorted_grid_cells, threshold):
    cumulative = 0
    for index, frequency in sorted_grid_cells.value_counts().sort_index().items():
        cumulative += index * frequency
        if cumulative >= threshold:
            return index
    return None

def accumulate_indices_until_threshold(df, threshold, seed=1234):
    random.seed(seed)
    # Count the occurrences of each grid_cell
    grid_cell_counts = df['grid_cell'].value_counts()

    # Sort the grid_cells by count in ascending order
    sorted_grid_cells = grid_cell_counts.sort_values()

    # Find the threshold index
    threshold_index = find_threshold_index(sorted_grid_cells, threshold)
    print(f'threshold_index is : {threshold_index}')
    
    # Group the DataFrame by 'grid_cell'
    grouped_df = df.groupby('grid_cell')

    accumulated_indices = []
    accumulated_count = 0
    all_remainings_indices = []

    # Iterate over sorted grid_cells and accumulate indices
    for grid_cell in sorted_grid_cells.index:
        group_indices = grouped_df.get_group(grid_cell).index.tolist()
        if len(group_indices) < threshold_index:
            accumulated_indices.extend(group_indices)
            accumulated_count += len(group_indices)
        elif len(group_indices) == threshold_index:
            all_remainings_indices.extend(group_indices)
        else:
            break


    # Calculate how many more indices we need to reach the threshold
    remaining_count = threshold - accumulated_count
    print(f'remaining is : {remaining_count}')

    # Randomly select the remaining indices from the current group
    accumulated_indices.extend(random.sample(all_remainings_indices, remaining_count))
    
    return accumulated_indices

def generate_cubic(adata, size, seed = 1234):
    
    scaler = StandardScaler()
    data_standardized = scaler.fit_transform(adata.X)
    X = data_standardized


    random.seed(seed)
    print(f'********* #Start# *********')
    start_time = time.time()

    N_components= min(adata.shape[1], 100)

    pca = PCA(n_components=N_components)
    pca.fit(X)
    X_pca = pca.transform(X)
    df = pd.DataFrame(X_pca[:, :N_components], columns=[f'PC{i+1}' for i in range(N_components)])


    feature_importances = set_min_to_two(pca)

    pca_bin_sample_(df, feature_importances)

    threshold = size  # Set your desired threshold
    samples = accumulate_indices_until_threshold(df, threshold, seed=seed)

    elapsed_time = time.time() - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

        
    return samples, elapsed_time



def generate_geo(adata, size, seed=1234):
    
    scaler = StandardScaler()
    data_standardized = scaler.fit_transform(adata.X)

    print(f'********* #Start# *********')
    np.random.seed(seed)
    start_time = time.time()

    # Compute PCs.
    from fbpca import pca
    k = adata.shape[1]
    U, s, Vt = pca(data_standardized, k=k) # E.g., 4 PCs.
    X_dimred = U[:, :k] * s[:k]
    # Now, you are ready to sketch!

    # Sketch.
    from geosketch import gs
    N = size # Number of samples to obtain from the data set.
    samples = gs(X_dimred, N, replace=False)

    # X_sketch = X_dimred[sketch_index]
    elapsed_time = time.time() - start_time
    print(f"Elapsed time: {elapsed_time} seconds")
        
    return samples, elapsed_time

import sys
import io
from contextlib import redirect_stdout

def generate_hopper(adata, size, seed=1234):
    
    scaler = StandardScaler()
    data_standardized = scaler.fit_transform(adata.X)

    from hopper.treehopper.hoppers import hopper, treehopper, PCATreePartition
    print(f'********* #Start# *********')

    np.random.seed(seed)
    start_time = time.time()

    scaler = StandardScaler()
    data_standardized = scaler.fit_transform(adata.X)
    #data_standardized = scaler.fit_transform(adata_shuffled.X)
    X = data_standardized

    N = size # Number of samples to obtain from the data set.
    with io.StringIO() as buf, redirect_stdout(buf):
        th = treehopper(data_standardized, partition=PCATreePartition, max_partition_size=1000)
        sketch = th.hop(size)

    samples = th.path[:size]

    # X_sketch = X_dimred[sketch_index]
    elapsed_time = time.time() - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

        
    return samples, elapsed_time


def generate_random(adata, size, seed=1234):

    print(f'********* #Start# *********')
    np.random.seed(seed)
    start_time = time.time()


    samples = np.random.randint(0, adata.shape[0], size=size)

    elapsed_time = time.time() - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

        
    return samples, elapsed_time




import argparse
import numpy as np
import pickle

# Assuming the sampling functions are defined elsewhere and imported
# from sampling_functions import generate_cubic, generate_geo, generate_random, generate_hopper

def main(ref, method, size, rep, seed):
    # Define a dictionary to store the results
    results = []
    
    address = os.path.join(PATH, f"{ref}/adata.h5ad")
    adata = sc.read_h5ad(address)
    adata.obs[label_key] = adata.obs[label_key].astype('category')
    adata.var.index = adata.var.index.astype('object')
    
    if isinstance(adata.X, scipy.sparse.csr_matrix):
        adata.X = adata.X.toarray()

    method_dict = {
    "cubic": (generate_cubic, {"adata": adata, "size": size, "seed": seed}),
    "hopper": (generate_hopper, {"adata": adata, "size": size, "seed": seed}),
    "random": (generate_random, {"adata": adata, "size": size, "seed": seed}),
    }
    
    if method in method_dict:
        func, args = method_dict[method]
        results = func(**args)
    else:
        print(f"No function associated with {method}")
        return
        
    output_address = os.path.join(PATH, f"{ref}/{method}/{size}/{rep}/results.pkl")
    
    with open(output_address, 'wb') as handle:
        pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run sampling methods in parallel for a given Reference, Method, Size, Replicate, and Seed.")
    parser.add_argument("--ref", type=int, required=True, help="Reference to process")
    parser.add_argument("--method", type=str, required=True, help="Method to process")
    parser.add_argument("--size", type=int, required=True, help="Size to process")
    parser.add_argument("--rep", type=int, required=True, help="Replicate to process")
    parser.add_argument("--seed", type=int, required=True, help="Seed to process")

    args = parser.parse_args()
    print("###############################")
    print("************ New run **********")
    
    # ref = 1
    # method = 'cubic'
    # size = 50000
    # rep = 0
    # seed = 6547
    
    # main(args.batch_id, args.size)
    # main(ref, method, size, rep, seed)
    main(args.ref, args.method, args.size, args.rep, args.seed)
