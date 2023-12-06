import pandas as pd
import os
import time
import numpy as np
from FACS_Sampling.preprocessing import process_dataframe, sample_and_save
import scanpy as sc



# Environment variable for file path
file_path_env = os.getenv('MY_FACS_DATA_PATH')
if file_path_env is None:
    raise ValueError("Environment variable 'MY_FACS_DATA_PATH' is not set")

input_file = os.path.join(file_path_env, 'sara_data', 'sara_raw.csv')
if not os.path.exists(input_file):
    raise FileNotFoundError(f"The file {input_file} does not exist")

# Standardized feature names (snake_case)
features = ['cd33_aligned', 'hla_dr_aligned', 'cd19_cd3_aligned', 'cd11b_aligned', 
            'cd66b_aligned', 'cd16_aligned', 'cd163_aligned', 'cd14_aligned', 
            'cd2_aligned', 'siglec8_aligned', 'cd13_aligned', 'cd45_aligned', 
            'cd141_aligned', 'cd15_aligned', 'cd123_aligned', 'cd11c_aligned', 
            'cd117_aligned', 'cd45ra_aligned', 'cd34_aligned', 'itgb7_aligned', 
            'cd88_cd89_aligned', 'fcer1a_aligned', 'unique_id', 'sex', 'age',
            'subtype', 'type', 'blastcount', 'survival_sorter', 'run','sample_id',
            'alignment_mc_aligned', 'flowsom_cluster', 'flowsom_metacluster','population']

output_file = os.path.join(file_path_env,'sara_data', 'sara_data_processed.csv')
output_file2 = os.path.join(file_path_env,'sara_data', 'sara_2M.csv')
output_file5 = os.path.join(file_path_env,'sara_data', 'sara_5M.csv')
output_file10 = os.path.join(file_path_env,'sara_data', 'sara_10M.csv')


start_time = time.time()  # Start timer
data = pd.read_csv(input_file, low_memory=False)

# Standardize data column names
data.columns = [col.lower().replace(' ', '_').replace('-', '_') for col in data.columns]

# Select only the required features
if not set(features).issubset(set(data.columns)):
    missing_features = set(features) - set(data.columns)
    raise KeyError(f"These features are missing in the data: {missing_features}")

selected_data = data[features].copy()  # Copy the selected features to a new DataFrame
del data

# Normalize the strings in the DataFrame
selected_data = process_dataframe(selected_data)

# Write the data to the output file
selected_data.to_csv(output_file, index=False)

np.random.seed(42)  # Set the seed for reproducibility
new_seed1 = np.random.randint(100000)
new_seed2 = np.random.randint(100000)
new_seed3 = np.random.randint(100000)

# Use the function for 2M rows
sample_and_save(selected_data, seed=new_seed1, num_samples=2000000, output_file=output_file2)

# Repeat the process for 5M rows
sample_and_save(selected_data, seed=new_seed2, num_samples=5000000, output_file=output_file5)

# Repeat the process for 10M rows
sample_and_save(selected_data, seed=new_seed3, num_samples=10000000, output_file=output_file10)


end_time = time.time()  # End timer
print(f"Iteration time: {end_time - start_time} seconds")
