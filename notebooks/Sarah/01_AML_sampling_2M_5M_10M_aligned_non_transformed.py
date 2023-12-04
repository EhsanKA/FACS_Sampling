import pandas as pd
import numpy as np
import os

# Get the file path from the environment variable
file_path_env = os.getenv('MY_FACS_DATA_PATH')

# If the environment variable is not set, raise an error
if file_path_env is None:
    raise ValueError("Environment variable 'MY_FACS_DATA_PATH' is not set")

# Construct the input file path
input_file = os.path.join(file_path_env, 'datasets', 'AML_samples_and_myeloid_aligned_non_normalized.csv')

# If the input file does not exist, raise an error
if not os.path.exists(input_file):
    raise FileNotFoundError(f"The file {input_file} does not exist")

# Construct the output file paths
output_file2 = os.path.join(file_path_env,'datasets', 'AML_data_2M.csv')
output_file5 = os.path.join(file_path_env,'datasets', 'AML_data_5M.csv')
output_file10 = os.path.join(file_path_env,'datasets', 'AML_data_10M.csv')

# Read the input file into a pandas DataFrame
data = pd.read_csv(input_file, low_memory=False)

def sample_and_save(data, seed, num_samples, output_file):
    np.random.seed(seed)
    rand_index = np.random.choice(data.index.values, num_samples, replace=False)
    df = data.iloc[rand_index]

    # Reset the index and keep the old index as a column
    df.reset_index(inplace=True)
    df.rename(columns={'index': 'old_index'}, inplace=True)

    # Set a new index
    df.reset_index(inplace=True)

    # Write the selected rows to the output file without writing the index
    df.to_csv(output_file, index=False)
    
# Use the function for 2M rows
sample_and_save(data, seed=235224, num_samples=2000000, output_file=output_file2)

# Repeat the process for 5M rows
sample_and_save(data, seed=235723, num_samples=5000000, output_file=output_file5)

# Repeat the process for 10M rows
sample_and_save(data, seed=1646274, num_samples=10000000, output_file=output_file10)

