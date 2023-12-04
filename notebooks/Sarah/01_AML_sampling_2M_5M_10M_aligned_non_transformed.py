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

# Set the random seed for reproducibility
seed = 235224
np.random.seed(seed)

# Randomly select 2M rows from the DataFrame
rand_index2 = np.random.choice(data.index.values, 2000000, replace=False)
df2 = data.iloc[rand_index2]

# Write the selected rows to the output file
df2.to_csv(output_file2)

# Repeat the process for 5M rows
seed = 235723
np.random.seed(seed)
rand_index5 = np.random.choice(data.index.values, 5000000, replace=False)
df5 = data.iloc[rand_index5]
df5.to_csv(output_file5)

# Repeat the process for 10M rows
seed = 1646274
np.random.seed(seed)
rand_index10 = np.random.choice(data.index.values, 10000000, replace=False)
df10 = data.iloc[rand_index10]
df10.to_csv(output_file10)