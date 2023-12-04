import pandas as pd
import os
import time

# Environment variable for file path
file_path_env = os.getenv('MY_FACS_DATA_PATH')
if file_path_env is None:
    raise ValueError("Environment variable 'MY_FACS_DATA_PATH' is not set")

input_file = os.path.join(file_path_env, 'AML_samples_and_myeloid_cells_only_SSC_FSC.csv')
if not os.path.exists(input_file):
    raise FileNotFoundError(f"The file {input_file} does not exist")

# Standardized feature names (snake_case)
features = ['cd33_aligned', 'hla_dr_aligned', 'cd19_cd3_aligned', 'cd11b_aligned', 
            'cd66b_aligned', 'cd16_aligned', 'cd163_aligned', 'cd14_aligned', 
            'cd2_aligned', 'siglec8_aligned', 'cd13_aligned', 'cd45_aligned', 
            'cd141_aligned', 'cd15_aligned', 'cd123_aligned', 'cd11c_aligned', 
            'cd117_aligned', 'cd45ra_aligned', 'cd34_aligned', 'itgb7_aligned', 
            'cd88_cd89_aligned', 'fcer1a_aligned', 'population']

output_file = os.path.join(file_path_env,'datasets', 
                           'AML_samples_and_myeloid_aligned_non_normalized.csv')

chunk_size = 50000  # Adjust based on your system's memory capacity
first_chunk = True

for chunk in pd.read_csv(input_file, chunksize=chunk_size, low_memory=False):
    start_time = time.time()  # Start timer

    # Standardize chunk column names
    chunk.columns = [col.lower().replace(' ', '_').replace('-', '_') for col in chunk.columns]

    # Select only the required features
    if not set(features).issubset(set(chunk.columns)):
        missing_features = set(features) - set(chunk.columns)
        raise KeyError(f"These features are missing in the data: {missing_features}")

    selected_data = chunk[features]

    # Write the chunk to the output file
    if first_chunk:
        selected_data.to_csv(output_file, mode='w', index=False)
        first_chunk = False
    else:
        selected_data.to_csv(output_file, mode='a', index=False, header=False)
    end_time = time.time()  # End timer
    print(f"Iteration time: {end_time - start_time} seconds")
