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

# Define a function to normalize the strings in the DataFrame
def process_dataframe(df):
    # Identify potential categorical columns (object dtype)
    categorical_cols = df.select_dtypes(include=['object']).columns

    # Convert potential categorical columns to categorical type
    for col in categorical_cols:
        df[col] = df[col].astype('category')

    # Standardize all string columns to snake_case format
    df.columns = [col.lower().replace(' ', '_').replace('-', '_') for col in df.columns]

    # Convert the values in the columns to snake_case
    for col in df.columns:
        if df[col].dtype.name == 'category':
            df[col] = df[col].apply(lambda x: str(x).lower().replace(' ', '_').replace('-', '_') if pd.notnull(x) else x)

    return df


# Define a function to sample the data and save it to a file
def sample_and_save(data, seed, num_samples, output_file):
    np.random.seed(seed)
    rand_index = np.random.choice(data.index.values, num_samples, replace=False)
    df = data.iloc[rand_index]

    # Reset the index and keep the old index as a column
    df.reset_index(inplace=True)
    df = df.copy()
    df.rename(columns={'index': 'old_index'}, inplace=True)

    # Set a new index
    df.reset_index(inplace=True)

    # Write the selected rows to the output file without writing the index
    df.to_csv(output_file, index=False)


# Process the DataFrame
data = process_dataframe(data)

# Use the function for 2M rows
sample_and_save(data, seed=235224, num_samples=2000000, output_file=output_file2)

# Repeat the process for 5M rows
sample_and_save(data, seed=235723, num_samples=5000000, output_file=output_file5)

# Repeat the process for 10M rows
sample_and_save(data, seed=1646274, num_samples=10000000, output_file=output_file10)

