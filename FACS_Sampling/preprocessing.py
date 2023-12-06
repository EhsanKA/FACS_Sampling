import numpy as np
import pandas as pd

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
    df.reset_index(drop=True, inplace=True)

    # Write the selected rows to the output file without writing the index
    df.to_csv(output_file, index=False)
