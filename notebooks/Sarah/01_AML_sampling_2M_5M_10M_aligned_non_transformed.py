import pandas as pd
# data = pd.read_csv("../../FACS_Sampling/datasets/AML_samples_and_myeloid_Aligned_non_normalized.csv")
data = pd.read_csv("FACS_Sampling/datasets/AML_samples_and_myeloid_Aligned_non_normalized.csv")


features = ['CD33_aligned', 'HLA-DR_aligned', 'CD19 CD3_aligned',
            'CD11b_aligned', 'CD66b_aligned', 'CD16_aligned', 'CD163_aligned',
            'CD14_aligned', 'CD2_aligned', 'Siglec8_aligned', 'CD13_aligned',
            'CD45_aligned', 'CD141_aligned', 'CD15_aligned', 'CD123_aligned',
            'CD11c_aligned', 'CD117_aligned', 'CD45RA_aligned', 'CD34_aligned',
            'ITGB7_aligned', 'CD88 CD89_aligned', 'FcER1A_aligned', 'Population']

import numpy as np
seed = 235224
np.random.seed(seed)
rand_index2 = np.random.choice(data.index.values, 2000000, replace=False)
df1 = data.iloc[rand_index2].drop(columns=["Unnamed: 0"])
df1.to_csv("../../FACS_Sampling/datasets/AML_Aligned_non_normalized_2M.csv")

seed = 235723
np.random.seed(seed)
rand_index5 = np.random.choice(data.index.values, 5000000, replace=False)
df5 = data.iloc[rand_index5].drop(columns=["Unnamed: 0"])
df5.to_csv("../../FACS_Sampling/datasets/AML_Aligned_non_normalized_5M.csv")

seed = 1646274
np.random.seed(seed)
rand_index10 = np.random.choice(data.index.values, 10000000, replace=False)
df10 = data.iloc[rand_index10].drop(columns=["Unnamed: 0"])
df10.to_csv("../../FACS_Sampling/datasets/AML_Aligned_non_normalized_10M.csv")