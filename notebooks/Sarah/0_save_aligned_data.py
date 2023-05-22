import pandas as pd
data = pd.read_csv("../../FACS_Sampling/datasets/AML_samples_and_myeloid_cells_only_SSC_FSC.csv")
features = ['CD33_aligned', 'HLA-DR_aligned', 'CD19 CD3_aligned',
            'CD11b_aligned', 'CD66b_aligned', 'CD16_aligned', 'CD163_aligned',
            'CD14_aligned', 'CD2_aligned', 'Siglec8_aligned', 'CD13_aligned',
            'CD45_aligned', 'CD141_aligned', 'CD15_aligned', 'CD123_aligned',
            'CD11c_aligned', 'CD117_aligned', 'CD45RA_aligned', 'CD34_aligned',
            'ITGB7_aligned', 'CD88 CD89_aligned', 'FcER1A_aligned', 'Population']
data[features].to_csv("../../FACS_Sampling/datasets/AML_samples_and_myeloid_Aligned_non_normalized.csv")