import pandas as pd
import numpy as np
from FACS_Sampling.methods.methods import bin_sample, sample_random, proportional_sampling, dist_sampling
from FACS_Sampling.preprocessing.utils import find_mutual_nn, clean
from sklearn import neighbors
from sklearn.metrics import classification_report


## bootstrap

import numpy as np
import pandas as pd

def bootstrap(ref_adata, label_key='labels', rep=5, seed=12345):
    # Set the seed for reproducibility
    np.random.seed(seed)

    # Initialize an empty DataFrame to store results
    groupped_df = pd.DataFrame(data=None, columns=["count", label_key, "method"])

    # Repeat the process for a specified number of times
    for i in range(rep):
        # Generate a new seed for each iteration
        new_seed = np.random.randint(100000)

        # Generate samples using different methods
        ps, _ = bin_sample(ref_adata, n_bins=14, s_size=20, seed=new_seed)
        rs = sample_random(ref_adata, s_size=ps.size, seed=new_seed)
        prop_s, _ = dist_sampling(ref_adata, n_bins=14, total_size=150, power=0.01, seed=new_seed)

        # Function to automate the process of creating a DataFrame from samples
        def automation(samples, method_type='Random'):
            # Get the count of each label
            a = ref_adata[samples].obs[label_key].value_counts().to_list()
            # Get the corresponding labels
            b = list(ref_adata[samples].obs[label_key].value_counts().keys())
            # Create a list with the method type
            c = [method_type for _ in a]
            # Create a DataFrame
            g_df = pd.DataFrame({"count": a, label_key: b, "method": c})
            return g_df

        # Create DataFrames for each method
        g_df1 = automation(rs, method_type='Random')
        g_df2 = automation(ps, method_type='FSBS')
        g_df3 = automation(prop_s, method_type='Proportional')

        # Concatenate all DataFrames
        groupped_df = pd.concat([groupped_df, g_df1, g_df2, g_df3])

    # Return the final DataFrame
    return groupped_df


## label transfer

def get_knn_classification_report(index_subset, ref_adata, label_key='labels'):
    # Subset the data and split into training and testing sets
    small_adata = ref_adata[index_subset]
    X_train, y_train = small_adata.X, small_adata.obs[label_key]
    
    mask = ref_adata.obs.index.isin(small_adata.obs.index)
    X_test, y_test = ref_adata[~mask].X, ref_adata[~mask].obs[label_key]

    # Initialize, train and predict using a KNN classifier
    clf = neighbors.KNeighborsClassifier(n_neighbors=5).fit(X_train, y_train)
    y_pred = clf.predict(X_test)

    # Generate a classification report and convert it to a DataFrame
    cr = classification_report(y_test, y_pred, output_dict=True)
    df_cr = pd.DataFrame(cr)
    
    # Return the DataFrame and the raw report
    return df_cr, cr


def get_mnn_classification_report(index_subset, ref_adata, label_key='labels'):
    small_adata = ref_adata[index_subset]
    bad_inds = small_adata.obs.index
    bad_df = ref_adata.obs.index.isin(bad_inds)

    X_train = small_adata.X
    y_train = small_adata.obs[label_key]
    X_test = ref_adata[~bad_df].X
    y_test = ref_adata[~bad_df].obs[label_key]

    y_pred = MNN_classifier(X_train, X_test, y_train, k=100)

    cr = classification_report(y_test, y_pred, output_dict=True)
    crc = classification_report(y_test, y_pred)
    df_cr = pd.DataFrame(cr)
    return df_cr, crc


def MNN_classifier(train_data, test_data, train_labels, k=100):
    mutual_1, mutual_2 = find_mutual_nn(train_data, test_data, k1=k, k2=k)
    dictionary = clean(mutual_1, mutual_2, train_labels)
    y_pred = []
    not_found = []

    clf = neighbors.KNeighborsClassifier(n_neighbors=5)
    clf.fit(train_data, train_labels)

    for i in range(test_data.shape[0]):
        if i in dictionary.keys():
            label = max(set(dictionary[i]), key=dictionary[i].count)
        else:
            label = clf.predict(test_data[i].reshape(1, -1))[0]
            not_found.append(i)
        y_pred.append(label)
    print("{} cells could not get MNN's in the other datasets.".format(len(not_found)))

    return y_pred


def repeat_classifier(ref_adata, label_key='labels', classifier=get_knn_classification_report, rep=5, seed=12345):
    np.random.seed(seed)
    output_random = []
    output_fsbs = []
    output_pbs = []

    for i in range(rep):
        new_seed = np.random.randint(10675)

        ps, _ = bin_sample(ref_adata, n_bins=14, s_size=20, seed=new_seed)
        rs = sample_random(ref_adata, s_size=ps.size, seed=new_seed)
        prop_s, _ = proportional_sampling(ref_adata, n_bins=14, total_size=150, power=0.01, seed=new_seed)

        output_random.append(classifier(rs, ref_adata, label_key=label_key)[0])
        output_fsbs.append(classifier(ps, ref_adata, label_key=label_key)[0])
        output_pbs.append(classifier(prop_s, ref_adata, label_key=label_key)[0])
    return output_random, output_fsbs, output_pbs
