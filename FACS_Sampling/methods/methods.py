import numpy as np
from scipy import stats
from numpy import linalg as LA



def sample_random(adata, s_size=1000, seed=12345):
    np.random.seed(seed)
    rand_index1 = np.random.choice(adata.obs.index, s_size, replace=False)

    return rand_index1


def sample_diagonal(adata, s_size=1000, seed=12345):
    np.random.seed(seed)

    x = adata.X
    prob = LA.norm(x, axis=1)
    s = sum(prob)
    prob = [float(i) / s for i in prob]
    rand_index2 = np.random.choice(adata.obs.index, s_size, p=prob, replace=False)

    return rand_index2


def bin_sample(adata, s_size=20, n_bins=10, seed=12345):
    np.random.seed(seed)
    total_indices = np.array([])
    for gene in adata.var_names.values:
        indices = np.array([])
        mini = adata[:, [gene]]
        std = np.std(np.array(mini.X))
        mean = np.mean(np.array(mini.X))

        measures = [-np.inf]
        var = (n_bins - 2) / 2
        for i in range(n_bins - 1):
            measures.append(mean - var * std)
            var -= 1
        measures.append(np.inf)

        for item in range(len(measures) - 1):
            measure1, measure2 = measures[item], measures[item + 1]
            local_indices = adata[np.array(mini.X > measure1) & np.array(mini.X < measure2)].obs.index.to_numpy()
            if len(local_indices) > s_size:
                indices = np.append(indices, np.random.choice(local_indices, s_size, replace=False))
                # if item==0 or item==len(measures)-2:
                #     indices = np.append(indices,local_indices)
                # else:
                #     indices = np.append(indices, np.random.choice(local_indices, s_size, replace=False))
            elif len(local_indices) > 0:
                indices = np.append(indices, local_indices)

        total_indices = np.append(total_indices, indices)

    return np.unique(total_indices), total_indices


def bin_sampleـpercentile(adata, total_size=100, n_bins=10, seed=12345):
    np.random.seed(seed)
    total_indices = np.array([])
    for gene in adata.var_names.values:
        indices = np.array([])
        mini = adata[:, [gene]]
        percent = 100 / n_bins

        min_v = mini.copy().X.min()
        max_v = mini.copy().X.max()
        measures = [np.percentile(mini.X, 0)]

        for i in range(n_bins):
            measures.append(np.percentile(mini.X, (i + 1) * percent))

        lenghts = [measures[i + 1] - measures[i] for i in range(len(measures) - 1)]

        prob = lenghts / sum(lenghts)
        prob *= total_size
        prob = [int(x) for x in prob]

        for item in range(len(measures) - 1):
            measure1, measure2 = measures[item], measures[item + 1]
            local_indices = adata[np.array(mini.X > measure1) & np.array(mini.X < measure2)].obs.index.to_numpy()
            if len(local_indices) < prob[item]:
                prob[item] = len(local_indices)
            indices = np.append(indices, np.random.choice(local_indices, prob[item], replace=False))

        total_indices = np.append(total_indices, indices)

    return np.unique(total_indices), total_indices


def bin_sampleـpercentile_power(adata, total_size=100, power=0.75, n_bins=10, seed=12345):
    np.random.seed(seed)
    total_indices = np.array([])
    for gene in adata.var_names.values:
        indices = np.array([])
        mini = adata[:, [gene]]
        percent = 100 / n_bins

        min_v = mini.copy().X.min()
        max_v = mini.copy().X.max()
        measures = [np.percentile(mini.X, 0)]

        for i in range(n_bins):
            measures.append(np.percentile(mini.X, (i + 1) * percent))

        lenghts = [measures[i + 1] - measures[i] for i in range(len(measures) - 1)]

        lenghts = np.power(lenghts, power)
        prob = lenghts / sum(lenghts)
        prob *= total_size
        prob = [int(x) for x in prob]

        for item in range(len(measures) - 1):
            measure1, measure2 = measures[item], measures[item + 1]
            local_indices = adata[np.array(mini.X > measure1) & np.array(mini.X < measure2)].obs.index.to_numpy()
            if len(local_indices) < prob[item]:
                prob[item] = len(local_indices)
            indices = np.append(indices, np.random.choice(local_indices, prob[item], replace=False))

            # TODO# This couple of lines have to be updated!
            # if len(local_indices) > s_size:
            #     indices = np.append(indices, np.random.choice(local_indices, prob[item], replace=False))
            # elif len(local_indices) > 0:
            #     indices = np.append(indices,local_indices)

        total_indices = np.append(total_indices, indices)

    return np.unique(total_indices), total_indices


def proportional_sampling(adata, total_size=100, n_bins=10, power=1.0, seed=12345):
    np.random.seed(seed)
    total_indices = np.array([])
    for gene in adata.var_names.values:
        indices = np.array([])
        mini = adata[:, [gene]]
        std = np.std(np.array(mini.X))
        mean = np.mean(np.array(mini.X))

        measures = [-np.inf]
        var= (n_bins-2)/2
        for i in range(n_bins-1):
            measures.append(mean-var*std)
            var-=1
        measures.append(np.inf)

        cell_counts = []
        for item in range(len(measures)-1):
            measure1, measure2 = measures[item], measures[item+1]
            local_indices = adata[np.array(mini.X > measure1) & np.array(mini.X < measure2)].obs.index.to_numpy()
            cell_counts.append(len(local_indices))

        cell_counts = np.power(cell_counts, power)
        prob = cell_counts / sum(cell_counts)
        prob *= total_size
        prob = [int(x) for x in prob]

        for item in range(len(measures)-1):
            measure1, measure2 = measures[item], measures[item+1]
            local_indices = adata[np.array(mini.X > measure1) & np.array(mini.X < measure2)].obs.index.to_numpy()
            if len(local_indices) < prob[item]:
                prob[item] = len(local_indices)
            indices = np.append(indices, np.random.choice(local_indices, prob[item], replace=False))

        total_indices = np.append(total_indices,indices)

    return np.unique(total_indices), total_indices

from scipy.spatial.distance import pdist

from scipy import stats
from scipy.spatial.distance import pdist, squareform


def dist_sampling(adata, sample_per_matrix=4, rng=1000, seed=12345):
    adata.obs = adata.obs.reset_index()
    np.random.seed(seed)
    all_indices = adata.obs.index.to_numpy().astype(int)
    np.random.shuffle(all_indices)

    start_i = 0
    if all_indices.shape[0] % rng != 0:
        end_of_range = int(all_indices.shape[0] / rng + 1)
    else:
        end_of_range = int(all_indices.shape[0] / rng)

    output = []
    for i in range(end_of_range):
        samples = all_indices[start_i:start_i + rng]
        start_i += rng

        x = adata.X[samples,]
        # x = stats.zscore(x, axis=1)
        dists = squareform(pdist(x, metric='euclidean'))
        th = np.max(dists)
        mask1 = dists > th
        while (mask1.sum(axis=1) > 0).sum() < 5 * sample_per_matrix:
            th = th / 2
            mask1 = dists > th

        mask = mask1

        fin_scores = np.multiply(dists, mask).sum(axis=0)
        t = samples[np.argsort(fin_scores)[0:sample_per_matrix]]
        output.extend(t)

    output = adata.obs.loc[output]['index'].values
    adata.obs.index = adata.obs['index']
    del adata.obs['index']
    return np.array(output)
