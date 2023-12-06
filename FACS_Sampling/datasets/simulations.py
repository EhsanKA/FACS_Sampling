from scipy.integrate import quad
from scipy.stats import gaussian_kde
import numpy as np
import scanpy as sc

def find_bounds_for_a_gene_distribution(gene_adata, freq_threshold=0.15):
    # Assign data
    data = np.array(gene_adata.X.reshape(-1))
    gene_name = gene_adata.var_names[0]

    # Create the KDE plot
    kde = gaussian_kde(data)
    x_eval = np.linspace(np.min(data), np.max(data), num=10000)
    y_eval = kde(x_eval)

    # Define the horizontal line
    y_intercept = freq_threshold

    # Find the x-coordinates where the KDE plot intersects the horizontal line
    intersections = [x for x, y in zip(x_eval, y_eval) if y >= y_intercept]

    # from math import log10, floor
    # def round_to_1(x):
    #     return round(x, -int(floor(log10(abs(x)))))

    unit = x_eval[1] - x_eval[0]

    def find_bounds(intersections, unit):
        bounds = []
        intersections.append(np.inf)
        dif_list = [intersections[i + 1] - intersections[i] for i in range(len(intersections) - 1)]
        min_val = intersections[0]
        for i in range(len(dif_list)):
            if dif_list[i] < 2 * unit:
                continue
            else:
                max_val = intersections[i]
                bounds.append((min_val, max_val))
                min_val = intersections[i + 1]

        return bounds

    return find_bounds(intersections, unit)

