import numpy as np
from sklearn.metrics import silhouette_score
from sklearn.neighbors import kneighbors_graph
from scipy.spatial import distance
import networkx as nx

def cell_type_matching_accuracy(m1_source_ct, m1_predict_ct, m2_source_ct, m2_predict_ct):
    """
    Calculate the cell-type prediction accuracy of cell-cell alignment.

    Parameters
    ----------
    m1_source_ct : list
        The original cell type labels of modality 1.
    m1_predict_ct : list
        The predicted cell type labels of modality 1.
    m2_source_ct : list
        The original cell type labels of modality 2.
    m2_predict_ct : list
        The predicted cell type labels of modality 2.

    Returns
    -------
    float
        The overall cell-type prediction accuracy (rounded to 4 decimals).
    """
    n1 = len(m1_source_ct)
    n2 = len(m2_source_ct)
    r1 = sum(m1_source_ct == m1_predict_ct)
    r2 = sum(m2_source_ct == m2_predict_ct)
    acc = (r1 + r2) / (n1 + n2)
    acc = round(acc, 4)
    return acc

def average_sihouette_width(embedding, cell_type_label):
    """
    Calculate the average silhouette score of integration performance.

    Parameters
    ----------
    embedding : np.ndarray of shape (n_samples, n_dims)
        The 2D (or low-dim) embedding array used for silhouette calculation.
    cell_type_label : array-like of shape (n_samples,)
        The cell type labels corresponding to each embedding point.

    Returns
    -------
    float
        The average silhouette score (rounded to 4 decimals).
    """
    sihouette_avg = silhouette_score(embedding, cell_type_label) 
    sihouette_avg = round(sihouette_avg, 4)
    return sihouette_avg

def feature_imputation_accuracy_corr(m1_feature, m2_aligned_feature1):
    """
    Calculate the feature imputation accuracy of the aligned feature profile using Pearson correlation.

    Parameters
    ----------
    m1_feature : np.ndarray of shape (n_samples, n_features)
        The original modality-1 feature matrix.
    m2_aligned_feature1 : np.ndarray of shape (n_samples, n_features)
        The imputed (aligned) modality-1 feature matrix obtained for modality 2.

    Returns
    -------
    float
        The average per-sample Pearson correlation (imputation accuracy).
    """
    assert m1_feature.shape == m2_aligned_feature1.shape
    n_samples = m1_feature.shape[0]
    corr_vec = np.zeros(n_samples)
    for i in range(n_samples):
        corr = np.corrcoef(m1_feature[i, :], m2_aligned_feature1[i, :])[0, 1]
        corr_vec[i] = round(corr, 4)
    impute_acc = np.mean(corr_vec)
    return impute_acc

def feature_imputation_rmse(m1_feature, m2_aligned_feature1):
    """
    Calculate the RMSE of the aligned feature profile.

    Parameters
    ----------
    m1_feature : np.ndarray of shape (n_samples, n_features)
        The original modality-1 feature matrix.
    m2_aligned_feature1 : np.ndarray of shape (n_samples, n_features)
        The imputed (aligned) modality-1 feature matrix obtained for modality 2.

    Returns
    -------
    float
        The root mean squared error between the original and imputed features.
    """
    assert m1_feature.shape == m2_aligned_feature1.shape
    error = m1_feature - m2_aligned_feature1
    squared_error = np.square(error)
    mean_squared_error = np.mean(squared_error)
    impute_rmse = np.sqrt(mean_squared_error)
    return impute_rmse

def uniFOSCTTM(m1_embedding, m2_embedding, true_matches_for_m2):
    """
    Calculate the proportion of samples closer than the true paired sample (uniFOSCTTM).

    Parameters
    ----------
    m1_embedding : np.ndarray of shape (n, d)
        Embedding of modality 1.
    m2_embedding : np.ndarray of shape (n, d)
        Embedding of modality 2.
    true_matches_for_m2 : array-like of length n
        Indices of the true matched cells in modality 1 for each cell in modality 2.

    Returns
    -------
    float
        The uniFOSCTTM score (rounded to 4 decimals).
    """
    distance_matrix = distance.cdist(m2_embedding, m1_embedding, metric = 'euclidean')
    n = len(true_matches_for_m2)
    vec = np.zeros(n)
    for idx, true_match in enumerate(true_matches_for_m2):
        true_distance = distance_matrix[idx, true_match]
        # Count how many cells in modality 1 are closer to cell idx in modality 2 than the true match
        closer_samples = np.sum(distance_matrix[idx, :] < true_distance)
        vec[idx] = closer_samples / distance_matrix.shape[1]
        prop = np.mean(vec)

    return round(prop, 4)


def calculate_graph_connectivity(data, labels, k=15):
    """
    Calculate the Graph Connectivity for each cell type in the dataset.

    Parameters
    ----------
    data : np.ndarray of shape (n_samples, n_features)
        The dataset where rows are samples and columns are features.
    labels : array-like of shape (n_samples,)
        The cell type labels for each sample.
    k : int, default=15
        Number of nearest neighbors to consider for each cell.

    Returns
    -------
    float
        The graph connectivity score averaged across cell types.
    """
    kng = kneighbors_graph(data, n_neighbors=k, mode='connectivity', include_self=False)
    G = nx.from_scipy_sparse_array(kng)

    unique_labels = np.unique(labels)
    M = len(unique_labels)
    sum_lcc_ratio = 0

    for label in unique_labels:
        indices = np.where(labels == label)[0]
        subG = G.subgraph(indices)
        largest_cc = max(nx.connected_components(subG), key=len)
        LCC_j = len(largest_cc)
        N_j = len(indices)
        lcc_ratio = LCC_j / N_j
        sum_lcc_ratio += lcc_ratio

    return sum_lcc_ratio / M
