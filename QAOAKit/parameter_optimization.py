import pickle
import numpy as np
from pathlib import Path
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV

from QAOAKit import get_full_qaoa_dataset_table

parameter_optimization_folder = Path(__file__).parent


def train_kde(p, n, n_jobs=1, bandwidth_range=np.logspace(-2, 1, 20)):
    """
    Fit KDE to optimized parameters for a given p <= 3 and n <= 9
    Resulting KDE can be used to sample optimized parameters for QAOA on MaxCut
    Follows the methodology of https://doi.org/10.1609/aaai.v34i03.5616

    Parameters
    ----------
    p : int
        Number of QAOA layers p
    n : int
        Number of nodes
        Optimal angles for all non-isomorphic graphs on n nodes are used
    n_jobs : int
        Number of jobs to use to perform cross validation
        for kernel bandwidth optimization
    bandwidth_range : list-like
        Values of bandwidth to consider

    Returns
    -------
    median, kde : tuple(np.array, sklearn.neighbors.KernelDensity)
        Tuple of median angles and fitted kernel density model
    """
    df = get_full_qaoa_dataset_table().reset_index().set_index("graph_id")
    df = df[(df["p_max"] == p) & (df["n"] == n)]
    df["average degree"] = df.apply(
        lambda row: 2 * row["G"].number_of_edges() / row["G"].number_of_nodes(), axis=1
    )

    if p == 1:
        data = df.apply(
            lambda x: np.hstack(
                [
                    np.array(x["gamma"])
                    / np.arctan(1 / np.sqrt(x["average degree"] - 1)),
                    np.array(x["beta"]),
                ]
            ),
            axis=1,
        ).values
    else:
        data = df.apply(
            lambda x: np.hstack(
                [
                    np.array(x["gamma"]) * np.sqrt(x["average degree"]),
                    np.array(x["beta"]),
                ]
            ),
            axis=1,
        ).values
    data = np.stack(data)
    median = np.median(data, axis=0)

    print(f"Fitting a KDE model on data of shape {data.shape}")
    # use grid search cross-validation to optimize the bandwidth
    params = {"bandwidth": bandwidth_range}
    grid = GridSearchCV(KernelDensity(), params, n_jobs=n_jobs)
    grid.fit(data)

    print(
        f"best bandwidth: {grid.best_estimator_.bandwidth}\n minimum tried {min(bandwidth_range)}\n maximum tried {max(bandwidth_range)}"
    )

    # use the best estimator to compute the kernel density estimate
    kde = grid.best_estimator_

    return median, kde


def get_median_pre_trained_kde(p):
    """
    Returns pre-fitted KDE with optimized parameters for a given p <= 3
    KDE is fitted on optimal parameters for all non-isomorphic graphs with n=9
    Resulting KDE can be used to sample optimized parameters for QAOA on MaxCut
    Follows the methodology of https://doi.org/10.1609/aaai.v34i03.5616

    Parameters
    ----------
    p : int
        Number of QAOA layers p
    n : int
        Number of nodes
        Optimal angles for all non-isomorphic graphs on n nodes are used

    Returns
    -------
    median, kde : tuple(np.array, sklearn.neighbors.KernelDensity)
        Tuple of median angles and fitted kernel density model
    """
    kde_path = Path(
        parameter_optimization_folder,
        f"../data/pretrained_models/kde_n=9_p={p}_large_bandwidth_range.p",
    )
    try:
        return pickle.load(open(kde_path, "rb"))
    except pickle.PickleError:
        raise pickle.PickleError(
            f"Failed to unpickle the pre-trained KDE at {kde_path}. Please re-train the model using QAOAKit.parameter_optimization.train_kde."
        )
