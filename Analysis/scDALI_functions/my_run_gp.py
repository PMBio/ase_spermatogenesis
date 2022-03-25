"""Wrapper function to run GP interpolation."""


import numpy as np
import sys
sys.path.append("/Users/jasper/Desktop/PhD/Projects/ASE_Spermatogenesis/Scripts/dali/")

import dali


from functools import partial

from dali.models import SparseGP
from dali.utils.run_model import *
from dali.utils.parallel import process_parallel
from dali.utils.matop import preprocess_clusters


def run_gp(
        A, D, cell_state,
        kernel='Linear',
        num_inducing=800,
        max_iter=2000,
        n_cores=1,
        return_prior_mean=False,
        return_full_cov=False,
        cluster_ids=None):
    """Train a Gaussian process model for each region.

    Infers posterior mean and variances for the allelic rate in each cell which
    can be used for visualization of cell-state-specific effects and effect
    size estimation.

    A, D are assumed to be n-by-d, where n is the number of cells and d the
    number of regions to test.

    Args:
        A: Alternative counts for each cell and region.
        D: Total counts for each cell and region.
        cell_state: Matrix of cell states, e.g. clusters or coordinates
            in a low-dimensional cell-state space.
        kernel: String describing kernel function. Options are 'RBF' and
            'Linear'.
        num_inducing: Number of inducing points for the sparse approximation.
            Induncing points reduce the computational complexity for fitting a
            GP model from O(n^3) to O(n + m^3) where m are the number of
            inducing points. For m=n, the model equals a classical GP.
        max_iter: Maximum number of iterations.
        n_cores: Number of cores to use.
        return_prior_mean: Return the prior mean.
        return_full_cov: Return full posterior covariances. Not recommended
            unless the number of cells or the number of regions are very small.
        cluster_ids: If not None, also return aggregate posterior estimates for
            each cluster by computing the mean and variance of the arithmetic
            mean.

    Returns:
        Dict with model model posterior and optional return values.
    """
    if A.shape != D.shape:
        raise ValueError('A and D need to be of the same shape.')
    if cell_state.shape[0] != A.shape[0]:
        raise ValueError('First dimensions of A and cell_state need to match.')
    if (cluster_ids is not None) and (cluster_ids.shape[1] != A.shape[0]):
        raise ValueError('First dimensions of A and cluster_ids need to match.')


    D = atleast_2d_column(D)
    A = atleast_2d_column(A)
    cell_state = atleast_2d_column(cell_state)

    n_cores = min(n_cores, D.shape[1])
    print('[dali] Running GP for %d regions on %d core(s) ... ' % (D.shape[1], n_cores), flush=True)

    show_progress = False if n_cores > 1 else True

    group_ids = None
    if cluster_ids is not None:
        group_ids = preprocess_clusters(cluster_ids)
    callbacks = [partial(
        callback_latent_posterior,
        E=cell_state,
        group_ids=cluster_ids,
        full_cov=return_full_cov)]
    if return_prior_mean:
        callbacks = callbacks + [callback_model_mean]

    f = partial(
        run_model,
        SparseGP,
        E=cell_state,
        init_kwargs={'kernel': kernel, 'num_inducing': num_inducing},
        fit_kwargs={'maxiter': max_iter},
        callbacks=callbacks,
        show_progress=show_progress)

    results = process_parallel(
            f,
            mat_dict={'A':A, 'D':D},
            n_cores=n_cores)

    out = dict()
    out['posterior_mean'] = np.asarray([r[0][0].flatten() for r in results]).T
    if return_full_cov:
        out['posterior_cov'] = [r[0][1] for r in results]
    else:
        out['posterior_var'] = np.asarray([r[0][1].flatten() for r in results]).T
    if cluster_ids is not None:
        out['posterior_cluster_mean'] = np.asarray([r[0][2].flatten() for r in results])
        out['posterior_cluster_var'] = np.asarray([r[0][3].flatten() for r in results])
    if return_prior_mean:
        out['prior_mean'] = [float(r[1]) for r in results]
    return out

