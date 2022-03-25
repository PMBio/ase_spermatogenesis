"""Wrapper function to run tests implemented in this package."""


from functools import partial

import numpy as np
import sys
sys.path.append("/Users/jasper/Desktop/PhD/Projects/ASE_Spermatogenesis/Scripts/dali/")

import dali

from dali.models import *
from dali.utils.run_model import run_model, callback_test
from dali.utils.parallel import process_parallel
from dali.utils.matop import atleast_2d_column



MODELS = {
    'daliBB': BetaBinom,
    'ttest': ClusterTTest,
    'meanBB': BetaBinomMean
}

def test_regions(
        A, D, cell_state=None,
        test='daliBB',
        null_mean=None,
        n_cores=1):
    """Run tests for each region.

    A, D are assumed to be n-by-d, where n is the number of cells and d the
    number of regions to test.

    Args:
        A: Alternative counts for each cell and region.
        D: Total counts for each cell and region.
        cell_state: Matrix of cell states, e.g. clusters or coordinates
            in a low-dimensional cell-state space. Ignored if test is
            meanBB.
        test: String indicating the test to run. Options are
            'daliBB' - a Beta-Binomial variance component score test to test
                for cell-state-specific variability of allelic rates. The
                cell_state matrix is required here.
            'ttest' - a one-vs-all t-test to test for inter-cluster variability
                in allelic rates. In this case cell_state needs to contain the
                cluster labels.
            'meanBB' - a Beta-Binomial test for deviation from a null mean to
                test for example if allelic ratios differ from .5. Does not
                require cell_state but null_mean to be specified.
        null_mean: Null mean for the meanBB test.
        n_cores: Number of cores to use.

    Returns:
        p-values for each region.
    """

    if A.shape != D.shape:
        raise ValueError('A and D need to be of the same shape.')
    if (cell_state is not None) and (cell_state.shape[0] != A.shape[0]):
        raise ValueError('First dimensions of A and cell_state need to match.')

    try:
        model = MODELS[test]
    except KeyError:
        msg = ('Model not recognized. Choices are '
            ', '.join(MODELS.keys()) + '.')
        raise ValueError(msg)

    if test in ['daliBB', 'ttest'] and cell_state is None:
        raise ValueError('%s requires cell_state to be specified' % test)
    if test == 'meanBB' and null_mean is None:
        raise ValueError('meanBB requires null_mean to be specified')

    D = atleast_2d_column(D)
    A = atleast_2d_column(A)
    cell_state = atleast_2d_column(cell_state)

    n_cores = min(n_cores, D.shape[1])
    #print('[dali] Testing %d regions on %d core(s) ... ' % (D.shape[1], n_cores), flush=True)

    init_kwargs = {} if null_mean is None else {'null_mean': null_mean}
    show_progress = False if n_cores > 1 else True

    f = partial(
            run_model,
            model,
            E=cell_state,
            init_kwargs = init_kwargs,
            callbacks=[callback_test],
            show_progress=show_progress)

    results = process_parallel(
            f,
            mat_dict={'A':A, 'D':D},
            n_cores=n_cores)

    return np.asarray([r[0] for r in results])



