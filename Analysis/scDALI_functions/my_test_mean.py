"""Wrapper function to run tests implemented in this package."""

from functools import partial

import numpy as np
import sys
sys.path.append("/Users/jasper/Desktop/PhD/Projects/ASE_Spermatogenesis/Scripts/dali/")

import dali
from dali.models import *
from dali.utils.matop import atleast_2d_column

def test_mean(a, d, mean_mu = 0.5):
	a = atleast_2d_column(a)
	d = atleast_2d_column(d)
	bbm_model = BetaBinomMean(a, d, mean_mu)
	bbm_model.fit()
	pval = bbm_model.test()
	results = [bbm_model.mu1, bbm_model.theta1, bbm_model.nll0, bbm_model.nll1, pval]
	return(results)
