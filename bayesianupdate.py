#!/usr/bin/env python
# coding: utf-8

## Bayesian Update module

# Import libraries
import numpy as np

# Define the rule for a iterative bayesian update step
def update(prior_mean, prior_sigma, measure_mean, measure_sigma, n_measures):
	# Update step 
	updated_mean = ((measure_sigma**2)*prior_mean + n_measures*measure_mean*(prior_sigma**2)) / (n_measures*(prior_sigma**2) + (measure_sigma**2))
	updated_sigma = (measure_sigma**2)*(prior_sigma**2) / (n_measures*(prior_sigma**2) + (measure_sigma**2))
	updated_sigma = np.sqrt(updated_sigma)
	# Return updated measures
	return updated_mean, updated_sigma
