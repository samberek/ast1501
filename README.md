# ast1501

## Estimating the mass of NGC 1407 using a bayesian model with Stan

We use a distribution function analysis to estimate the enclosed mass of NGC 1407, using globular clusters as the tracer population. Details of the choice of distribution are below.

The Stan code for the model is in [single_hernquist.stan](single_hernquist.stan). This model is a single Hernquist model, with one Hernquist function modeling both baryonic and dark matter. 

The Stan code is run in python, using pystan. This is found in [pystan_ngc1407.ipynb](pystan_ngc1407.ipynb). 

The last notebook, [prelim_1407.ipynb](prelim_1407.ipynb), contains some preliminary analysis of the globular cluster data, as well as a plot of the distribution function used in the Stan model. 

