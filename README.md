# ast1501

## Estimating the mass of NGC 1407 using a bayesian model with Stan

We use a distribution function analysis to estimate the enclosed mass of NGC 1407, using globular clusters as the tracer population. Details of the choice of distribution are below.

The Stan code for the model is in [single_hernquist.stan](single_hernquist.stan). This model is a single Hernquist model, with one Hernquist function modeling both baryonic and dark matter. 

The Stan code is run in python, using pystan. This is found in [pystan_ngc1407.ipynb](pystan_ngc1407.ipynb). 

The notebook, [prelim_1407.ipynb](prelim_1407.ipynb) contains some preliminary analysis of the globular cluster data, as well as a plot of the distribution function used in the Stan model. 

The globular cluster data is in [NGC1407_GC.cat](NGC1407_GC.cat). Radial velocity measurements are in km/s, while radius values (projected distance to the center of the galaxy) are in arcmin. 


The Stan code is tested with simulated data in [simulated_data.stan](simulated_data.stan) and run in [simulated_data.ipynb](simulated_data.ipynb). Full, 3D positions and velocities of globular clusters are generated assuming a total galaxy mass of logM = 12.5. As of now, this simulated data also causes the Stan code to diverge most of the time, suggesting that the problem is with the code and not the data. 