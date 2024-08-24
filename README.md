# Code for "Algorithms for Models with Intractable Normalizing Functions" in Handbook of Markov Chain Monte Carlo
Authors: Murali Haran, Bokgyeong Kang, and Jaewoo Park

We provide instructions for implementing some algorithms for models with intractable normalizing functions and testing the quality of samples. 


## Required packages:
The code has been tested with R version 4.2.2, "Innocent and Trusting."  The following R packages must be installed before the code will run successfully:

- `Rcpp`
- `RcppArmadillo`
- `sitmo`
- `akima`
- `mpcmp`
- `qrng`
- `ngspatial`
- `foreach`
- `coda`
- `batchmeans`
- `tidyverse`

The `mpcmp` can be installed from [GitHub](https://github.com/thomas-fung/mpcmp):
```s
require(devtools)
devtools::install_github("thomas-fung/mpcmp")
```

## Spatial mean-COMP regression model
Before running any code, ensure the required R packages have been installed. Set the R working directory to `/scomp`.

### Pre-computation for spline approximation
`/scomp/appx.R`
- Generate a number of particles over $\Psi = [\log(0.01), \log(\mu_{\max})] \times [0.01, \nu_{\max}]$ and find $\lambda(\mu, \nu)$ for each particle using the Newtwon-Raphson algorithm
- This step corresponds to Section 3.1 of the manuscript
- All components are saved in the file `/scomp/appx/set.RData`

### Simulate data
`/scomp/data.R`
- Generate data from the spatial mean-COMP regression model
- All components are saved in the file `/scomp/data/sim.RData`

### Fit the model and summarize the results
`/scomp/fitCOMP.R`
- Fit the spatial mean-COMP regression model to the simulated dataset
- Posterior samples are saved in the file `/scomp/fit/simCOMP.RData`
- Obtain summary statistics for model parameters


## Spatial mean-ZICOMP regression model
Before running any code, ensure the required R packages have been installed. Set the R working directory to `/szicomp`.

### Pre-computation for spline approximation
`/szicomp/appx.R`
- Generate a number of particles over $\Psi = [\log(0.01), \log(\mu_{\max})] \times [0.01, \nu_{\max}]$ and find $\lambda(\mu, \nu)$ for each particle using the Newtwon-Raphson algorithm
- This step corresponds to Section 3.1 of the manuscript
- All components are saved in the file `/szicomp/appx/set.RData`

### Simulate data
`/szicomp/data.R`
- Generate data from the spatial mean-COMP regression model
- All components are saved in the file `/szicomp/data/sim.RData`

### Fit the model and summarize the results
`/szicomp/fitZICOMP.R`
- Fit the spatial mean-ZICOMP regression model to the simulated dataset
- Posterior samples are saved in the file `/szicomp/fit/simZICOMP.RData`
- Obtain summary statistics for model parameters

