# Code for "Algorithms for Models with Intractable Normalizing Functions" in Handbook of Markov Chain Monte Carlo
Authors: Murali Haran, Bokgyeong Kang, and Jaewoo Park

We provide instructions for implementing some algorithms for models with intractable normalizing functions and testing the quality of samples. 


## 1. A Potts model
Before running any code, ensure the required R packages have been installed. Set the R working directory to `/potts`.

### Required packages
The code has been tested with R version 4.2.2, "Innocent and Trusting."  The following R packages must be installed before the code will run successfully:

- `Rcpp`
- `RcppArmadillo`
- `qrng`
- `foreach`
- `coda`
- `batchmeans`
- `tidyverse`
- `potts`
- `DiceKriging`
- `DiceDesign`
- `fields`


### Simulate data
`/potts/data.R`

- Simulate a 30 by 30 dataset
- All components are saved in `/potts/data/`


### Run DMH, ABC-MCMC, and LikEm algorithms
`/potts/postSamp_dmh.R` `/potts/postSamp_abcmcmc.R` `/potts/postSamp_likem.R`

- Generate posterior samples for model parameters using the algorithms 
- Posterior samples are saved in `/potts/postSamp/`


### Testing quality of samples

- `aux.R`: Sample a number of particles over a parameter space and generate auxiliary variables for each particle. The auxiliary variables are saved in `\potts\aux\`
- `appx_dmh.R` `appx_abcmcmc.R` `appx_likem.R`: Approximate the posterior's score function $u(\theta)$ and the half-vectorization of the sum of sensitivity and variability matrices $d(\theta) = vech[J(\theta) + H(\theta)]$ for each posterior sample. The approximations are saved in `\potts\appx\`
- `acd_dmh.R` `acd_abcmcmc.R` `acd_likem.R`: Compute an approximate curvature diagnostic (ACD) for each posterior sample path. The diagnostic values are saved in `\potts\acd\`



## 2. An exponential random graph model (ERGM)
Before running any code, ensure the required R packages have been installed. Set the R working directory to `/ergm`.

### Required packages
The code has been tested with R version 4.2.2, "Innocent and Trusting."  The following R packages must be installed before the code will run successfully:

- `Rcpp`
- `RcppArmadillo`
- `qrng`
- `foreach`
- `coda`
- `batchmeans`
- `tidyverse`
- `ergm`
- `fields`
- `Bergm`
- `Matrix`
- `rlist`
- `MASS`
- `spatstat`
- `statmod`

The `Bergm` works with versions 3.7.1 and 3.8.0 of `ergm` which can be installed as follows:
```s
require(devtools)
install_version("ergm", version = "3.8.0", repos = "http://cran.us.r-project.org")
```

### Load data
`/ergm/data.R`

- Load the Florentine marriage dataset (Breiger and Pattison, 1986)
- All components are saved in `/ergm/data/`


### Run ALR, DMH, ABC-MCMC, and VI algorithms
`/ergm/postSamp_alr.R` `/ergm/postSamp_dmh.R` `/ergm/postSamp_abcmcmc.R` `/ergm/postSamp_vi.R`

- Generate posterior samples for model parameters using the algorithms 
- Posterior samples are saved in `/ergm/postSamp/`


### Testing quality of samples

- `aux.R`: Sample a number of particles over a parameter space and generate auxiliary variables for each particle. The auxiliary variables are saved in `\ergm\aux\`
- `appx_alr.R` `appx_dmh.R` `appx_abcmcmc.R` `appx_vi.R`: Approximate the posterior's score function $u(\theta)$ and the half-vectorization of the sum of sensitivity and variability matrices $d(\theta) = vech[J(\theta) + H(\theta)]$ for each posterior sample. The approximations are saved in `\ergm\appx\`
- `acd_alr.R` `acd_dmh.R` `acd_abcmcmc.R` `acd_vi.R`: Compute an approximate curvature diagnostic (ACD) for each posterior sample path. The diagnostic values are saved in `\ergm\acd\`
