# Multiple-output calibration of a honeycomb seal via on-site surrogates

This repository contains the code and documentation to implement the fully integrated inference for calibration parameter by combining all outputs through optimization and full Bayes in Figure 8 of the paper. 

### What is in this directory? 


* **multi_oss_cali.R**: combine all 16 outputs in honeycomb in modularized optization to infer unknown calibration parameter under Beta (2, 2) prior for U.
* **multi_oss_bayes.R**: combine all 16 outputs in honeycomb in fully Bayesian MCMC to infer unknown calibration parameter under Beta (2, 2) prior for U.

Note: both code require high performance computing environments with multiple threads. 


### Who do I talk to? ###

* This repository is actively maintained by Jiangeng Huang <huangj@vt.edu>


