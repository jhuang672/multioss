# Multiple-output calibration of a honeycomb seal via on-site surrogates

This repository contains the data and R code to implement the fully integrated inference for calibration parameter by combining all outputs through the fully Bayesian approach. The sampled posterior distribution was marginally plotted in diagonal and lower triangle of Figure 8 in the paper. 

### What is in this directory? 

Data: 
* **field.csv**: the full field data of honeycomb physical experiments, including inputs and multiple outputs at multiple frequencies.
* **simulation.csv**: the full simulation data from the isotseal computer model, including inputs, parameters, and multiple outputs at multiple frequencies. 

Code: 

* **bayes.R**: combine all 16 outputs in honeycomb in fully Bayesian MCMC to infer unknown calibration parameter under Beta (2, 2) prior for U.

Note: this R code requires high performance computing environments with Intel's Math Kernel Library (MKL) setup with multiple threads. 


### Who do I talk to? ###

* This repository is actively maintained by Jiangeng Huang <huangj@vt.edu>.


