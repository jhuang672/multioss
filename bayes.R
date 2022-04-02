
########################################################################################################
## Fully Bayesian calibration using on-site surrogates on all 16 outputs together
## Combining 16 frequencies into PCs for each output
## Implement the Bayesian joint inference method described in Section 5.1
## Reproduce the fully Bayesian integrated results in Figure 8
########################################################################################################

## Content of this code: this fully Bayesian implementation has been broken down into 4 steps
##
## ---------------------------------------------------------------------------------------
##                      Step 1. Setup before posterior sampling for U
## 
##                      Step 1.1 Read and organize data
##                      Step 1.2 PC rotation on whole U space for all observed discrepancy
##                      Step 1.3 Build new on-site surrogates on the rotated new space
##                      Step 1.4 Pre-calculation for on-site kernels and their inverses
##                      Step 1.5 Train and check the fitted discrepancies for multiple PCs 
## ---------------------------------------------------------------------------------------
##                      Step 2. Calculate the posterior: initial evaluation
## 
##                      Step 2.1. inverse part
##                      Step 2.2. determinant part
##                      Step 2.3. log posterior: put all parts together
## ---------------------------------------------------------------------------------------
##                      Step 3. MCMC Sampling 
##
## ---------------------------------------------------------------------------------------
##                      Step 4. Collect samples and save objects  
##
## ---------------------------------------------------------------------------------------


## ------------------------------------------------------------------------
##                      Step 1. Setup before posterior sampling for U
## ------------------------------------------------------------------------

## Set up required packages 
library(laGP)        
library(plgp)        
library(factoextra)  

## ------------------------------------------------------------------------
##                      Step 1.1 Read and organize data:
## ------------------------------------------------------------------------


# Read in field data: including all inputs Xf and multiple outputs 
Dfield <- read.csv("field.csv")

# Scale field responses: 
# field responses Cdir and Ccross need to be scaled: 
Dfield$Cxx_2 <- Dfield$Cxx_2 * 10^3
Dfield$Cxx_5 <- Dfield$Cxx_5 * 10^3
Dfield$Cxx_7 <- Dfield$Cxx_7 * 10^3
Dfield$Cxx_9 <- Dfield$Cxx_9 * 10^3
Dfield$Cxx_11 <- Dfield$Cxx_11 * 10^3
Dfield$Cxx_13 <- Dfield$Cxx_13 * 10^3

Dfield$c_xy_2 <- Dfield$c_xy_2 * 10^3
Dfield$c_xy_5 <- Dfield$c_xy_5 * 10^3
Dfield$c_xy_7 <- Dfield$c_xy_7 * 10^3
Dfield$c_xy_9 <- Dfield$c_xy_9 * 10^3
Dfield$c_xy_11 <- Dfield$c_xy_11 * 10^3
Dfield$c_xy_13 <- Dfield$c_xy_13 * 10^3


# Read in isotseal data: including all inputs X, U, and multiple outputs 
Dsim <- read.csv("simulation.cvs") # 1000 maximinLHS design for all the same U at each Xfield locations
Dsim <- Dsim[, -1] # remove the redundant first column

# scale all simulated outputs: 
## Kdir: 
Dsim$Kdir   <- Dsim$Kdir / 10 ^ 6  # scale response to field unit
Dsim$Kdir5  <- Dsim$Kdir5 / 10 ^ 6  # scale response to field unit
Dsim$Kdir7  <- Dsim$Kdir7 / 10 ^ 6  # scale response to field unit
Dsim$Kdir9  <- Dsim$Kdir9 / 10 ^ 6  # scale response to field unit
Dsim$Kdir11 <- Dsim$Kdir11 / 10 ^ 6  # scale response to field unit
Dsim$Kdir13 <- Dsim$Kdir13 / 10 ^ 6  # scale response to field unit

## Kcross: 
Dsim$Kcross   <- Dsim$Kcross / 10 ^ 6  # scale response to field unit
Dsim$Kcross5  <- Dsim$Kcross5 / 10 ^ 6  # scale response to field unit
Dsim$Kcross7  <- Dsim$Kcross7 / 10 ^ 6  # scale response to field unit
Dsim$Kcross9  <- Dsim$Kcross9 / 10 ^ 6  # scale response to field unit
Dsim$Kcross11 <- Dsim$Kcross11 / 10 ^ 6  # scale response to field unit
Dsim$Kcross13 <- Dsim$Kcross13 / 10 ^ 6  # scale response to field unit

## Cdir: 
Dsim$Cdir   <- Dsim$Cdir / 10 ^ 3  # scale response to field unit
Dsim$Cdir5  <- Dsim$Cdir5 / 10 ^ 3  # scale response to field unit
Dsim$Cdir7  <- Dsim$Cdir7 / 10 ^ 3  # scale response to field unit
Dsim$Cdir9  <- Dsim$Cdir9 / 10 ^ 3  # scale response to field unit
Dsim$Cdir11 <- Dsim$Cdir11 / 10 ^ 3  # scale response to field unit
Dsim$Cdir13 <- Dsim$Cdir13 / 10 ^ 3  # scale response to field unit

## Cdir: 
Dsim$Ccross   <- Dsim$Ccross / 10 ^ 3  # scale response to field unit
Dsim$Ccross5  <- Dsim$Ccross5 / 10 ^ 3  # scale response to field unit
Dsim$Ccross7  <- Dsim$Ccross7 / 10 ^ 3  # scale response to field unit
Dsim$Ccross9  <- Dsim$Ccross9 / 10 ^ 3  # scale response to field unit
Dsim$Ccross11 <- Dsim$Ccross11 / 10 ^ 3  # scale response to field unit
Dsim$Ccross13 <- Dsim$Ccross13 / 10 ^ 3  # scale response to field unit


## Organize multiple outputs: 

# for  "Kdir"
## for simulation:
ym_1_2 <- Dsim$Kdir
ym_1_5 <- Dsim$Kdir5
ym_1_9 <- Dsim$Kdir9
ym_1_11 <- Dsim$Kdir11

## for field outputs: 
YF_1 <- data.frame(yf_1_2 = Dfield$Kxx_2, yf_1_5 = Dfield$Kxx_5, yf_1_9 = Dfield$Kxx_9, yf_1_11 = Dfield$Kxx_11)

# for "Kcross"
## for simulation: 
ym_2_2 <- Dsim$Kcross
ym_2_5 <- Dsim$Kcross5
ym_2_9 <- Dsim$Kcross9
ym_2_11 <- Dsim$Kcross11

## for field outputs: 
YF_2 <- data.frame(yf_2_2 = Dfield$k_xy_2, yf_2_5 = Dfield$k_xy_5, yf_2_9 = Dfield$k_xy_9, yf_2_11 = Dfield$k_xy_11)


# for "Cdir"
## for simulation: 
ym_3_2 <- Dsim$Cdir
ym_3_5 <- Dsim$Cdir5
ym_3_9 <- Dsim$Cdir9
ym_3_11 <- Dsim$Cdir11

## For field outputs: 
YF_3 <- data.frame(yf_3_2 = Dfield$Cxx_2, yf_3_5 = Dfield$Cxx_5, yf_3_9 = Dfield$Cxx_9, yf_3_11 = Dfield$Cxx_11)


# for "Ccross"
## for simulation: 
ym_4_2 <- Dsim$Ccross
ym_4_5 <- Dsim$Ccross5
ym_4_9 <- Dsim$Ccross9
ym_4_11 <- Dsim$Ccross11

## for field outputs: 
YF_4 <- data.frame(yf_4_2 = Dfield$c_xy_2, yf_4_5 = Dfield$c_xy_5, yf_4_9 = Dfield$c_xy_9, yf_4_11 = Dfield$c_xy_11)


# create an index list for NAs in each set of 1000 run simulation data:
list_na_index <- list()

for (i in 1:292){  
  list_na_index[[i]] <- which(is.na(ym_2_2[1:1000 + 1000 * (i-1)]) == FALSE)
}   


# Set up X and U names: 
xnames <- c("Speed", "Pin", "Pout", "Temp", "Swirl", "viscosity", "gamma", "zeta", "Cin", "Cout", "D", "L", "Cd")
unames <- c("ffactor_Ns", "ffactor_Ms","ffactor_Nr", "ffactor_Mr")

# Set up U limits
eps <- sqrt(.Machine$double.eps)
lower_limits <- c(0.0028, -.15, .0021, -.2956)
upper_limits <- c(0.1, -eps, .0746, -eps) 
urange <- rbind(lower_limits, upper_limits)
colnames(urange) <- unames


# Scale both inputs X and U into (0, 1) unit cube: 
XUori <- Dsim[, 1:17]
Xori <- Dfield[, xnames]
Ulimits <- rbind(lower_limits, upper_limits)
colnames(Ulimits) <- unames
limits <- cbind(apply(Dfield[c(xnames)], 2, range), Ulimits)
XU <- matrix(NA, nrow = nrow(XUori), ncol = ncol(XUori))


# scale XU for isotseal: only 1000 runs need to be scaled due to the same design for U
for (j in 1:ncol(XUori)) {
  XU[ , j] <- (XUori[ , j] - limits[1, j]) / (limits[2, j] - limits[1, j]) 
}


# scale U for isotseal: each 292 has exactly the same 1000 Us: 
Uori <- Dsim[1:1000, 14:17]
U <- matrix(NA, nrow = nrow(Uori), ncol = ncol(Uori))


# scale XU for isotseal: only 1000 runs need to be scaled 
for (j in 1:4) {
  U[ , j] <- (Uori[ , j] - limits[1, 13+j]) / (limits[2, 13+j] - limits[1, 13+j]) 
}

## Check the space filling of parameter U: 
summary(U)
pairs(U)

# scale X for field data
X <- matrix(NA, nrow = nrow(Xori), ncol = ncol(Xori))
for (j in 1:ncol(Xori)) {
  X[ , j] <- (Xori[ , j] - limits[1, j]) / (limits[2, j] - limits[1, j]) 
}


## ---------------------------------------------------------------------------------------
##                      Step 1.2 PC rotation on whole U space for all observed discrepancy
## ---------------------------------------------------------------------------------------


## For a fair comparison on the likelihood for optimization: 
## fit a HUGE principal component on the observed discrepancy  ( YF - YM ) on the whole U space first! 

YM_1 <- cbind(ym_1_2, ym_1_5, ym_1_9, ym_1_11)
YM_2 <- cbind(ym_2_2, ym_2_5, ym_2_9, ym_2_11)
YM_3 <- cbind(ym_3_2, ym_3_5, ym_3_9, ym_3_11)
YM_4 <- cbind(ym_4_2, ym_4_5, ym_4_9, ym_4_11)



## Repeat 1000 YF for each field data to obtain discrepancies:

# Rep rows for a matrix: 
rep.row<-function(x,n){
  matrix(rep(x,each=n), nrow = n * nrow(x), ncol = ncol(x))
}

YF_1_rep <- rep.row(as.matrix(YF_1), 1000)
YF_2_rep <- rep.row(as.matrix(YF_2), 1000)
YF_3_rep <- rep.row(as.matrix(YF_3), 1000)
YF_4_rep <- rep.row(as.matrix(YF_4), 1000)

## Save all the observed bias: 
BIAS_1 <- YF_1_rep - YM_1
BIAS_2 <- YF_2_rep - YM_2
BIAS_3 <- YF_3_rep - YM_3
BIAS_4 <- YF_4_rep - YM_4

colnames(BIAS_1) <- c("y2", "y5", "y9", "y11")
colnames(BIAS_2) <- c("y2", "y5", "y9", "y11")
colnames(BIAS_3) <- c("y2", "y5", "y9", "y11")
colnames(BIAS_4) <- c("y2", "y5", "y9", "y11")

# Remove missing: 
BIAS_1 <- na.omit(BIAS_1)
BIAS_2 <- na.omit(BIAS_2)
BIAS_3 <- na.omit(BIAS_3)
BIAS_4 <- na.omit(BIAS_4)

# Perform PCA on the frequencies for all 292,000 observed discrepancy: 
Bias_pr_1 <- prcomp(BIAS_1, scale = TRUE)
Bias_pr_2 <- prcomp(BIAS_2, scale = TRUE)
Bias_pr_3 <- prcomp(BIAS_3, scale = TRUE)
Bias_pr_4 <- prcomp(BIAS_4, scale = TRUE)

# Check the PCA on each outputs: 
summary(Bias_pr_1)
summary(Bias_pr_2)
summary(Bias_pr_3)
summary(Bias_pr_4)

fviz_eig(Bias_pr_1, main = paste("Scree plot of observed bias for kdir"))
fviz_eig(Bias_pr_2, main = paste("Scree plot of observed bias for kcross"))
fviz_eig(Bias_pr_3, main = paste("Scree plot of observed bias for cdir"))
fviz_eig(Bias_pr_4, main = paste("Scree plot of observed bias for ccross"))



## ---------------------------------------------------------------------------------------
##                      Step 1.3 Build on-site surrogates on the rotated new space: 
## ---------------------------------------------------------------------------------------

n <- 1000  

# Train MLEs for the hyperparameters: 
MLEs_Kdir_pc <- data.frame(matrix(NA, nrow=292, ncol= 6))
names(MLEs_Kdir_pc) <- c("theta1", "theta2", "theta3", "theta4", "g", "itera")

MLEs_Kcross_pc <- data.frame(matrix(NA, nrow=292, ncol= 6))
names(MLEs_Kcross_pc) <- c("theta1", "theta2", "theta3", "theta4", "g", "itera")

MLEs_cdir_pc <- data.frame(matrix(NA, nrow=292, ncol= 6))
names(MLEs_cdir_pc) <- c("theta1", "theta2", "theta3", "theta4", "g", "itera")

MLEs_ccross_pc <- data.frame(matrix(NA, nrow=292, ncol= 6))
names(MLEs_ccross_pc) <- c("theta1", "theta2", "theta3", "theta4", "g", "itera")


# Train on-site surrogates on the first PCs and save the MLEs of length-scale and nuggets
for (i in 1:292) {
  
  ## fit new OSSs on first principal components: 
  ## set up the inputs: 
  U_i <- U[list_na_index[[i]], ]
  
  ## set up the outputs: due to missing 
  ## kdir: 
  yi_1_2 <- ym_1_2[1:1000 + 1000 * (i-1)][list_na_index[[i]]] 
  yi_1_5 <- ym_1_5[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_1_9 <- ym_1_9[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_1_11 <- ym_1_11[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  YM_1_i <- cbind(yi_1_2, yi_1_5, yi_1_9, yi_1_11)
  colnames(YM_1_i) <- c("y2", "y5", "y9", "y11")  ## name each column for PCA rotation
  ym_1_i_pc1 <- predict(Bias_pr_1, YM_1_i)[,1]      ## center, standardize, and rotate all yM into the new space
  
  ## kcross:
  yi_2_2 <- ym_2_2[1:1000 + 1000 * (i-1)][list_na_index[[i]]] 
  yi_2_5 <- ym_2_5[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_2_9 <- ym_2_9[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_2_11 <- ym_2_11[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  YM_2_i <- cbind(yi_2_2, yi_2_5, yi_2_9, yi_2_11)
  colnames(YM_2_i) <- c("y2", "y5", "y9", "y11")  ## name each column for PCA rotation
  ym_2_i_pc1 <- predict(Bias_pr_2, YM_2_i)[,1]      ## center, standardize, and rotate all yM into the new space
  
  ## cdir:
  yi_3_2 <- ym_3_2[1:1000 + 1000 * (i-1)][list_na_index[[i]]] 
  yi_3_5 <- ym_3_5[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_3_9 <- ym_3_9[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_3_11 <- ym_3_11[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  YM_3_i <- cbind(yi_3_2, yi_3_5, yi_3_9, yi_3_11)
  colnames(YM_3_i) <- c("y2", "y5", "y9", "y11")  ## name each column for PCA rotation
  ym_3_i_pc1 <- predict(Bias_pr_3, YM_3_i)[,1]      ## center, standardize, and rotate all yM into the new space
  
  ## ccross:
  yi_4_2 <- ym_4_2[1:1000 + 1000 * (i-1)][list_na_index[[i]]] 
  yi_4_5 <- ym_4_5[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_4_9 <- ym_4_9[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_4_11 <- ym_4_11[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  YM_4_i <- cbind(yi_4_2, yi_4_5, yi_4_9, yi_4_11)
  colnames(YM_4_i) <- c("y2", "y5", "y9", "y11")  ## name each column for PCA rotation
  ym_4_i_pc1 <- predict(Bias_pr_4, YM_4_i)[,1]      ## center, standardize, and rotate all yM into the new space
  
  ## set up the priors
  da <- darg(list(mle=TRUE, max=1, min=eps), U_i, samp.size = nrow(U_i)); da # keep all length-scales less than 1 
  ga_1 <- garg(list(mle=TRUE, max=1), ym_1_i_pc1); ga_1    ##  
  ga_2 <- garg(list(mle=TRUE, max=1), ym_2_i_pc1); ga_2    ##  
  ga_3 <- garg(list(mle=TRUE, max=1), ym_3_i_pc1); ga_3    ##  
  ga_4 <- garg(list(mle=TRUE, max=1), ym_4_i_pc1); ga_4    ##  
  
  ## fit and save the new GP fits: 
  gpsep_1_pc1 <- newGPsep(U_i, ym_1_i_pc1, d=rep(da$start, ncol(U_i)), g=ga_1$start, dK=TRUE)
  gpsep_2_pc1 <- newGPsep(U_i, ym_2_i_pc1, d=rep(da$start, ncol(U_i)), g=ga_2$start, dK=TRUE)
  gpsep_3_pc1 <- newGPsep(U_i, ym_3_i_pc1, d=rep(da$start, ncol(U_i)), g=ga_3$start, dK=TRUE)
  gpsep_4_pc1 <- newGPsep(U_i, ym_4_i_pc1, d=rep(da$start, ncol(U_i)), g=ga_4$start, dK=TRUE)
  
  mle_1_pc_i <- mleGPsep(gpsep_1_pc1 , param = "both", tmin=c(da$min, ga_1$min), tmax=c(da$max, ga_1$max))#$theta # MLE for both theta and g
  mle_2_pc_i <- mleGPsep(gpsep_2_pc1 , param = "both", tmin=c(da$min, ga_2$min), tmax=c(da$max, ga_2$max))#$theta # MLE for both theta and g
  mle_3_pc_i <- mleGPsep(gpsep_3_pc1 , param = "both", tmin=c(da$min, ga_3$min), tmax=c(da$max, ga_3$max))#$theta # MLE for both theta and g
  mle_4_pc_i <- mleGPsep(gpsep_4_pc1 , param = "both", tmin=c(da$min, ga_4$min), tmax=c(da$max, ga_4$max))#$theta # MLE for both theta and g
  
  MLEs_Kdir_pc[i, ] <- c(mle_1_pc_i$theta, mle_1_pc_i$its)
  MLEs_Kcross_pc[i, ] <- c(mle_2_pc_i$theta, mle_2_pc_i$its)
  MLEs_cdir_pc[i, ] <- c(mle_3_pc_i$theta, mle_3_pc_i$its)
  MLEs_ccross_pc[i, ] <- c(mle_4_pc_i$theta, mle_4_pc_i$its)
  
  ## Clean up 
  deleteGPsep(gpsep_1_pc1)
  deleteGPsep(gpsep_2_pc1)
  deleteGPsep(gpsep_3_pc1)
  deleteGPsep(gpsep_4_pc1)
  
  print(i) ## check number of sites finished 
}

## check the MLEs: 
summary(MLEs_Kdir_pc)
summary(MLEs_Kcross_pc)
summary(MLEs_cdir_pc)
summary(MLEs_ccross_pc)


## ---------------------------------------------------------------------------------------
##                      Step 1.4 Pre-calculation for on-site matrices and their inverses
## ---------------------------------------------------------------------------------------

## BIG pre-calculation for MCMC: this part only needs to be pre-calculated once
# each of 292 locations have a different length due to a different rate of missing

# Pre-calculate the quantities in the posterior: 
CY_pc_1 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1
CY_pc_2 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1
CY_pc_3 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1
CY_pc_4 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1

Vi_pc_1 <- list()  # combined Vi
Vi_pc_2 <- list()  # combined Vi
Vi_pc_3 <- list()  # combined Vi
Vi_pc_4 <- list()  # combined Vi

Vi_inv_pc_1 <- list()  # combined Vi_inv
Vi_inv_pc_2 <- list()  # combined Vi_inv
Vi_inv_pc_3 <- list()  # combined Vi_inv
Vi_inv_pc_4 <- list()  # combined Vi_inv

Tau2hat_1 <- matrix(NA, nrow = 292, ncol = 1)
Tau2hat_2 <- matrix(NA, nrow = 292, ncol = 1)
Tau2hat_3 <- matrix(NA, nrow = 292, ncol = 1)
Tau2hat_4 <- matrix(NA, nrow = 292, ncol = 1)

## Pre-calculations: 

for (i in 1:292) {
  
  ## set up the outputs: due to missing 
  ## kdir: 
  yi_1_2 <- ym_1_2[1:1000 + 1000 * (i-1)][list_na_index[[i]]] 
  yi_1_5 <- ym_1_5[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_1_9 <- ym_1_9[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_1_11 <- ym_1_11[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  YM_1_i <- cbind(yi_1_2, yi_1_5, yi_1_9, yi_1_11)
  colnames(YM_1_i) <- c("y2", "y5", "y9", "y11")  
  ym_1_i_pc1 <- predict(Bias_pr_1, YM_1_i)[,1]     
  
  ## kcross:
  yi_2_2 <- ym_2_2[1:1000 + 1000 * (i-1)][list_na_index[[i]]] 
  yi_2_5 <- ym_2_5[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_2_9 <- ym_2_9[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_2_11 <- ym_2_11[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  YM_2_i <- cbind(yi_2_2, yi_2_5, yi_2_9, yi_2_11)
  colnames(YM_2_i) <- c("y2", "y5", "y9", "y11")  
  ym_2_i_pc1 <- predict(Bias_pr_2, YM_2_i)[,1]      
  
  ## cdir:
  yi_3_2 <- ym_3_2[1:1000 + 1000 * (i-1)][list_na_index[[i]]] 
  yi_3_5 <- ym_3_5[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_3_9 <- ym_3_9[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_3_11 <- ym_3_11[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  YM_3_i <- cbind(yi_3_2, yi_3_5, yi_3_9, yi_3_11)
  colnames(YM_3_i) <- c("y2", "y5", "y9", "y11")  
  ym_3_i_pc1 <- predict(Bias_pr_3, YM_3_i)[,1]      
  
  ## ccross:
  yi_4_2 <- ym_4_2[1:1000 + 1000 * (i-1)][list_na_index[[i]]] 
  yi_4_5 <- ym_4_5[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_4_9 <- ym_4_9[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_4_11 <- ym_4_11[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  YM_4_i <- cbind(yi_4_2, yi_4_5, yi_4_9, yi_4_11)
  colnames(YM_4_i) <- c("y2", "y5", "y9", "y11")  
  ym_4_i_pc1 <- predict(Bias_pr_4, YM_4_i)[,1]      
  
  ## Fixed covariance part: 
  Cm_pc_1  <- covar.sep(U[list_na_index[[i]], ],  d = MLEs_Kdir_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) # 1000 by 1000
  Cm_pc_2  <- covar.sep(U[list_na_index[[i]], ],  d = MLEs_Kcross_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) # 1000 by 1000
  Cm_pc_3  <- covar.sep(U[list_na_index[[i]], ],  d = MLEs_cdir_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) # 1000 by 1000
  Cm_pc_4  <- covar.sep(U[list_na_index[[i]], ],  d = MLEs_ccross_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) # 1000 by 1000
  
  ## Inverse of on-site designs: 
  Cm_inv_pc_1 <- solve(Cm_pc_1)
  Cm_inv_pc_2 <- solve(Cm_pc_2)
  Cm_inv_pc_3 <- solve(Cm_pc_3)
  Cm_inv_pc_4 <- solve(Cm_pc_4)
  
  ## The fixed whole part: solve(Sigma(U, U)) %*% ym: only need for bias term fitting
  CY_pc_1[[i]] <- Cm_inv_pc_1 %*% ym_1_i_pc1      # 1000 by 1 for all 292 runs
  CY_pc_2[[i]] <- Cm_inv_pc_2 %*% ym_2_i_pc1      # 1000 by 1 for all 292 runs
  CY_pc_3[[i]] <- Cm_inv_pc_3 %*% ym_3_i_pc1      # 1000 by 1 for all 292 runs
  CY_pc_4[[i]] <- Cm_inv_pc_4 %*% ym_4_i_pc1      # 1000 by 1 for all 292 runs
  
  
  ## calculate and save tau2hat's in this step: 
  tau2hat_pc_1 <- drop(t(ym_1_i_pc1) %*%  Cm_inv_pc_1  %*% ym_1_i_pc1 / length(ym_1_i_pc1))   # mle of tau^2
  tau2hat_pc_2 <- drop(t(ym_2_i_pc1) %*%  Cm_inv_pc_2  %*% ym_2_i_pc1 / length(ym_2_i_pc1))   # mle of tau^2
  tau2hat_pc_3 <- drop(t(ym_3_i_pc1) %*%  Cm_inv_pc_3  %*% ym_3_i_pc1 / length(ym_3_i_pc1))   # mle of tau^2
  tau2hat_pc_4 <- drop(t(ym_4_i_pc1) %*%  Cm_inv_pc_4  %*% ym_4_i_pc1 / length(ym_4_i_pc1))   # mle of tau^2
 
  Tau2hat_1[i, ] <- tau2hat_pc_1   # save each local tau 
  Tau2hat_2[i, ] <- tau2hat_pc_2   # save each local tau 
  Tau2hat_3[i, ] <- tau2hat_pc_3   # save each local tau 
  Tau2hat_4[i, ] <- tau2hat_pc_4   # save each local tau 
  
  ##  Vis
  Vi_pc_1[[i]] <-  tau2hat_pc_1 * Cm_pc_1              # save Vis for V(d), each is around 1000 by 1000
  Vi_pc_2[[i]] <-  tau2hat_pc_2 * Cm_pc_2              # save Vis for V(d), each is around 1000 by 1000
  Vi_pc_3[[i]] <-  tau2hat_pc_3 * Cm_pc_3              # save Vis for V(d), each is around 1000 by 1000
  Vi_pc_4[[i]] <-  tau2hat_pc_4 * Cm_pc_4              # save Vis for V(d), each is around 1000 by 1000
  
  ## Vi_inv 
  Vi_inv_pc_1[[i]] <- solve(Vi_pc_1[[i]])  
  Vi_inv_pc_2[[i]] <- solve(Vi_pc_2[[i]])  
  Vi_inv_pc_3[[i]] <- solve(Vi_pc_3[[i]])  
  Vi_inv_pc_4[[i]] <- solve(Vi_pc_4[[i]])  
}


## ---------------------------------------------------------------------------------------
##                      Step  1.5 Train and check the fitted discrepancies for multiple PCs 
## ---------------------------------------------------------------------------------------


# generate priors for discrepancy term: 
da <- darg(list(mle=TRUE, max=5, min=eps), X, samp.size = nrow(X)); da  # also have stronger regularization on bias term

## prior for nugget effects: need to be careful in principal component space 
## ga <- garg(list(mle=TRUE, max=1), Bias_pr$x[,1]); ga
ga_1 <- garg(list(mle=TRUE, max=1), Bias_pr_1$x[,1]); ga_1
ga_2 <- garg(list(mle=TRUE, max=1), Bias_pr_2$x[,1]); ga_2
ga_3 <- garg(list(mle=TRUE, max=1), Bias_pr_3$x[,1]); ga_3
ga_4 <- garg(list(mle=TRUE, max=1), Bias_pr_4$x[,1]); ga_4

## load one good value of u estimated from optimization: 
u <- c(0.5667774, 0.8118366, 0.5064965, 0.4891309)

# fit for the separate bias term in PC space for each outputs:
yM_pc_1 <- rep(NA, 292)
yM_pc_2 <- rep(NA, 292)
yM_pc_3 <- rep(NA, 292)
yM_pc_4 <- rep(NA, 292)

# loop over 292 Xf:
for (i in 1:292){
  CX_pc_1 <- covar.sep(u, U[list_na_index[[i]], ], d = MLEs_Kdir_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
  yM_pc_1[i] <- CX_pc_1 %*% CY_pc_1[[i]]  # predict from ym for means  
  
  CX_pc_2 <- covar.sep(u, U[list_na_index[[i]], ], d = MLEs_Kcross_pc[i, 1:4], g = MLEs_Kcross_pc[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
  yM_pc_2[i] <- CX_pc_2 %*% CY_pc_2[[i]]  # predict from ym for means  
  
  CX_pc_3 <- covar.sep(u, U[list_na_index[[i]], ], d = MLEs_cdir_pc[i, 1:4], g = MLEs_cdir_pc[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
  yM_pc_3[i] <- CX_pc_3 %*% CY_pc_3[[i]]  # predict from ym for means  
  
  CX_pc_4 <- covar.sep(u, U[list_na_index[[i]], ], d = MLEs_ccross_pc[i, 1:4], g = MLEs_ccross_pc[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
  yM_pc_4[i] <- CX_pc_4 %*% CY_pc_4[[i]]  # predict from ym for means  
  
}


## rotate in the PC spaces: 
colnames(YF_1) <- c("y2", "y5", "y9", "y11")
colnames(YF_2) <- c("y2", "y5", "y9", "y11")
colnames(YF_3) <- c("y2", "y5", "y9", "y11")
colnames(YF_4) <- c("y2", "y5", "y9", "y11")
yF_pc_1 <- predict(Bias_pr_1, YF_1)[, 1]
yF_pc_2 <- predict(Bias_pr_2, YF_2)[, 1]
yF_pc_3 <- predict(Bias_pr_3, YF_3)[, 1]
yF_pc_4 <- predict(Bias_pr_4, YF_4)[, 1]

## use a fixed global PCA on every new U values for fair comparison: 
bias_u_1 <-  yF_pc_1 - yM_pc_1    ##  observed bias given U value
bias_u_2 <-  yF_pc_2 - yM_pc_2    ##  observed bias given U value
bias_u_3 <-  yF_pc_3 - yM_pc_3    ##  observed bias given U value
bias_u_4 <-  yF_pc_4 - yM_pc_4    ##  observed bias given U value


## Fit GPs on the first principal component for all outputs:
bhat_1 <- newGPsep(X, bias_u_1 , d=rep(da$start, ncol(X)), g=ga_1$start, dK=TRUE)
bhat_2 <- newGPsep(X, bias_u_2 , d=rep(da$start, ncol(X)), g=ga_2$start, dK=TRUE)
bhat_3 <- newGPsep(X, bias_u_3 , d=rep(da$start, ncol(X)), g=ga_3$start, dK=TRUE)
bhat_4 <- newGPsep(X, bias_u_4 , d=rep(da$start, ncol(X)), g=ga_4$start, dK=TRUE)

mle_b_1 <- mleGPsep(bhat_1, param="both", tmin=c(da$min, ga_1$min), tmax=c(da$max, ga_1$max), ab=c(da$ab, ga_1$ab), maxit=1000)
mle_b_2 <- mleGPsep(bhat_2, param="both", tmin=c(da$min, ga_2$min), tmax=c(da$max, ga_2$max), ab=c(da$ab, ga_2$ab), maxit=1000)
mle_b_3 <- mleGPsep(bhat_3, param="both", tmin=c(da$min, ga_3$min), tmax=c(da$max, ga_3$max), ab=c(da$ab, ga_3$ab), maxit=1000)
mle_b_4 <- mleGPsep(bhat_4, param="both", tmin=c(da$min, ga_4$min), tmax=c(da$max, ga_4$max), ab=c(da$ab, ga_4$ab), maxit=1000)

# regularize with a marginal beta(2, 2) prior on U:
Uprior <-  - sum(dbeta(u, 2, 2, log = TRUE))  # negative log likelihood to minimize: be careful with being positive in MCMC

## Check negative log-likelihood and priors: 
- llikGPsep(bhat_1, dab=da$ab, gab=ga_1$ab) +  Uprior
- llikGPsep(bhat_2, dab=da$ab, gab=ga_2$ab) +  Uprior
- llikGPsep(bhat_3, dab=da$ab, gab=ga_3$ab) +  Uprior
- llikGPsep(bhat_4, dab=da$ab, gab=ga_4$ab) +  Uprior


# estimate the tau2bhat for all the bias terms: 
Cd_1 <- covar.sep(X, X, d=mle_b_1$theta[1:13], g=mle_b_1$theta[14])  # this part need nugget g
Cd_inv_1 <- solve(Cd_1)
taub2hat_1 <- drop(t(bias_u_1) %*% Cd_inv_1  %*% (bias_u_1 ) / length(bias_u_1)); taub2hat_1
Vd_1 <- taub2hat_1 * Cd_1 # the covariance of the bias term keeps the same with new proposed theta

Cd_2 <- covar.sep(X, X, d=mle_b_2$theta[1:13], g=mle_b_2$theta[14])  # this part need nugget g
Cd_inv_2 <- solve(Cd_2)
taub2hat_2 <- drop(t(bias_u_2) %*% Cd_inv_2  %*% (bias_u_2 ) / length(bias_u_2)); taub2hat_2
Vd_2 <- taub2hat_2 * Cd_2 # the covariance of the bias term keeps the same with new proposed theta

Cd_3 <- covar.sep(X, X, d=mle_b_3$theta[1:13], g=mle_b_3$theta[14])  # this part need nugget g
Cd_inv_3 <- solve(Cd_3)
taub2hat_3 <- drop(t(bias_u_3) %*% Cd_inv_3  %*% (bias_u_3 ) / length(bias_u_3)); taub2hat_3
Vd_3 <- taub2hat_3 * Cd_3 # the covariance of the bias term keeps the same with new proposed theta

Cd_4 <- covar.sep(X, X, d=mle_b_4$theta[1:13], g=mle_b_4$theta[14])  # this part need nugget g
Cd_inv_4 <- solve(Cd_4)
taub2hat_4 <- drop(t(bias_u_4) %*% Cd_inv_4  %*% (bias_u_4 ) / length(bias_u_4)); taub2hat_4
Vd_4 <- taub2hat_4 * Cd_4 # the covariance of the bias term keeps the same with new proposed theta


bpred_pc_1 <- predGPsep(bhat_1, X, lite=FALSE)  ## in the PCA space
bpred_pc_2 <- predGPsep(bhat_2, X, lite=FALSE)  ## in the PCA space
bpred_pc_3 <- predGPsep(bhat_3, X, lite=FALSE)  ## in the PCA space
bpred_pc_4 <- predGPsep(bhat_4, X, lite=FALSE)  ## in the PCA space


## Kdir: 
## check the emulation for observed bias in first PC space:    
b_pc1_pred = data.frame(
  mean   =   bpred_pc_1$mean,
  var = bpred_pc_1$Sigma,
  yF = bias_u_1)
b_pc1_pred$se <- sqrt(diag(bpred_pc_1$Sigma))

par(mfrow = c(1, 1))
ggplot(b_pc1_pred,               
       aes(x = yF,
           y = mean)) +
  geom_errorbar(aes(ymin = mean - 1.96 * se,
                    ymax = mean + 1.96 * se),
                width=0.1) +
  geom_point(shape = 19,
             size  = 2, 
             color ="red") +
  geom_abline(col = "blue") + 
  ggtitle("Prediction on bias in 1st PC space (95% Credible Interval)") +
  labs( x="Observed bias", y = "Predicted bias") + 
  theme(plot.title = element_text(size = 50),  
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20))  +
  theme_bw() 

## Kcross: 
## check the emulation for observed bias in first PC space:    
b_pc1_pred = data.frame(
  mean   =   bpred_pc_2$mean,
  var = bpred_pc_2$Sigma,
  yF = bias_u_2)
b_pc1_pred$se <- sqrt(diag(bpred_pc_2$Sigma))

par(mfrow = c(1, 1))
ggplot(b_pc1_pred,               
       aes(x = yF,
           y = mean)) +
  geom_errorbar(aes(ymin = mean - 1.96 * se,
                    ymax = mean + 1.96 * se),
                width=0.1) +
  geom_point(shape = 19,
             size  = 2, 
             color ="red") +
  geom_abline(col = "blue") + 
  ggtitle("Prediction on bias in 1st PC space (95% Credible Interval)") +
  labs( x="Observed bias", y = "Predicted bias") + 
  theme(plot.title = element_text(size = 50),  
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20))  +
  theme_bw() 


## Cdir:  
## check the emulation for observed bias in first PC space:    
b_pc1_pred = data.frame(
  mean   =   bpred_pc_3$mean,
  var = bpred_pc_3$Sigma,
  yF = bias_u_3)
b_pc1_pred$se <- sqrt(diag(bpred_pc_3$Sigma))

par(mfrow = c(1, 1))
ggplot(b_pc1_pred,               
       aes(x = yF,
           y = mean)) +
  geom_errorbar(aes(ymin = mean - 1.96 * se,
                    ymax = mean + 1.96 * se),
                width=0.1) +
  geom_point(shape = 19,
             size  = 2, 
             color ="red") +
  geom_abline(col = "blue") + 
  ggtitle("Prediction on bias in 1st PC space (95% Credible Interval)") +
  labs( x="Observed bias", y = "Predicted bias") + 
  theme(plot.title = element_text(size = 50),  
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20))  +
  theme_bw() 

## Ccross:  
## check the emulation for observed bias in first PC space:    
b_pc1_pred = data.frame(
  mean   =   bpred_pc_4$mean,
  var = bpred_pc_4$Sigma,
  yF = bias_u_4)
b_pc1_pred$se <- sqrt(diag(bpred_pc_4$Sigma))

par(mfrow = c(1, 1))
ggplot(b_pc1_pred,               
       aes(x = yF,
           y = mean)) +
  geom_errorbar(aes(ymin = mean - 1.96 * se,
                    ymax = mean + 1.96 * se),
                width=0.1) +
  geom_point(shape = 19,
             size  = 2, 
             color ="red") +
  geom_abline(col = "blue") + 
  ggtitle("Prediction on bias in 1st PC space (95% Credible Interval)") +
  labs( x="Observed bias", y = "Predicted bias") + 
  theme(plot.title = element_text(size = 50),  
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20))  +
  theme_bw() 


### Pre-calculate the 1st PC for all 4 yms: rotate all at once and use directly in MCMC

colnames(YM_1) <- c("y2", "y5", "y9", "y11")  ## name it for PCA rotation
ym_pc_1  <- predict(Bias_pr_1, YM_1)
ym_pc1 <- ym_pc_1[, 1]

colnames(YM_2) <- c("y2", "y5", "y9", "y11")  ## name it for PCA rotation
ym_pc_2  <- predict(Bias_pr_2, YM_2)
ym_pc2 <- ym_pc_2[, 1]

colnames(YM_3) <- c("y2", "y5", "y9", "y11")  ## name it for PCA rotation
ym_pc_3  <- predict(Bias_pr_3, YM_3)
ym_pc3 <- ym_pc_3[, 1]

colnames(YM_4) <- c("y2", "y5", "y9", "y11")  ## name it for PCA rotation
ym_pc_4  <- predict(Bias_pr_4, YM_4)
ym_pc4 <- ym_pc_4[, 1]


## Pre-calculate the 1st PC for all 4 yfs: rotate all at once and use directly 

colnames(YF_1) <- c("y2", "y5", "y9", "y11")
yf_pc_1  <- predict(Bias_pr_1, YF_1)
yf_pc1_1 <- yf_pc_1[, 1]

colnames(YF_2) <- c("y2", "y5", "y9", "y11")
yf_pc_2  <- predict(Bias_pr_2, YF_2)
yf_pc1_2 <- yf_pc_2[, 1]

colnames(YF_3) <- c("y2", "y5", "y9", "y11")
yf_pc_3  <- predict(Bias_pr_3, YF_3)
yf_pc1_3 <- yf_pc_3[, 1]

colnames(YF_4) <- c("y2", "y5", "y9", "y11")
yf_pc_4  <- predict(Bias_pr_4, YF_4)
yf_pc1_4 <- yf_pc_4[, 1]


## -------------------------------------------------------------------------------
##                      Step 2. Calculate the posterior: initial evaluation
## -------------------------------------------------------------------------------

## start with the initial posterior: for U0 
u0 <- u

## ------------------------------------------------------------------------
##                      Step 2.1. inverse part: 
## ------------------------------------------------------------------------

## Vo_inv_(Vob)^T: updated because Vob is a function of u
## get Vob^T first: this is changing with u at each iteration
Vob_pc_1 <- list()  ## combined in 1st PC space
Vob_pc_2 <- list()  ## combined in 1st PC space
Vob_pc_3 <- list()  ## combined in 1st PC space
Vob_pc_4 <- list()  ## combined in 1st PC space

for (i in 1:292) {   
  Cob_u_i_pc_1 <-  covar.sep(u0, U[list_na_index[[i]], ],  d = MLEs_Kdir_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) # 1 by 1000
  Cob_u_i_pc_2 <-  covar.sep(u0, U[list_na_index[[i]], ],  d = MLEs_Kcross_pc[i, 1:4], g = MLEs_Kcross_pc[i, 5]) # 1 by 1000
  Cob_u_i_pc_3 <-  covar.sep(u0, U[list_na_index[[i]], ],  d = MLEs_cdir_pc[i, 1:4], g = MLEs_cdir_pc[i, 5]) # 1 by 1000
  Cob_u_i_pc_4 <-  covar.sep(u0, U[list_na_index[[i]], ],  d = MLEs_ccross_pc[i, 1:4], g = MLEs_ccross_pc[i, 5]) # 1 by 1000
  
  ## Combined: 
  Vob_pc_1[[i]] <-  t(Tau2hat_1[i, 1]  * Cob_u_i_pc_1)
  Vob_pc_2[[i]] <-  t(Tau2hat_2[i, 1]  * Cob_u_i_pc_2)
  Vob_pc_3[[i]] <-  t(Tau2hat_3[i, 1]  * Cob_u_i_pc_3)
  Vob_pc_4[[i]] <-  t(Tau2hat_4[i, 1]  * Cob_u_i_pc_4)
}

Vo_inv_Vob_t_pc_1 <- list()
Vo_inv_Vob_t_pc_2 <- list()
Vo_inv_Vob_t_pc_3 <- list()
Vo_inv_Vob_t_pc_4 <- list()

for (i in 1:292) {
  Vo_inv_Vob_t_pc_1[[i]] <- Vi_inv_pc_1[[i]] %*% Vob_pc_1[[i]]   # each part is 1000 by 1
  Vo_inv_Vob_t_pc_2[[i]] <- Vi_inv_pc_2[[i]] %*% Vob_pc_2[[i]]   # each part is 1000 by 1
  Vo_inv_Vob_t_pc_3[[i]] <- Vi_inv_pc_3[[i]] %*% Vob_pc_3[[i]]   # each part is 1000 by 1
  Vo_inv_Vob_t_pc_4[[i]] <- Vi_inv_pc_4[[i]] %*% Vob_pc_4[[i]]   # each part is 1000 by 1
}

## C_in_u: also updated with u

##  Vob_Vo_inv_Vbo: 292 by 292 diagonal matrix
Vob_Vo_inv_Vbo_pc_1 <- matrix(0, nrow=292, ncol=292)
Vob_Vo_inv_Vbo_pc_2 <- matrix(0, nrow=292, ncol=292)
Vob_Vo_inv_Vbo_pc_3 <- matrix(0, nrow=292, ncol=292)
Vob_Vo_inv_Vbo_pc_4 <- matrix(0, nrow=292, ncol=292)

for (i in 1:292) {
  Vob_Vo_inv_Vbo_pc_1[i, i] <- t(Vob_pc_1[[i]]) %*% Vi_inv_pc_1[[i]] %*% Vob_pc_1[[i]]   # each part is a scalar
  Vob_Vo_inv_Vbo_pc_2[i, i] <- t(Vob_pc_2[[i]]) %*% Vi_inv_pc_2[[i]] %*% Vob_pc_2[[i]]   # each part is a scalar
  Vob_Vo_inv_Vbo_pc_3[i, i] <- t(Vob_pc_3[[i]]) %*% Vi_inv_pc_3[[i]] %*% Vob_pc_3[[i]]   # each part is a scalar
  Vob_Vo_inv_Vbo_pc_4[i, i] <- t(Vob_pc_4[[i]]) %*% Vi_inv_pc_4[[i]] %*% Vob_pc_4[[i]]   # each part is a scalar
}

## Vb remains the same with different u values: 
vI_nf_1 <- diag(0, 292)
vI_nf_2 <- diag(0, 292)
vI_nf_3 <- diag(0, 292)
vI_nf_4 <- diag(0, 292)

for (i in 1:292){  
  vI_nf_1[i, i]  <- covar.sep(u,  t(u), d = MLEs_Kdir_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) * Tau2hat_1[i, 1]    # this part has nugget g
  vI_nf_2[i, i]  <- covar.sep(u,  t(u), d = MLEs_Kcross_pc[i, 1:4], g = MLEs_Kcross_pc[i, 5]) * Tau2hat_2[i, 1]    # this part has nugget g
  vI_nf_3[i, i]  <- covar.sep(u,  t(u), d = MLEs_cdir_pc[i, 1:4], g = MLEs_cdir_pc[i, 5]) * Tau2hat_3[i, 1]    # this part has nugget g
  vI_nf_4[i, i]  <- covar.sep(u,  t(u), d = MLEs_ccross_pc[i, 1:4], g = MLEs_ccross_pc[i, 5]) * Tau2hat_4[i, 1]    # this part has nugget g
}

Vb_pc_1 <-  vI_nf_1   +  Vd_1
Vb_pc_2 <-  vI_nf_2   +  Vd_2
Vb_pc_3 <-  vI_nf_3   +  Vd_3
Vb_pc_4 <-  vI_nf_4   +  Vd_4

C_u_1 <-  Vb_pc_1 - Vob_Vo_inv_Vbo_pc_1
C_u_2 <-  Vb_pc_2 - Vob_Vo_inv_Vbo_pc_2
C_u_3 <-  Vb_pc_3 - Vob_Vo_inv_Vbo_pc_3
C_u_4 <-  Vb_pc_4 - Vob_Vo_inv_Vbo_pc_4

C_u_inv_1 <- solve(C_u_1)
C_u_inv_2 <- solve(C_u_2)
C_u_inv_3 <- solve(C_u_3)
C_u_inv_4 <- solve(C_u_4)

c1_1 <- t(yf_pc1_1) %*% C_u_inv_1 %*% (yf_pc1_1) ; c1_1
c1_2 <- t(yf_pc1_2) %*% C_u_inv_2 %*% (yf_pc1_2) ; c1_2
c1_3 <- t(yf_pc1_3) %*% C_u_inv_3 %*% (yf_pc1_3) ; c1_3
c1_4 <- t(yf_pc1_4) %*% C_u_inv_4 %*% (yf_pc1_4) ; c1_4

### fixed part: c2 = t(y - m) A11_inv (y-m), only need to be calculated once. 
c2_1 <- rep(NA, 292)
c2_2 <- rep(NA, 292)
c2_3 <- rep(NA, 292)
c2_4 <- rep(NA, 292)

for (i in 1:292) {
  yi_pc1_1 <- ym_pc1[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_pc1_2 <- ym_pc2[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_pc1_3 <- ym_pc3[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_pc1_4 <- ym_pc4[1:1000 + 1000 * (i-1)][list_na_index[[i]]]

  c2_1[i] <- t( yi_pc1_1 ) %*% Vi_inv_pc_1[[i]]  %*%  yi_pc1_1 
  c2_2[i] <- t( yi_pc1_2 ) %*% Vi_inv_pc_2[[i]]  %*%  yi_pc1_2 
  c2_3[i] <- t( yi_pc1_3 ) %*% Vi_inv_pc_3[[i]]  %*%  yi_pc1_3 
  c2_4[i] <- t( yi_pc1_4 ) %*% Vi_inv_pc_4[[i]]  %*%  yi_pc1_4 
}

c2_1 <- sum(c2_1); c2_1 
c2_2 <- sum(c2_2); c2_2 
c2_3 <- sum(c2_3); c2_3 
c2_4 <- sum(c2_4); c2_4 

### updating parts: 
c3_1 <- rep(NA, 292)
c3_2 <- rep(NA, 292)
c3_3 <- rep(NA, 292)
c3_4 <- rep(NA, 292)

c6_1 <- rep(NA, 292)
c6_2 <- rep(NA, 292)
c6_3 <- rep(NA, 292)
c6_4 <- rep(NA, 292)

for (i in 1:292) {
  ## calculate t(y) A11_inv A_12 in lists
  yi_pc1_1 <- ym_pc1[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_pc1_2 <- ym_pc2[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_pc1_3 <- ym_pc3[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_pc1_4 <- ym_pc4[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  
  c3_1[i] <- t(yi_pc1_1) %*% Vo_inv_Vob_t_pc_1[[i]] 
  c3_2[i] <- t(yi_pc1_2) %*% Vo_inv_Vob_t_pc_2[[i]] 
  c3_3[i] <- t(yi_pc1_3) %*% Vo_inv_Vob_t_pc_3[[i]] 
  c3_4[i] <- t(yi_pc1_4) %*% Vo_inv_Vob_t_pc_4[[i]] 
  
  c6_1[i] <- t(Vo_inv_Vob_t_pc_1[[i]]) %*% yi_pc1_1
  c6_2[i] <- t(Vo_inv_Vob_t_pc_2[[i]]) %*% yi_pc1_2
  c6_3[i] <- t(Vo_inv_Vob_t_pc_3[[i]]) %*% yi_pc1_3
  c6_4[i] <- t(Vo_inv_Vob_t_pc_4[[i]]) %*% yi_pc1_4
  
}

c4_1 <- t(c3_1) %*% C_u_inv_1 %*% c3_1  ## the whole updating part for the 1st quadratic form
c4_2 <- t(c3_2) %*% C_u_inv_2 %*% c3_2  ## the whole updating part for the 1st quadratic form
c4_3 <- t(c3_3) %*% C_u_inv_3 %*% c3_3  ## the whole updating part for the 1st quadratic form
c4_4 <- t(c3_4) %*% C_u_inv_4 %*% c3_4  ## the whole updating part for the 1st quadratic form

c7_1 <- rep(NA, 292) ## 
c7_2 <- rep(NA, 292) ## 
c7_3 <- rep(NA, 292) ## 
c7_4 <- rep(NA, 292) ## 

for (i in 1:292) {
  c7_1[i] <- t(C_u_inv_1[, i]) %*% c6_1
  c7_2[i] <- t(C_u_inv_2[, i]) %*% c6_2
  c7_3[i] <- t(C_u_inv_3[, i]) %*% c6_3
  c7_4[i] <- t(C_u_inv_4[, i]) %*% c6_4
}

## c5 is the whole t(z) V21 y mixed z y quadratic part: 
c5_1 <- - t(yf_pc1_1) %*% c7_1
c5_2 <- - t(yf_pc1_2) %*% c7_2
c5_3 <- - t(yf_pc1_3) %*% c7_3
c5_4 <- - t(yf_pc1_4) %*% c7_4


## ------------------------------------------------------------------------
##                      Step 2.2. determinant part: 
## ------------------------------------------------------------------------

## fixed part: log(det(A11))
c8_1 <- rep(NA, 292)
c8_2 <- rep(NA, 292)
c8_3 <- rep(NA, 292)
c8_4 <- rep(NA, 292)

for (i in 1:292){
  c8_1[i] <- determinant(Vi_pc_1[[i]], logarithm=TRUE)$modulus
  c8_2[i] <- determinant(Vi_pc_2[[i]], logarithm=TRUE)$modulus
  c8_3[i] <- determinant(Vi_pc_3[[i]], logarithm=TRUE)$modulus
  c8_4[i] <- determinant(Vi_pc_4[[i]], logarithm=TRUE)$modulus
}

c9_1 <- sum(c8_1); c9_1
c9_2 <- sum(c8_2); c9_2
c9_3 <- sum(c8_3); c9_3
c9_4 <- sum(c8_4); c9_4

## updated part: log(det(C_u))
c10_1 <- determinant(C_u_1, logarithm=TRUE)$modulus; c10_1
c10_2 <- determinant(C_u_2, logarithm=TRUE)$modulus; c10_2
c10_3 <- determinant(C_u_3, logarithm=TRUE)$modulus; c10_3
c10_4 <- determinant(C_u_4, logarithm=TRUE)$modulus; c10_4

## ------------------------------------------------------------------------
##                      Step 2.3. log posterior: put all parts together
## ------------------------------------------------------------------------

# log likelihood: 
log_like <- (-1/2 *(c9_1 + c10_1) - (c1_1 + c2_1 + c4_1 + 2* c5_1) / 2 
             -1/2 *(c9_2 + c10_2) - (c1_2 + c2_2 + c4_2 + 2* c5_2) / 2 
             -1/2 *(c9_3 + c10_3) - (c1_3 + c2_3 + c4_3 + 2* c5_3) / 2 
             -1/2 *(c9_4 + c10_4) - (c1_4 + c2_4 + c4_4 + 2* c5_4) / 2 ) 
print(log_like)

# log prior on theta: 
logprior <-  sum(dbeta(u0, 2, 2, log = TRUE))  # beta(2, 2) prior

# log posterior
logpost_u <- logprior + log_like
print(logpost_u)



## ---------------------------------------------------------
##                      Step 3. MCMC Sampling 
## ---------------------------------------------------------

# number of total iterations
its <- 10   ## Try 10 iterations
## its <- 10^5 # actual number sampled in the paper

Us <- matrix(NA, nrow=its, ncol=4)
Us[1, ] <- u0

logposts <- rep(NA, its)
logposts[1] <- logpost_u

# use step sizes for different theta:
steps <- c(.1, .1, .1, .1) 

## Metropolis-within-Gibbs sampler: 
for (iter in 2:its){
  ut <- Us[iter - 1, ]
  
  ## proposal:  
  # Gibbs sampling for each dimension:
  u_cand <- ut
  remainder <- iter %% 4 +1
  u_cand[remainder] <- ut[remainder] + rnorm(1, 0, steps[remainder]) 
  
  # compute the log posterior: 
  ## u_cand enter the covariance matrix here: through V_ob only
  for (i in 1:292) {   
    Cob_u_i_pc_1 <-  covar.sep(u_cand, U[list_na_index[[i]], ],  d = MLEs_Kdir_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) # 1 by 1000
    Cob_u_i_pc_2 <-  covar.sep(u_cand, U[list_na_index[[i]], ],  d = MLEs_Kcross_pc[i, 1:4], g = MLEs_Kcross_pc[i, 5]) # 1 by 1000
    Cob_u_i_pc_3 <-  covar.sep(u_cand, U[list_na_index[[i]], ],  d = MLEs_cdir_pc[i, 1:4], g = MLEs_cdir_pc[i, 5]) # 1 by 1000
    Cob_u_i_pc_4 <-  covar.sep(u_cand, U[list_na_index[[i]], ],  d = MLEs_ccross_pc[i, 1:4], g = MLEs_ccross_pc[i, 5]) # 1 by 1000
    
    Vob_pc_1[[i]] <-  t(Tau2hat_1[i, 1]  * Cob_u_i_pc_1)
    Vob_pc_2[[i]] <-  t(Tau2hat_2[i, 1]  * Cob_u_i_pc_2)
    Vob_pc_3[[i]] <-  t(Tau2hat_3[i, 1]  * Cob_u_i_pc_3)
    Vob_pc_4[[i]] <-  t(Tau2hat_4[i, 1]  * Cob_u_i_pc_4)
  }
  
  Vo_inv_Vob_t_pc_1 <- list()
  Vo_inv_Vob_t_pc_2 <- list()
  Vo_inv_Vob_t_pc_3 <- list()
  Vo_inv_Vob_t_pc_4 <- list()
  
  for (i in 1:292) {
    Vo_inv_Vob_t_pc_1[[i]] <- Vi_inv_pc_1[[i]] %*% Vob_pc_1[[i]]   # each part is 1000 by 1
    Vo_inv_Vob_t_pc_2[[i]] <- Vi_inv_pc_2[[i]] %*% Vob_pc_2[[i]]   # each part is 1000 by 1
    Vo_inv_Vob_t_pc_3[[i]] <- Vi_inv_pc_3[[i]] %*% Vob_pc_3[[i]]   # each part is 1000 by 1
    Vo_inv_Vob_t_pc_4[[i]] <- Vi_inv_pc_4[[i]] %*% Vob_pc_4[[i]]   # each part is 1000 by 1
  }
  
  Vob_Vo_inv_Vbo_pc_1 <- matrix(0, nrow=292, ncol=292)
  Vob_Vo_inv_Vbo_pc_2 <- matrix(0, nrow=292, ncol=292)
  Vob_Vo_inv_Vbo_pc_3 <- matrix(0, nrow=292, ncol=292)
  Vob_Vo_inv_Vbo_pc_4 <- matrix(0, nrow=292, ncol=292)
  
  for (i in 1:292) {
    Vob_Vo_inv_Vbo_pc_1[i, i] <- t(Vob_pc_1[[i]]) %*% Vi_inv_pc_1[[i]] %*% Vob_pc_1[[i]]   # each part is a scalar
    Vob_Vo_inv_Vbo_pc_2[i, i] <- t(Vob_pc_2[[i]]) %*% Vi_inv_pc_2[[i]] %*% Vob_pc_2[[i]]   # each part is a scalar
    Vob_Vo_inv_Vbo_pc_3[i, i] <- t(Vob_pc_3[[i]]) %*% Vi_inv_pc_3[[i]] %*% Vob_pc_3[[i]]   # each part is a scalar
    Vob_Vo_inv_Vbo_pc_4[i, i] <- t(Vob_pc_4[[i]]) %*% Vi_inv_pc_4[[i]] %*% Vob_pc_4[[i]]   # each part is a scalar
  }
  
  C_u_1 <-  Vb_pc_1 - Vob_Vo_inv_Vbo_pc_1
  C_u_2 <-  Vb_pc_2 - Vob_Vo_inv_Vbo_pc_2
  C_u_3 <-  Vb_pc_3 - Vob_Vo_inv_Vbo_pc_3
  C_u_4 <-  Vb_pc_4 - Vob_Vo_inv_Vbo_pc_4
  
  C_u_inv_1 <- solve(C_u_1)
  C_u_inv_2 <- solve(C_u_2)
  C_u_inv_3 <- solve(C_u_3)
  C_u_inv_4 <- solve(C_u_4)
  
  ## The quadratic form for yM and yF: 
  c1_1 <- t(yf_pc1_1) %*% C_u_inv_1 %*% (yf_pc1_1) 
  c1_2 <- t(yf_pc1_2) %*% C_u_inv_2 %*% (yf_pc1_2) 
  c1_3 <- t(yf_pc1_3) %*% C_u_inv_3 %*% (yf_pc1_3) 
  c1_4 <- t(yf_pc1_4) %*% C_u_inv_4 %*% (yf_pc1_4) 

  ### updating parts: 
  c3_1 <- rep(NA, 292)
  c3_2 <- rep(NA, 292)
  c3_3 <- rep(NA, 292)
  c3_4 <- rep(NA, 292)
  
  c6_1 <- rep(NA, 292)
  c6_2 <- rep(NA, 292)
  c6_3 <- rep(NA, 292)
  c6_4 <- rep(NA, 292)
  
  for (i in 1:292) {
    ## calculate t(y) A11_inv A_12 by list
    yi_pc1_1 <- ym_pc1[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
    yi_pc1_2 <- ym_pc2[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
    yi_pc1_3 <- ym_pc3[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
    yi_pc1_4 <- ym_pc4[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
    
    c3_1[i] <- t(yi_pc1_1) %*% Vo_inv_Vob_t_pc_1[[i]] 
    c3_2[i] <- t(yi_pc1_2) %*% Vo_inv_Vob_t_pc_2[[i]] 
    c3_3[i] <- t(yi_pc1_3) %*% Vo_inv_Vob_t_pc_3[[i]] 
    c3_4[i] <- t(yi_pc1_4) %*% Vo_inv_Vob_t_pc_4[[i]] 
    
    ## t(z)V21y = c5 part: 
    c6_1[i] <- t(Vo_inv_Vob_t_pc_1[[i]]) %*% yi_pc1_1
    c6_2[i] <- t(Vo_inv_Vob_t_pc_2[[i]]) %*% yi_pc1_2
    c6_3[i] <- t(Vo_inv_Vob_t_pc_3[[i]]) %*% yi_pc1_3
    c6_4[i] <- t(Vo_inv_Vob_t_pc_4[[i]]) %*% yi_pc1_4
  }

  c4_1 <- t(c3_1) %*% C_u_inv_1 %*% c3_1  ## the whole updating part for the 1st quadratic form
  c4_2 <- t(c3_2) %*% C_u_inv_2 %*% c3_2  ## the whole updating part for the 1st quadratic form
  c4_3 <- t(c3_3) %*% C_u_inv_3 %*% c3_3  ## the whole updating part for the 1st quadratic form
  c4_4 <- t(c3_4) %*% C_u_inv_4 %*% c3_4  ## the whole updating part for the 1st quadratic form

  c7_1 <- rep(NA, 292) ## 
  c7_2 <- rep(NA, 292) ## 
  c7_3 <- rep(NA, 292) ## 
  c7_4 <- rep(NA, 292) ## 
  
  for (i in 1:292) {
    c7_1[i] <- t(C_u_inv_1[, i]) %*% c6_1
    c7_2[i] <- t(C_u_inv_2[, i]) %*% c6_2
    c7_3[i] <- t(C_u_inv_3[, i]) %*% c6_3
    c7_4[i] <- t(C_u_inv_4[, i]) %*% c6_4
  }
  
  c5_1 <- - t(yf_pc1_1) %*% c7_1
  c5_2 <- - t(yf_pc1_2) %*% c7_2
  c5_3 <- - t(yf_pc1_3) %*% c7_3
  c5_4 <- - t(yf_pc1_4) %*% c7_4
  
  ## updated part: log(det(C_u))
  c10_1 <- determinant(C_u_1, logarithm=TRUE)$modulus; c10_1
  c10_2 <- determinant(C_u_2, logarithm=TRUE)$modulus; c10_2
  c10_3 <- determinant(C_u_3, logarithm=TRUE)$modulus; c10_3
  c10_4 <- determinant(C_u_4, logarithm=TRUE)$modulus; c10_4
  
  # log likelihood: 
  log_like <- (-1/2 *(c9_1 + c10_1) - (c1_1 + c2_1 + c4_1 + 2* c5_1) / 2 
               -1/2 *(c9_2 + c10_2) - (c1_2 + c2_2 + c4_2 + 2* c5_2) / 2 
               -1/2 *(c9_3 + c10_3) - (c1_3 + c2_3 + c4_3 + 2* c5_3) / 2 
               -1/2 *(c9_4 + c10_4) - (c1_4 + c2_4 + c4_4 + 2* c5_4) / 2 ) 
  
  # log prior on theta: 
  logprior <-  sum(dbeta(u_cand, 2, 2, log = TRUE))  # beta(2, 2) prior
  
  # log posterior
  logpost_u_cand <- logprior + log_like
  
  # accept (or not) proposal for theta:
  if (log(runif(1)) < logpost_u_cand - logposts[iter-1]) {
    Us[iter, ] <- u_cand
    logposts[iter] <- logpost_u_cand
  }
  else {
    Us[iter, ] <- ut
    logposts[iter] <- logposts[iter-1]
  }
  print(iter)
  print(Us[iter, ])
  print(logposts[iter])
}

## -----------------------------------------------------------------
##                      Step 4. Collect samples and save objects  
## -----------------------------------------------------------------

## Save all objects: 
save(list=c("MLEs_Kdir_pc", "MLEs_Kcross_pc", "MLEs_cdir_pc", "MLEs_ccross_pc", 
            "mle_b_1", "mle_b_2","mle_b_3","mle_b_4",
            "Us",  "logposts", 
            "Tau2hat_1", "Tau2hat_2", "Tau2hat_3", "Tau2hat_4",
            "tau2hat_pc_1", "tau2hat_pc_2", "tau2hat_pc_3", "tau2hat_pc_4", 
            "urange", "steps"), file=paste("full_Bayes.RData", sep=""))
