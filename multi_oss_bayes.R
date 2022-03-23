
########################################################################################################
## Fully Bayesian calibration using on-site surrogates on all 16 outputs together
## 16 -- PC --> 4 ----> 1
## Combining frequencies into PCs for each out put
## On-site surrogates built from "multi_oss_emu.R"
## Full Bayes calibration: using a zero-mean Gaussian process to model the discrepency
## Re-fit on-site surrogates in the new first PC space for MCMC
########################################################################################################


## Set up: 
library(laGP)        # for discrepancy term GP
library(plgp)        # for covariance matrix
library(nloptr)      # for optimization with optim
library(factoextra)  # for PCA
library(corrplot)    # for PCA visual 
library(ggpubr)      # grab PCA plots 
library(lhs)
library(corrplot)
# library(expm)        # for square root matrix (no need, also very slow)


## ------------------------------------------------------------------------
##                       1. Setup before posterior sampling for U
## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
##                       1.1 Read and organize data:
## ------------------------------------------------------------------------


xnames <- c("Speed", "Pin", "Pout", "Temp", "Swirl", "viscosity", "gamma", "zeta", "Cin", "Cout", "D", "L", "Cd")
unames <- c("ffactor_Ns", "ffactor_Ms","ffactor_Nr", "ffactor_Mr")
eps <- sqrt(.Machine$double.eps)

# Set up X and U: 
lower_limits <- c(0.0028, -.15, .0021, -.2956)
upper_limits <- c(0.1, -eps, .0746, -eps) # ff_mr and ff_ms are only physically meaningful being negative
urange <- rbind(lower_limits, upper_limits)
colnames(urange) <- unames

# field data
Dfield <- read.csv("../../data/Global_File_StaticValues_MirkoFULL_rev15.csv")

# field responses Cdir and Ccross need to be scaled
# scale field response to a 0 to 100 scale: 
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

summary(Dfield)


# load all isotseal data: 
# all inputs with 4 outputs at freq = 2: 
D_freq_2 <- read.csv("../../data/U_LHS_1000_2.csv") # 1000 maximinLHS design for all the same U at each Xfield locations

# all outputs with 4 outputs at freqs 5, 7, 9, 11, 13: 
D_freq_5_13 <- read.csv("../../data/Multi_outputs_5_to_13_OSS_U_maximinLHS_1000.cvs") # 1000 maximinLHS design for all the same U at each Xfield locations


# Clean and scale all outputs: a silly way 
D_freq_2 <- D_freq_2[, -1]    # remove the first index column
D_freq_5_13 <- D_freq_5_13[, -1]    # remove the first index column

## Combine all 6 frequencies: 
Disot <- cbind(D_freq_2, D_freq_5_13) 
summary(Disot)

rm(D_freq_2)
rm(D_freq_5_13)



# load the on-site surrogates: 
load("multi_emu.RData", verbose = TRUE)


## scale all outputs: 

# Kdir: 
Disot$Kdir   <- Disot$Kdir / 10 ^ 6  # scale response to field unit
Disot$Kdir5  <- Disot$Kdir5 / 10 ^ 6  # scale response to field unit
Disot$Kdir7  <- Disot$Kdir7 / 10 ^ 6  # scale response to field unit
Disot$Kdir9  <- Disot$Kdir9 / 10 ^ 6  # scale response to field unit
Disot$Kdir11 <- Disot$Kdir11 / 10 ^ 6  # scale response to field unit
Disot$Kdir13 <- Disot$Kdir13 / 10 ^ 6  # scale response to field unit

# Kcross: 
Disot$Kcross   <- Disot$Kcross / 10 ^ 6  # scale response to field unit
Disot$Kcross5  <- Disot$Kcross5 / 10 ^ 6  # scale response to field unit
Disot$Kcross7  <- Disot$Kcross7 / 10 ^ 6  # scale response to field unit
Disot$Kcross9  <- Disot$Kcross9 / 10 ^ 6  # scale response to field unit
Disot$Kcross11 <- Disot$Kcross11 / 10 ^ 6  # scale response to field unit
Disot$Kcross13 <- Disot$Kcross13 / 10 ^ 6  # scale response to field unit

# Cdir: 
Disot$Cdir   <- Disot$Cdir / 10 ^ 3  # scale response to field unit
Disot$Cdir5  <- Disot$Cdir5 / 10 ^ 3  # scale response to field unit
Disot$Cdir7  <- Disot$Cdir7 / 10 ^ 3  # scale response to field unit
Disot$Cdir9  <- Disot$Cdir9 / 10 ^ 3  # scale response to field unit
Disot$Cdir11 <- Disot$Cdir11 / 10 ^ 3  # scale response to field unit
Disot$Cdir13 <- Disot$Cdir13 / 10 ^ 3  # scale response to field unit

# Cdir: 
Disot$Ccross   <- Disot$Ccross / 10 ^ 3  # scale response to field unit
Disot$Ccross5  <- Disot$Ccross5 / 10 ^ 3  # scale response to field unit
Disot$Ccross7  <- Disot$Ccross7 / 10 ^ 3  # scale response to field unit
Disot$Ccross9  <- Disot$Ccross9 / 10 ^ 3  # scale response to field unit
Disot$Ccross11 <- Disot$Ccross11 / 10 ^ 3  # scale response to field unit
Disot$Ccross13 <- Disot$Ccross13 / 10 ^ 3  # scale response to field unit


## Check both YM and YF are in the same scales: 
summary(Disot)
boxplot(Disot[, c(19:22, 24:27, 29:32, 34:37, 39:42, 44:47)], main ="all simulated outputs ym")

summary(Dfield)
boxplot(Dfield[, 21:44], main ="all field outputs yf")


## for  "Kdir"

## simulation:
ym_1_2 <- Disot$Kdir
ym_1_5 <- Disot$Kdir5
ym_1_9 <- Disot$Kdir9
ym_1_11 <- Disot$Kdir11

## set up matrix of field outputs: 
YF_1 <- data.frame(yf_1_2 = Dfield$Kxx_2, yf_1_5 = Dfield$Kxx_5, yf_1_9 = Dfield$Kxx_9, yf_1_11 = Dfield$Kxx_11)

# set up oss mles: 
OSS_1 <- list(oss2 = MLEs$Kdir2, oss5 = MLEs$Kdir5, oss9 = MLEs$Kdir9, oss11 = MLEs$Kdir11)


## for "Kcross"

## simulation: 
ym_2_2 <- Disot$Kcross
ym_2_5 <- Disot$Kcross5
ym_2_9 <- Disot$Kcross9
ym_2_11 <- Disot$Kcross11

## set up matrix of field outputs: 
YF_2 <- data.frame(yf_2_2 = Dfield$k_xy_2, yf_2_5 = Dfield$k_xy_5, yf_2_9 = Dfield$k_xy_9, yf_2_11 = Dfield$k_xy_11)

# set up oss mles: 
OSS_2 <- list(oss2 = MLEs$Kcross2, oss5 = MLEs$Kcross5, oss9 = MLEs$Kcross9, oss11 = MLEs$Kcross11)


## for "Cdir"
ym_3_2 <- Disot$Cdir
ym_3_5 <- Disot$Cdir5
ym_3_9 <- Disot$Cdir9
ym_3_11 <- Disot$Cdir11

## set up matrix of field outputs: 
YF_3 <- data.frame(yf_3_2 = Dfield$Cxx_2, yf_3_5 = Dfield$Cxx_5, yf_3_9 = Dfield$Cxx_9, yf_3_11 = Dfield$Cxx_11)

# set up oss mles: 
OSS_3 <- list(oss2 = MLEs$Cdir2, oss5 = MLEs$Cdir5, oss9 = MLEs$Cdir9, oss11 = MLEs$Cdir11)


## for "Ccross"
ym_4_2 <- Disot$Ccross
ym_4_5 <- Disot$Ccross5
ym_4_9 <- Disot$Ccross9
ym_4_11 <- Disot$Ccross11

## set up matrix of field outputs: 
YF_4 <- data.frame(yf_4_2 = Dfield$c_xy_2, yf_4_5 = Dfield$c_xy_5, yf_4_9 = Dfield$c_xy_9, yf_4_11 = Dfield$c_xy_11)

# set up oss mles: 
OSS_4 <- list(oss2 = MLEs$Ccross2, oss5 = MLEs$Ccross5, oss9 = MLEs$Ccross9, oss11 = MLEs$Ccross11)



# create an index list for NAs in each set of 1000 run isotseal data:
list_na_index <- list()

for (i in 1:292){  
  list_na_index[[i]] <- which(is.na(ym_2_2[1:1000 + 1000 * (i-1)]) == FALSE)
}   


# Scale both inputs X and U: 
XUori <- Disot[, 1:17]
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
Uori <- Disot[1:1000, 14:17]
U <- matrix(NA, nrow = nrow(Uori), ncol = ncol(Uori))


# scale XU for isotseal: only 1000 runs need to be scaled 
for (j in 1:4) {
  U[ , j] <- (Uori[ , j] - limits[1, 13+j]) / (limits[2, 13+j] - limits[1, 13+j]) 
}

summary(U)
pairs(U)
cor(U)  # random maximin LHS design

# scale X for field data
X <- matrix(NA, nrow = nrow(Xori), ncol = ncol(Xori))
for (j in 1:ncol(Xori)) {
  X[ , j] <- (Xori[ , j] - limits[1, j]) / (limits[2, j] - limits[1, j]) 
}


## ---------------------------------------------------------------------------------------
##                       1.2 PC rotation on whole U space for all observed discrepancy
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



BIAS_1 <- YF_1_rep - YM_1
BIAS_2 <- YF_2_rep - YM_2
BIAS_3 <- YF_3_rep - YM_3
BIAS_4 <- YF_4_rep - YM_4


colnames(BIAS_1) <- c("y2", "y5", "y9", "y11")
colnames(BIAS_2) <- c("y2", "y5", "y9", "y11")
colnames(BIAS_3) <- c("y2", "y5", "y9", "y11")
colnames(BIAS_4) <- c("y2", "y5", "y9", "y11")


summary(BIAS_1)
summary(BIAS_2)
summary(BIAS_3)
summary(BIAS_4)


# Remove missing: 
BIAS_1 <- na.omit(BIAS_1)
BIAS_2 <- na.omit(BIAS_2)
BIAS_3 <- na.omit(BIAS_3)
BIAS_4 <- na.omit(BIAS_4)


corrplot.mixed(cor(BIAS_1), lower = "number", upper = "ellipse",  tl.cex=0.75)
corrplot.mixed(cor(BIAS_2), lower = "number", upper = "ellipse",  tl.cex=0.75)
corrplot.mixed(cor(BIAS_3), lower = "number", upper = "ellipse",  tl.cex=0.75)
corrplot.mixed(cor(BIAS_4), lower = "number", upper = "ellipse",  tl.cex=0.75)



# Perform PCA on the frequencies for all 292,000 observed discrepancy: 
Bias_pr_1 <- prcomp(BIAS_1, scale = TRUE)
Bias_pr_2 <- prcomp(BIAS_2, scale = TRUE)
Bias_pr_3 <- prcomp(BIAS_3, scale = TRUE)
Bias_pr_4 <- prcomp(BIAS_4, scale = TRUE)

summary(Bias_pr_1)
summary(Bias_pr_2)
summary(Bias_pr_3)
summary(Bias_pr_4)


fviz_eig(Bias_pr_1, main = paste("Scree plot of observed bias for kdir"))
fviz_eig(Bias_pr_2, main = paste("Scree plot of observed bias for kcross"))
fviz_eig(Bias_pr_3, main = paste("Scree plot of observed bias for cdir"))
fviz_eig(Bias_pr_4, main = paste("Scree plot of observed bias for ccross"))


## Check the center, scale, and rotation constants: 
c_pc <- Bias_pr_1$center; c_pc  ## center vector for all responses
c_pc <- Bias_pr_2$center; c_pc  ## center vector for all responses
c_pc <- Bias_pr_3$center; c_pc  ## center vector for all responses
c_pc <- Bias_pr_4$center; c_pc  ## center vector for all responses

s_pc <- Bias_pr_1$scale; s_pc   ## scale vector for all responses
s_pc <- Bias_pr_2$scale; s_pc   ## scale vector for all responses
s_pc <- Bias_pr_3$scale; s_pc   ## scale vector for all responses
s_pc <- Bias_pr_4$scale; s_pc   ## scale vector for all responses

## just check: 
w1 <- Bias_pr_1$rotation[, 1]; w1## rotation vector for the first principal component
mu_m <- sum( - c_pc*w1/s_pc ); mu_m    ## Grand mean for all y^M after C&S and rotation

w1 <- Bias_pr_2$rotation[, 1]; w1## rotation vector for the first principal component
mu_m <- sum( - c_pc*w1/s_pc ); mu_m    ## Grand mean for all y^M after C&S and rotation

w1 <- Bias_pr_3$rotation[, 1]; w1## rotation vector for the first principal component
mu_m <- sum( - c_pc*w1/s_pc ); mu_m    ## Grand mean for all y^M after C&S and rotation

w1 <- Bias_pr_4$rotation[, 1]; w1## rotation vector for the first principal component
mu_m <- sum( - c_pc*w1/s_pc ); mu_m    ## Grand mean for all y^M after C&S and rotation


## ---------------------------------------------------------------------------------------
##                       1.3 Build new on-site surrogates on the rotated new space: 
## ---------------------------------------------------------------------------------------

n <- 1000  

# For kdir: 
MLEs_Kdir_pc <- data.frame(matrix(NA, nrow=292, ncol= 6))
names(MLEs_Kdir_pc) <- c("theta1", "theta2", "theta3", "theta4", "g", "itera")

MLEs_Kcross_pc <- data.frame(matrix(NA, nrow=292, ncol= 6))
names(MLEs_Kcross_pc) <- c("theta1", "theta2", "theta3", "theta4", "g", "itera")


MLEs_cdir_pc <- data.frame(matrix(NA, nrow=292, ncol= 6))
names(MLEs_cdir_pc) <- c("theta1", "theta2", "theta3", "theta4", "g", "itera")


MLEs_ccross_pc <- data.frame(matrix(NA, nrow=292, ncol= 6))
names(MLEs_ccross_pc) <- c("theta1", "theta2", "theta3", "theta4", "g", "itera")


## data.frame to collect out-of-sample prediction performance:
## prederror_Kdir_pc  <- data.frame(matrix(NA, nrow=292, ncol=7 ))
## names(prederror_Kdir_pc) <- c("Min", "Q1", "Median", "Mean", "Q3", "Max", "RMSE")



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
  
  
  ## Check out of sample prediction: 
  
  #p_Kdir_pc <- predGPsep( gpsep_pc1 , as.matrix(Utest0))
  
  #prederror_Kdir_pc[i, ] <- c(summary(p_Kdir_pc$mean - ytest_Kdir_pc), sqrt(mean((p_Kdir_pc$mean - ytest_Kdir_pc)^2)))
  
  
  ## Clean up 
  deleteGPsep(gpsep_1_pc1)
  deleteGPsep(gpsep_2_pc1)
  deleteGPsep(gpsep_3_pc1)
  deleteGPsep(gpsep_4_pc1)
  
  print(i)
}


summary(MLEs_Kdir_pc)
summary(MLEs_Kcross_pc)
summary(MLEs_cdir_pc)
summary(MLEs_ccross_pc)


## ---------------------------------------------------------------------------------------
##                       1.4 Pre-calculation for on-site kernels and their inverses
## ---------------------------------------------------------------------------------------

## BIG pre-calculation for MCMC: this part only needs to be pre-calculated once

# each of 292 locations have a different length due to a different rate of missing

## create CY() for all frequencies:  no need for posterior/ only need for prediction 
## CY_pc <- list()  # create a list of pre-calculated values, most of them are 1000 by 1

CY_pc_1 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1
CY_pc_2 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1
CY_pc_3 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1
CY_pc_4 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1



# each of 292 locations have a different length due to a different rate of missing
## Vi_pc <- list()  # combined Vi

Vi_pc_1 <- list()  # combined Vi
Vi_pc_2 <- list()  # combined Vi
Vi_pc_3 <- list()  # combined Vi
Vi_pc_4 <- list()  # combined Vi


## Vi_inv_pc <- list()  # combined Vi_inv

Vi_inv_pc_1 <- list()  # combined Vi_inv
Vi_inv_pc_2 <- list()  # combined Vi_inv
Vi_inv_pc_3 <- list()  # combined Vi_inv
Vi_inv_pc_4 <- list()  # combined Vi_inv


## Tau2hat <- matrix(NA, nrow = 292, ncol = 1)

Tau2hat_1 <- matrix(NA, nrow = 292, ncol = 1)
Tau2hat_2 <- matrix(NA, nrow = 292, ncol = 1)
Tau2hat_3 <- matrix(NA, nrow = 292, ncol = 1)
Tau2hat_4 <- matrix(NA, nrow = 292, ncol = 1)



## Precalculate: 

for (i in 1:292) {
  
  ## set up the outputs: 
  ## yi_2 <- ym2[1:1000 + 1000 * (i-1)][list_na_index[[i]]] 
  ## yi_5 <- ym5[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  ## yi_9 <- ym9[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  ## yi_11 <- ym11[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  
  ## YM_i <- cbind(yi_2, yi_5, yi_9, yi_11)
  ## colnames(YM_i) <- c("y2", "y5", "y9", "y11")  ## name each column for PCA rotation
  ## ym_i_pc1 <- predict(Bias_pr, YM_i)[,1]         ## center, standardize, and rotate all yM into the new space
  
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
  
  
  ## Fixed covariance part: 
  ## Cm_pc  <- covar.sep(U[list_na_index[[i]], ],  d = MLEs_Kdir_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) # 1000 by 1000
  
  Cm_pc_1  <- covar.sep(U[list_na_index[[i]], ],  d = MLEs_Kdir_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) # 1000 by 1000
  Cm_pc_2  <- covar.sep(U[list_na_index[[i]], ],  d = MLEs_Kcross_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) # 1000 by 1000
  Cm_pc_3  <- covar.sep(U[list_na_index[[i]], ],  d = MLEs_cdir_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) # 1000 by 1000
  Cm_pc_4  <- covar.sep(U[list_na_index[[i]], ],  d = MLEs_ccross_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) # 1000 by 1000
  
  
  ## Inverse of on-site designs: 
  ## Cm_inv_pc <- solve(Cm_pc)
  Cm_inv_pc_1 <- solve(Cm_pc_1)
  Cm_inv_pc_2 <- solve(Cm_pc_2)
  Cm_inv_pc_3 <- solve(Cm_pc_3)
  Cm_inv_pc_4 <- solve(Cm_pc_4)
  
  
  ## The fixed whole part: solve(Sigma(U, U)) %*% ym: only need for bias term fitting
  ## CY_pc[[i]] <- Cm_inv_pc %*% ym_i_pc1      # 1000 by 1 for all 292 runs
  
  CY_pc_1[[i]] <- Cm_inv_pc_1 %*% ym_1_i_pc1      # 1000 by 1 for all 292 runs
  CY_pc_2[[i]] <- Cm_inv_pc_2 %*% ym_2_i_pc1      # 1000 by 1 for all 292 runs
  CY_pc_3[[i]] <- Cm_inv_pc_3 %*% ym_3_i_pc1      # 1000 by 1 for all 292 runs
  CY_pc_4[[i]] <- Cm_inv_pc_4 %*% ym_4_i_pc1      # 1000 by 1 for all 292 runs
  
  
  ## save tau2hat's in this step: 
  ## this is important!!! Need to estimate tau2hat in the transformed space: basically multiply by (wj1/sj)^2 
  
  ## tau2hat_pc <- drop(t(ym_i_pc1) %*%  Cm_inv_pc  %*% ym_i_pc1 / length(ym_i_pc1))   # mle of tau^2
  
  tau2hat_pc_1 <- drop(t(ym_1_i_pc1) %*%  Cm_inv_pc_1  %*% ym_1_i_pc1 / length(ym_1_i_pc1))   # mle of tau^2
  tau2hat_pc_2 <- drop(t(ym_2_i_pc1) %*%  Cm_inv_pc_2  %*% ym_2_i_pc1 / length(ym_2_i_pc1))   # mle of tau^2
  tau2hat_pc_3 <- drop(t(ym_3_i_pc1) %*%  Cm_inv_pc_3  %*% ym_3_i_pc1 / length(ym_3_i_pc1))   # mle of tau^2
  tau2hat_pc_4 <- drop(t(ym_4_i_pc1) %*%  Cm_inv_pc_4  %*% ym_4_i_pc1 / length(ym_4_i_pc1))   # mle of tau^2

  ## Tau2hat[i, ] <- tau2hat_pc   # save each local tau 
  
  Tau2hat_1[i, ] <- tau2hat_pc_1   # save each local tau 
  Tau2hat_2[i, ] <- tau2hat_pc_2   # save each local tau 
  Tau2hat_3[i, ] <- tau2hat_pc_3   # save each local tau 
  Tau2hat_4[i, ] <- tau2hat_pc_4   # save each local tau 
    
  
  ## All Vis
  ## Vi_pc[[i]] <-  tau2hat_pc * Cm_pc              # save Vis for V(d), each is around 1000 by 1000
  Vi_pc_1[[i]] <-  tau2hat_pc_1 * Cm_pc_1              # save Vis for V(d), each is around 1000 by 1000
  Vi_pc_2[[i]] <-  tau2hat_pc_2 * Cm_pc_2              # save Vis for V(d), each is around 1000 by 1000
  Vi_pc_3[[i]] <-  tau2hat_pc_3 * Cm_pc_3              # save Vis for V(d), each is around 1000 by 1000
  Vi_pc_4[[i]] <-  tau2hat_pc_4 * Cm_pc_4              # save Vis for V(d), each is around 1000 by 1000
  
  
  ## All Vi_inv 
  ## Vi_inv_pc[[i]] <- solve(Vi_pc[[i]])  
  Vi_inv_pc_1[[i]] <- solve(Vi_pc_1[[i]])  
  Vi_inv_pc_2[[i]] <- solve(Vi_pc_2[[i]])  
  Vi_inv_pc_3[[i]] <- solve(Vi_pc_3[[i]])  
  Vi_inv_pc_4[[i]] <- solve(Vi_pc_4[[i]])  
  
}


# check these lists: 
sapply(Vi_pc_1, dim)
sapply(Vi_pc_2, dim)
sapply(Vi_pc_3, dim)
sapply(Vi_pc_4, dim)


sapply(Vi_inv_pc_1, dim)
sapply(Vi_inv_pc_2, dim)
sapply(Vi_inv_pc_3, dim)
sapply(Vi_inv_pc_4, dim)


# just check: not needed for MCMC; but they are really small and do not have to remove after using
sapply(CY_pc_1, dim)
sapply(CY_pc_2, dim)
sapply(CY_pc_3, dim)
sapply(CY_pc_4, dim)


summary(Tau2hat_1)   # scale parameters in the original space
summary(Tau2hat_2)   # scale parameters in the original space
summary(Tau2hat_3)   # scale parameters in the original space
summary(Tau2hat_4)   # scale parameters in the original space


##  Check fitted on-site surrogates hyperparameters: need tau2hats for this checking

summary(MLEs_Kdir_pc)
summary(MLEs_Kcross_pc)
summary(MLEs_cdir_pc)
summary(MLEs_ccross_pc)


## Length-scales: 
boxplot(MLEs_Kdir_pc[, 1:4], las=2, horizontal = TRUE, main = "Boxplots of fitted OSSs length-scales for the YM_PC1")
boxplot(MLEs_Kcross_pc[, 1:4], las=2, horizontal = TRUE, main = "Boxplots of fitted OSSs length-scales for the YM_PC1")
boxplot(MLEs_cdir_pc[, 1:4], las=2, horizontal = TRUE, main = "Boxplots of fitted OSSs length-scales for the YM_PC1")
boxplot(MLEs_ccross_pc[, 1:4], las=2, horizontal = TRUE, main = "Boxplots of fitted OSSs length-scales for the YM_PC1")



## Nugget effect: mostly noise free: 
boxplot(MLEs_Kdir_pc[, 5], main="Boxplot of fitted OSSs nuggets g for the combination of 4 frequencies")
boxplot(MLEs_Kcross_pc[, 5], main="Boxplot of fitted OSSs nuggets g for the combination of 4 frequencies")
boxplot(MLEs_cdir_pc[, 5], main="Boxplot of fitted OSSs nuggets g for the combination of 4 frequencies")
boxplot(MLEs_ccross_pc[, 5], main="Boxplot of fitted OSSs nuggets g for the combination of 4 frequencies")


## Check overall noise level by combining the nuggets and scales: the left noise in yM are very small now!!
Noise <- MLEs_Kdir_pc[, 5] * Tau2hat_1
matplot(Noise,  type = c("l"), pch=1, col = 1, xlab="Site index", ylab="Estimated on-site noise level", 
        main = "Fitted OSSs noise level g*Tau^2 for Kdir (in PC)")

Noise <- MLEs_Kcross_pc[, 5] * Tau2hat_2
matplot(Noise,  type = c("l"), pch=1, col = 1, xlab="Site index", ylab="Estimated on-site noise level", 
        main = "Fitted OSSs noise level g*Tau^2 for Kcross (in PC)")


Noise <- MLEs_cdir_pc[, 5] * Tau2hat_3
matplot(Noise,  type = c("l"), pch=1, col = 1, xlab="Site index", ylab="Estimated on-site noise level", 
        main = "Fitted OSSs noise level g*Tau^2 for Cdir (in PC)")

Noise <- MLEs_ccross_pc[, 5] * Tau2hat_4
matplot(Noise,  type = c("l"), pch=1, col = 1, xlab="Site index", ylab="Estimated on-site noise level", 
        main = "Fitted OSSs noise level g*Tau^2 for Ccross (in PC)")





## ---------------------------------------------------------------------------------------
##                       1.5 Check fitted discrepancy in first PC space: 
## ---------------------------------------------------------------------------------------

## ------------------------------------------------------------------------
# generate priors for discrepancy term: 
da <- darg(list(mle=TRUE, max=5, min=eps), X, samp.size = nrow(X)); da  # also have stronger regularization on bias term

## prior for nugget effects: need to be careful in principal component space -----------------------------
## ga <- garg(list(mle=TRUE, max=1), Bias_pr$x[,1]); ga
ga_1 <- garg(list(mle=TRUE, max=1), Bias_pr_1$x[,1]); ga_1
ga_2 <- garg(list(mle=TRUE, max=1), Bias_pr_2$x[,1]); ga_2
ga_3 <- garg(list(mle=TRUE, max=1), Bias_pr_3$x[,1]); ga_3
ga_4 <- garg(list(mle=TRUE, max=1), Bias_pr_4$x[,1]); ga_4

## load the optimization results: 
load("RData/all_outputs_cali_opti_betaprior.RData", verbose=TRUE)


# best fitted U value in [0, 1]^4 from optimization: 
## sol1 <- sol[which.min(sol$value),]; sol1

sol1 <- sol[order(sol$value),][2, ]  ## use the second one
u <- as.numeric(sol1[5:8]) ; u

# fit for the bias term in PC space:
# only one u at a time in this case: 
## yM_pc <- rep(NA, 292)

yM_pc_1 <- rep(NA, 292)
yM_pc_2 <- rep(NA, 292)
yM_pc_3 <- rep(NA, 292)
yM_pc_4 <- rep(NA, 292)


# loop over 292 Xfield
for (i in 1:292){
##   CX_pc <- covar.sep(u, U[list_na_index[[i]], ], d = MLEs_Kdir_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
##   yM_pc[i] <- CX_pc %*% CY_pc[[i]]  # predict from ym for means
  CX_pc_1 <- covar.sep(u, U[list_na_index[[i]], ], d = MLEs_Kdir_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
  yM_pc_1[i] <- CX_pc_1 %*% CY_pc_1[[i]]  # predict from ym for means  
  
  CX_pc_2 <- covar.sep(u, U[list_na_index[[i]], ], d = MLEs_Kcross_pc[i, 1:4], g = MLEs_Kcross_pc[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
  yM_pc_2[i] <- CX_pc_2 %*% CY_pc_2[[i]]  # predict from ym for means  
  
  CX_pc_3 <- covar.sep(u, U[list_na_index[[i]], ], d = MLEs_cdir_pc[i, 1:4], g = MLEs_cdir_pc[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
  yM_pc_3[i] <- CX_pc_3 %*% CY_pc_3[[i]]  # predict from ym for means  
  
  CX_pc_4 <- covar.sep(u, U[list_na_index[[i]], ], d = MLEs_ccross_pc[i, 1:4], g = MLEs_ccross_pc[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
  yM_pc_4[i] <- CX_pc_4 %*% CY_pc_4[[i]]  # predict from ym for means  
  
}


# regularize with a marginal beta(2, 2) prior on U:
Uprior <-  - sum(dbeta(u, 2, 2, log = TRUE))  # negative log likelihood to minimize: be careful with being positive in MCMC

# prior option 2: weak wrong prior (far away fro the truth)
Uprior <-  - (dbeta(u[1], 1.5, 2.5, log = TRUE) + dbeta(u[2], 1.5, 2.5, log = TRUE) + 
                dbeta(u[3], 2.5, 1.5, log = TRUE) + dbeta(u[4], 2.5, 1.5, log = TRUE))  # negative log likelihood to minimize

# prior option 3: weak prescient prior (close to the truth)
Uprior <-  - (dbeta(u[1], 2.5, 1.5, log = TRUE) + dbeta(u[2], 2.5, 1.5, log = TRUE) + 
                dbeta(u[3], 1.5, 2.5, log = TRUE) + dbeta(u[4], 1.5, 2.5, log = TRUE))  # negative log likelihood to minimize


## colnames(YF) <- c("y2", "y5", "y9", "y11")
## yF_pc <- predict(Bias_pr, YF)[, 1]

colnames(YF_1) <- c("y2", "y5", "y9", "y11")
colnames(YF_2) <- c("y2", "y5", "y9", "y11")
colnames(YF_3) <- c("y2", "y5", "y9", "y11")
colnames(YF_4) <- c("y2", "y5", "y9", "y11")
yF_pc_1 <- predict(Bias_pr_1, YF_1)[, 1]
yF_pc_2 <- predict(Bias_pr_2, YF_2)[, 1]
yF_pc_3 <- predict(Bias_pr_3, YF_3)[, 1]
yF_pc_4 <- predict(Bias_pr_4, YF_4)[, 1]


## use a fixed global PCA on every new U values for fair comparison: 
## bias_u <-  yF_pc - yM_pc    ##  observed bias given U value
bias_u_1 <-  yF_pc_1 - yM_pc_1    ##  observed bias given U value
bias_u_2 <-  yF_pc_2 - yM_pc_2    ##  observed bias given U value
bias_u_3 <-  yF_pc_3 - yM_pc_3    ##  observed bias given U value
bias_u_4 <-  yF_pc_4 - yM_pc_4    ##  observed bias given U value



## Visualize observed discrepancy at u_0
matplot(bias_u_1, 
        type = c("l"), pch=1, col = 1, xlab="Site index", ylab="Observed discrepancy", 
        main = "Observed yF-yM(u0) in 1st P.C. space across site index")

matplot(data.frame(bias_u_1, bias_u_2, bias_u_3, bias_u_4),
        type = c("l"), pch=1, col = 1, xlab="Site index", ylab="Observed discrepancy", 
        main = "Observed yF-yM(u0) in 1st P.C. space across site index")

boxplot(data.frame(bias_u_1, bias_u_2, bias_u_3, bias_u_4))


## Fit GP on the first principal component of all biases 
## bhat <- newGPsep(X, bias_u , d=rep(da$start, ncol(X)), g=ga$start, dK=TRUE)

bhat_1 <- newGPsep(X, bias_u_1 , d=rep(da$start, ncol(X)), g=ga_1$start, dK=TRUE)
bhat_2 <- newGPsep(X, bias_u_2 , d=rep(da$start, ncol(X)), g=ga_2$start, dK=TRUE)
bhat_3 <- newGPsep(X, bias_u_3 , d=rep(da$start, ncol(X)), g=ga_3$start, dK=TRUE)
bhat_4 <- newGPsep(X, bias_u_4 , d=rep(da$start, ncol(X)), g=ga_4$start, dK=TRUE)



#if(ga$mle) cmle <- jmleGPsep(bhat, drange=c(da$min, da$max), grange=c(ga$min, ga$max), dab=da$ab, gab=ga$ab)
## mle_b <- mleGPsep(bhat, param="both", tmin=c(da$min, ga$min), tmax=c(da$max, ga$max), ab=c(da$ab, ga$ab), maxit=1000)

mle_b_1 <- mleGPsep(bhat_1, param="both", tmin=c(da$min, ga_1$min), tmax=c(da$max, ga_1$max), ab=c(da$ab, ga_1$ab), maxit=1000)
mle_b_2 <- mleGPsep(bhat_2, param="both", tmin=c(da$min, ga_2$min), tmax=c(da$max, ga_2$max), ab=c(da$ab, ga_2$ab), maxit=1000)
mle_b_3 <- mleGPsep(bhat_3, param="both", tmin=c(da$min, ga_3$min), tmax=c(da$max, ga_3$max), ab=c(da$ab, ga_3$ab), maxit=1000)
mle_b_4 <- mleGPsep(bhat_4, param="both", tmin=c(da$min, ga_4$min), tmax=c(da$max, ga_4$max), ab=c(da$ab, ga_4$ab), maxit=1000)



## Check negative log-likelihood: 
## - llikGPsep(bhat, dab=da$ab, gab=ga$ab) +  Uprior
- llikGPsep(bhat_1, dab=da$ab, gab=ga_1$ab) +  Uprior
- llikGPsep(bhat_2, dab=da$ab, gab=ga_2$ab) +  Uprior
- llikGPsep(bhat_3, dab=da$ab, gab=ga_3$ab) +  Uprior
- llikGPsep(bhat_4, dab=da$ab, gab=ga_4$ab) +  Uprior


print(mle_b_1)
print(mle_b_2)
print(mle_b_3)
print(mle_b_4)


# estimate the tau2bhat for bias term: 
## Cd <- covar.sep(X, X, d=mle_b$theta[1:13], g=mle_b$theta[14])  # this part need nugget g
## Cd_inv <- solve(Cd)
## taub2hat <- drop(t(bias_u) %*% Cd_inv  %*% (bias_u ) / length(bias_u)); taub2hat
## Vd <- taub2hat * Cd # the covariance of the bias term keeps the same with new proposed theta

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



## Discrepancy appears to be more stable
par(mfrow = c(1, 2))
boxplot(bias_u_1, main=paste("Observed bias, tau^2=", round(taub2hat_1, 3)))
qqnorm(bias_u_1)
qqline(bias_u_1)

par(mfrow = c(1, 2))
boxplot(bias_u_4, main=paste("Observed bias, tau^2=", round(taub2hat_4, 3)))
qqnorm(bias_u_4)
qqline(bias_u_4)

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
## colnames(YM) <- c("y2", "y5", "y9", "y11")  ## name it for PCA rotation
## ym_pc  <- predict(Bias_pr, YM)
## ym_pc1 <- ym_pc[, 1]


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
## colnames(YF) <- c("y2", "y5", "y9", "y11")
## yf_pc  <- predict(Bias_pr, YF)
## yf_pc1 <- yf_pc[, 1]

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





## ------------------------------------------------------------------------
##    2. Posterior of U conditioning on the estimated hyperparameters:
## ------------------------------------------------------------------------


## start with the initial posterior: for U0 

u0 <- u


## ------------------------------------------------------------------------
##   # 2.1. inverse part: 
## ------------------------------------------------------------------------



## Now, need to do this calculation for 4 times with multiple GPs on frequencies:

## Vo_inv_(Vob)^T: updated because Vob is a function of u

## get Vob^T first: this is changing with u at each iteration

## Vob_pc <- list()  ## combined in 1st PC space
Vob_pc_1 <- list()  ## combined in 1st PC space
Vob_pc_2 <- list()  ## combined in 1st PC space
Vob_pc_3 <- list()  ## combined in 1st PC space
Vob_pc_4 <- list()  ## combined in 1st PC space



for (i in 1:292) {   
  ## Cob_u_i_pc <-  covar.sep(u0, U[list_na_index[[i]], ],  d = MLEs_Kdir_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) # 1 by 1000
  
  Cob_u_i_pc_1 <-  covar.sep(u0, U[list_na_index[[i]], ],  d = MLEs_Kdir_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) # 1 by 1000
  Cob_u_i_pc_2 <-  covar.sep(u0, U[list_na_index[[i]], ],  d = MLEs_Kcross_pc[i, 1:4], g = MLEs_Kcross_pc[i, 5]) # 1 by 1000
  Cob_u_i_pc_3 <-  covar.sep(u0, U[list_na_index[[i]], ],  d = MLEs_cdir_pc[i, 1:4], g = MLEs_cdir_pc[i, 5]) # 1 by 1000
  Cob_u_i_pc_4 <-  covar.sep(u0, U[list_na_index[[i]], ],  d = MLEs_ccross_pc[i, 1:4], g = MLEs_ccross_pc[i, 5]) # 1 by 1000
  
  ## Combined: 
  ## Vob_pc[[i]] <-  t(Tau2hat[i, 1]  * Cob_u_i_pc)
  
  Vob_pc_1[[i]] <-  t(Tau2hat_1[i, 1]  * Cob_u_i_pc_1)
  Vob_pc_2[[i]] <-  t(Tau2hat_2[i, 1]  * Cob_u_i_pc_2)
  Vob_pc_3[[i]] <-  t(Tau2hat_3[i, 1]  * Cob_u_i_pc_3)
  Vob_pc_4[[i]] <-  t(Tau2hat_4[i, 1]  * Cob_u_i_pc_4)
  
}

sapply(Vob_pc_1, dim)
sapply(Vob_pc_2, dim)
sapply(Vob_pc_3, dim)
sapply(Vob_pc_4, dim)


## Vo_inv_(Vob)^T
## Vo_inv_Vob_t_pc <- list()
Vo_inv_Vob_t_pc_1 <- list()
Vo_inv_Vob_t_pc_2 <- list()
Vo_inv_Vob_t_pc_3 <- list()
Vo_inv_Vob_t_pc_4 <- list()

for (i in 1:292) {
  ## Vo_inv_Vob_t_pc[[i]] <- Vi_inv_pc[[i]] %*% Vob_pc[[i]]   # each part is 1000 by 1
  Vo_inv_Vob_t_pc_1[[i]] <- Vi_inv_pc_1[[i]] %*% Vob_pc_1[[i]]   # each part is 1000 by 1
  Vo_inv_Vob_t_pc_2[[i]] <- Vi_inv_pc_2[[i]] %*% Vob_pc_2[[i]]   # each part is 1000 by 1
  Vo_inv_Vob_t_pc_3[[i]] <- Vi_inv_pc_3[[i]] %*% Vob_pc_3[[i]]   # each part is 1000 by 1
  Vo_inv_Vob_t_pc_4[[i]] <- Vi_inv_pc_4[[i]] %*% Vob_pc_4[[i]]   # each part is 1000 by 1
}


sapply(Vo_inv_Vob_t_pc_1, dim)
sapply(Vo_inv_Vob_t_pc_2, dim)
sapply(Vo_inv_Vob_t_pc_3, dim)
sapply(Vo_inv_Vob_t_pc_4, dim)




## C_in_u: also updated with u

##  Vob_Vo_inv_Vbo: 292 by 292 diagonal matrix
## Vob_Vo_inv_Vbo_pc <- matrix(0, nrow=292, ncol=292)
Vob_Vo_inv_Vbo_pc_1 <- matrix(0, nrow=292, ncol=292)
Vob_Vo_inv_Vbo_pc_2 <- matrix(0, nrow=292, ncol=292)
Vob_Vo_inv_Vbo_pc_3 <- matrix(0, nrow=292, ncol=292)
Vob_Vo_inv_Vbo_pc_4 <- matrix(0, nrow=292, ncol=292)



for (i in 1:292) {
  ## Vob_Vo_inv_Vbo_pc[i, i] <- t(Vob_pc[[i]]) %*% Vi_inv_pc[[i]] %*% Vob_pc[[i]]   # each part is a scalar
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
  ## vI_nf[i, i]  <- covar.sep(u,  t(u), d = MLEs_Kdir_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) * Tau2hat[i, 1]    # this part has nugget g

  vI_nf_1[i, i]  <- covar.sep(u,  t(u), d = MLEs_Kdir_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) * Tau2hat_1[i, 1]    # this part has nugget g
  vI_nf_2[i, i]  <- covar.sep(u,  t(u), d = MLEs_Kcross_pc[i, 1:4], g = MLEs_Kcross_pc[i, 5]) * Tau2hat_2[i, 1]    # this part has nugget g
  vI_nf_3[i, i]  <- covar.sep(u,  t(u), d = MLEs_cdir_pc[i, 1:4], g = MLEs_cdir_pc[i, 5]) * Tau2hat_3[i, 1]    # this part has nugget g
  vI_nf_4[i, i]  <- covar.sep(u,  t(u), d = MLEs_ccross_pc[i, 1:4], g = MLEs_ccross_pc[i, 5]) * Tau2hat_4[i, 1]    # this part has nugget g
   }


summary(diag(vI_nf_1))
summary(diag(vI_nf_2))
summary(diag(vI_nf_3))
summary(diag(vI_nf_4))

summary(diag(Vd_1))
summary(diag(Vd_2))
summary(diag(Vd_3))
summary(diag(Vd_4))


vI_nf_1[1:10, 1:10]
diag(vI_nf_1)


## Vb_pc <-  vI_nf   +  Vd
Vb_pc_1 <-  vI_nf_1   +  Vd_1
Vb_pc_2 <-  vI_nf_2   +  Vd_2
Vb_pc_3 <-  vI_nf_3   +  Vd_3
Vb_pc_4 <-  vI_nf_4   +  Vd_4


## C_u <-  Vb_pc - Vob_Vo_inv_Vbo_pc
C_u_1 <-  Vb_pc_1 - Vob_Vo_inv_Vbo_pc_1
C_u_2 <-  Vb_pc_2 - Vob_Vo_inv_Vbo_pc_2
C_u_3 <-  Vb_pc_3 - Vob_Vo_inv_Vbo_pc_3
C_u_4 <-  Vb_pc_4 - Vob_Vo_inv_Vbo_pc_4


## C_u_inv <- solve(C_u)
C_u_inv_1 <- solve(C_u_1)
C_u_inv_2 <- solve(C_u_2)
C_u_inv_3 <- solve(C_u_3)
C_u_inv_4 <- solve(C_u_4)


## The quadratic form for yM and yF: be careful with a shift for mean now: 

# t(d)V_inv(d)


## t(z) C_u_inv (z) = c1
## c1 <- t(yf_pc1) %*% C_u_inv %*% (yf_pc1) ; c1
c1_1 <- t(yf_pc1_1) %*% C_u_inv_1 %*% (yf_pc1_1) ; c1_1
c1_2 <- t(yf_pc1_2) %*% C_u_inv_2 %*% (yf_pc1_2) ; c1_2
c1_3 <- t(yf_pc1_3) %*% C_u_inv_3 %*% (yf_pc1_3) ; c1_3
c1_4 <- t(yf_pc1_4) %*% C_u_inv_4 %*% (yf_pc1_4) ; c1_4


### fixed part: c2 = t(y - m) A11_inv (y-m), only need to be calculated once. 
## c2 <- rep(NA, 292)
c2_1 <- rep(NA, 292)
c2_2 <- rep(NA, 292)
c2_3 <- rep(NA, 292)
c2_4 <- rep(NA, 292)



for (i in 1:292) {
  
  ## simply use the precalculated pc 1 directly: 
  
  
  ## yi_pc1 <- ym_pc1[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_pc1_1 <- ym_pc1[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_pc1_2 <- ym_pc2[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_pc1_3 <- ym_pc3[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_pc1_4 <- ym_pc4[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  
  
  ## c2[i] <- t( yi_pc1 ) %*% Vi_inv_pc[[i]]  %*%  yi_pc1 
  c2_1[i] <- t( yi_pc1_1 ) %*% Vi_inv_pc_1[[i]]  %*%  yi_pc1_1 
  c2_2[i] <- t( yi_pc1_2 ) %*% Vi_inv_pc_2[[i]]  %*%  yi_pc1_2 
  c2_3[i] <- t( yi_pc1_3 ) %*% Vi_inv_pc_3[[i]]  %*%  yi_pc1_3 
  c2_4[i] <- t( yi_pc1_4 ) %*% Vi_inv_pc_4[[i]]  %*%  yi_pc1_4 
}

## c2 <- sum(c2); c2 # just the length of y since tau is estimated by t(y) Cin y / length(y)
c2_1 <- sum(c2_1); c2_1 
c2_2 <- sum(c2_2); c2_2 
c2_3 <- sum(c2_3); c2_3 
c2_4 <- sum(c2_4); c2_4 


### updating part: 
## c3 <- rep(NA, 292)
c3_1 <- rep(NA, 292)
c3_2 <- rep(NA, 292)
c3_3 <- rep(NA, 292)
c3_4 <- rep(NA, 292)

## c6 <- rep(NA, 292)
c6_1 <- rep(NA, 292)
c6_2 <- rep(NA, 292)
c6_3 <- rep(NA, 292)
c6_4 <- rep(NA, 292)


for (i in 1:292) {
  ## calculate t(y) A11_inv A_12 by list
  ## yi_pc1 <- ym_pc1[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_pc1_1 <- ym_pc1[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_pc1_2 <- ym_pc2[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_pc1_3 <- ym_pc3[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  yi_pc1_4 <- ym_pc4[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  
  
  ## c3[i] <- t(yi_pc1) %*% Vo_inv_Vob_t_pc[[i]] 
  c3_1[i] <- t(yi_pc1_1) %*% Vo_inv_Vob_t_pc_1[[i]] 
  c3_2[i] <- t(yi_pc1_2) %*% Vo_inv_Vob_t_pc_2[[i]] 
  c3_3[i] <- t(yi_pc1_3) %*% Vo_inv_Vob_t_pc_3[[i]] 
  c3_4[i] <- t(yi_pc1_4) %*% Vo_inv_Vob_t_pc_4[[i]] 
  
  
  ## t(z)V21y = c5 part: 
  ## c6[i] <- t(Vo_inv_Vob_t_pc[[i]]) %*% yi_pc1
  c6_1[i] <- t(Vo_inv_Vob_t_pc_1[[i]]) %*% yi_pc1_1
  c6_2[i] <- t(Vo_inv_Vob_t_pc_2[[i]]) %*% yi_pc1_2
  c6_3[i] <- t(Vo_inv_Vob_t_pc_3[[i]]) %*% yi_pc1_3
  c6_4[i] <- t(Vo_inv_Vob_t_pc_4[[i]]) %*% yi_pc1_4
  
  }

## c4 <- t(c3) %*% C_u_inv %*% c3  ## the whole updating part for the 1st quadratic form
c4_1 <- t(c3_1) %*% C_u_inv_1 %*% c3_1  ## the whole updating part for the 1st quadratic form
c4_2 <- t(c3_2) %*% C_u_inv_2 %*% c3_2  ## the whole updating part for the 1st quadratic form
c4_3 <- t(c3_3) %*% C_u_inv_3 %*% c3_3  ## the whole updating part for the 1st quadratic form
c4_4 <- t(c3_4) %*% C_u_inv_4 %*% c3_4  ## the whole updating part for the 1st quadratic form

## c7 <- rep(NA, 292) ## 
c7_1 <- rep(NA, 292) ## 
c7_2 <- rep(NA, 292) ## 
c7_3 <- rep(NA, 292) ## 
c7_4 <- rep(NA, 292) ## 


for (i in 1:292) {
  ## c7[i] <- t(C_u_inv[, i]) %*% c6
  c7_1[i] <- t(C_u_inv_1[, i]) %*% c6_1
  c7_2[i] <- t(C_u_inv_2[, i]) %*% c6_2
  c7_3[i] <- t(C_u_inv_3[, i]) %*% c6_3
  c7_4[i] <- t(C_u_inv_4[, i]) %*% c6_4
  }

## c5 is the whole t(z) V21 y mixed z y quadratic part: 
## c5 <- - t(yf_pc1) %*% c7

c5_1 <- - t(yf_pc1_1) %*% c7_1
c5_2 <- - t(yf_pc1_2) %*% c7_2
c5_3 <- - t(yf_pc1_3) %*% c7_3
c5_4 <- - t(yf_pc1_4) %*% c7_4



## ------------------------------------------------------------------------
##   2.2. determinant part: 
## ------------------------------------------------------------------------

## fixed part: log(det(A11))
## c8 <- rep(NA, 292)

c8_1 <- rep(NA, 292)
c8_2 <- rep(NA, 292)
c8_3 <- rep(NA, 292)
c8_4 <- rep(NA, 292)


for (i in 1:292){
##   c8[i] <- determinant(Vi_pc[[i]], logarithm=TRUE)$modulus
  c8_1[i] <- determinant(Vi_pc_1[[i]], logarithm=TRUE)$modulus
  c8_2[i] <- determinant(Vi_pc_2[[i]], logarithm=TRUE)$modulus
  c8_3[i] <- determinant(Vi_pc_3[[i]], logarithm=TRUE)$modulus
  c8_4[i] <- determinant(Vi_pc_4[[i]], logarithm=TRUE)$modulus
  
}


## c9 <- sum(c8); c9
c9_1 <- sum(c8_1); c9_1
c9_2 <- sum(c8_2); c9_2
c9_3 <- sum(c8_3); c9_3
c9_4 <- sum(c8_4); c9_4



## updated part: log(det(C_u))
## c10 <- determinant(C_u, logarithm=TRUE)$modulus; c10
c10_1 <- determinant(C_u_1, logarithm=TRUE)$modulus; c10_1
c10_2 <- determinant(C_u_2, logarithm=TRUE)$modulus; c10_2
c10_3 <- determinant(C_u_3, logarithm=TRUE)$modulus; c10_3
c10_4 <- determinant(C_u_4, logarithm=TRUE)$modulus; c10_4



## ------------------------------------------------------------------------
##   2.3. log posterior: put all parts together
## ------------------------------------------------------------------------


# log likelihood: 
## log_like <- -1/2 *(c9 + c10) - (c1 + c2+ c4 + 2* c5) / 2; log_like
log_like <- (-1/2 *(c9_1 + c10_1) - (c1_1 + c2_1 + c4_1 + 2* c5_1) / 2 
            -1/2 *(c9_2 + c10_2) - (c1_2 + c2_2 + c4_2 + 2* c5_2) / 2 
            -1/2 *(c9_3 + c10_3) - (c1_3 + c2_3 + c4_3 + 2* c5_3) / 2 
            -1/2 *(c9_4 + c10_4) - (c1_4 + c2_4 + c4_4 + 2* c5_4) / 2 ) 

log_like


# log prior on theta: 
## logprior <-  sum(dbeta(u0, 2, 2, log = TRUE))  # beta(2, 2) prior
#logprior <-  sum(dunif(theta0, 0, 1, log = TRUE))  # uniform(0, 1) prior

# prior option 2: weak wrong prior (far away fro the truth)
logprior <-   (dbeta(u0[1], 1.5, 2.5, log = TRUE) + dbeta(u0[2], 1.5, 2.5, log = TRUE) + 
                dbeta(u0[3], 2.5, 1.5, log = TRUE) + dbeta(u0[4], 2.5, 1.5, log = TRUE)) 


# prior option 3: prescient prior
logprior <-   (dbeta(u0[1], 2.5, 1.5, log = TRUE) + dbeta(u0[2], 2.5, 1.5, log = TRUE) + 
                 dbeta(u0[3], 1.5, 2.5, log = TRUE) + dbeta(u0[4], 1.5, 2.5, log = TRUE)) 

# log posterior
logpost_u <- logprior + log_like; logpost_u



## ---------------------------------------------------------
##      3. MCMC Sampling loop: 
## ---------------------------------------------------------


its <- 100000 # number of total iterations

Us <- matrix(NA, nrow=its, ncol=4)
Us[1, ] <- u0

logposts <- rep(NA, its)
logposts[1] <- logpost_u

# use different step sizes for different theta: 
#steps <- c(.2, .2, .2, .2) # use a different step size for different theta
steps <- c(.1, .1, .1, .1) # use a different step size for different theta

for (iter in 2:its){
  
  ut <- Us[iter - 1, ]
  
  ## proposal:  
  # Gibbs sampling for each dimension:
  u_cand <- ut
  remainder <- iter %% 4 +1
  u_cand[remainder] <- ut[remainder] + rnorm(1, 0, steps[remainder]) 
  
  #theta_cand[remainder] <- runif(1)
  
  # compute the log posterior: 
  
  ## Vob_pc1 part: updated
  
  #Vob_2 <- list()
  #Vob_5 <- list()
  #Vob_9 <- list()
  #Vob_11 <- list()
  #Vob_pc1 <- list()  ## combined in 1st PC space
  
  ## u_cand enter the covariance matrix here: through V_ob only!!!
  for (i in 1:292) {   
 ## Cob_u_i_2 <-  covar.sep(u_cand, U[list_na_index[[i]], ], d = MLEs_Kdir_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) #* (w1[1]/s_pc[1])^2# 1 by 1000
 ## Vob_pc[[i]] <-  t( Tau2hat[i, 1]  * Cob_u_i_2 ) 
    
    Cob_u_i_pc_1 <-  covar.sep(u_cand, U[list_na_index[[i]], ],  d = MLEs_Kdir_pc[i, 1:4], g = MLEs_Kdir_pc[i, 5]) # 1 by 1000
    Cob_u_i_pc_2 <-  covar.sep(u_cand, U[list_na_index[[i]], ],  d = MLEs_Kcross_pc[i, 1:4], g = MLEs_Kcross_pc[i, 5]) # 1 by 1000
    Cob_u_i_pc_3 <-  covar.sep(u_cand, U[list_na_index[[i]], ],  d = MLEs_cdir_pc[i, 1:4], g = MLEs_cdir_pc[i, 5]) # 1 by 1000
    Cob_u_i_pc_4 <-  covar.sep(u_cand, U[list_na_index[[i]], ],  d = MLEs_ccross_pc[i, 1:4], g = MLEs_ccross_pc[i, 5]) # 1 by 1000
    
    Vob_pc_1[[i]] <-  t(Tau2hat_1[i, 1]  * Cob_u_i_pc_1)
    Vob_pc_2[[i]] <-  t(Tau2hat_2[i, 1]  * Cob_u_i_pc_2)
    Vob_pc_3[[i]] <-  t(Tau2hat_3[i, 1]  * Cob_u_i_pc_3)
    Vob_pc_4[[i]] <-  t(Tau2hat_4[i, 1]  * Cob_u_i_pc_4)
    
    }
  
  
  ## Vo_inv_Vob_t_pc <- list()
  
  ## for (i in 1:292) {
  ##   Vo_inv_Vob_t_pc[[i]] <- Vi_inv_pc[[i]] %*% Vob_pc[[i]]   # each part is 1000 by 1
  ## }
  
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
  ## Vob_Vo_inv_Vbo_pc <- matrix(0, nrow=292, ncol=292)
  
  ## for (i in 1:292) {
  ##  Vob_Vo_inv_Vbo_pc[i, i] <- t(Vob_pc[[i]]) %*% Vi_inv_pc[[i]] %*% Vob_pc[[i]]   # each part is a scalar
  ## }
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
  
  
  ## C_u <-  Vb_pc - Vob_Vo_inv_Vbo_pc
  ## C_u_inv <- solve(C_u)
  
  ## C_u <-  Vb_pc - Vob_Vo_inv_Vbo_pc
  C_u_1 <-  Vb_pc_1 - Vob_Vo_inv_Vbo_pc_1
  C_u_2 <-  Vb_pc_2 - Vob_Vo_inv_Vbo_pc_2
  C_u_3 <-  Vb_pc_3 - Vob_Vo_inv_Vbo_pc_3
  C_u_4 <-  Vb_pc_4 - Vob_Vo_inv_Vbo_pc_4
  
  ## C_u_inv <- solve(C_u)
  C_u_inv_1 <- solve(C_u_1)
  C_u_inv_2 <- solve(C_u_2)
  C_u_inv_3 <- solve(C_u_3)
  C_u_inv_4 <- solve(C_u_4)
  
  ## The quadratic form for yM and yF: be careful with a shift for mean now: 
  
  # t(d)V_inv(d)
  ## t(z) C_u_inv (z) = c1
  ## c1 <- t(yf_pc1) %*% C_u_inv %*% (yf_pc1)
  
  c1_1 <- t(yf_pc1_1) %*% C_u_inv_1 %*% (yf_pc1_1) 
  c1_2 <- t(yf_pc1_2) %*% C_u_inv_2 %*% (yf_pc1_2) 
  c1_3 <- t(yf_pc1_3) %*% C_u_inv_3 %*% (yf_pc1_3) 
  c1_4 <- t(yf_pc1_4) %*% C_u_inv_4 %*% (yf_pc1_4) 
  
  
  
  ### c2 is constant of u
  
  ### updating part: 
  ## c3 <- rep(NA, 292)
  ## c6 <- rep(NA, 292)
  
  ## for (i in 1:292) {
    ## calculate t(y) A11_inv A_12 by list
  ##   yi_pc1 <- ym_pc1[1:1000 + 1000 * (i-1)][list_na_index[[i]]]
  ##   c3[i] <- t( yi_pc1 ) %*% Vo_inv_Vob_t_pc[[i]] 
    
    ## t(z)V21y = c5 part: 
  ##   c6[i] <- t(Vo_inv_Vob_t_pc[[i]]) %*%  yi_pc1 
  ## }
  
  ### updating part: 
  ## c3 <- rep(NA, 292)
  c3_1 <- rep(NA, 292)
  c3_2 <- rep(NA, 292)
  c3_3 <- rep(NA, 292)
  c3_4 <- rep(NA, 292)
  
  ## c6 <- rep(NA, 292)
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
    
    
    ## c3[i] <- t(yi_pc1) %*% Vo_inv_Vob_t_pc[[i]] 
    c3_1[i] <- t(yi_pc1_1) %*% Vo_inv_Vob_t_pc_1[[i]] 
    c3_2[i] <- t(yi_pc1_2) %*% Vo_inv_Vob_t_pc_2[[i]] 
    c3_3[i] <- t(yi_pc1_3) %*% Vo_inv_Vob_t_pc_3[[i]] 
    c3_4[i] <- t(yi_pc1_4) %*% Vo_inv_Vob_t_pc_4[[i]] 
    
    
    ## t(z)V21y = c5 part: 
    ## c6[i] <- t(Vo_inv_Vob_t_pc[[i]]) %*% yi_pc1
    c6_1[i] <- t(Vo_inv_Vob_t_pc_1[[i]]) %*% yi_pc1_1
    c6_2[i] <- t(Vo_inv_Vob_t_pc_2[[i]]) %*% yi_pc1_2
    c6_3[i] <- t(Vo_inv_Vob_t_pc_3[[i]]) %*% yi_pc1_3
    c6_4[i] <- t(Vo_inv_Vob_t_pc_4[[i]]) %*% yi_pc1_4
    
  }
  
  
  ## c4 <- t(c3) %*% C_u_inv %*% c3  ## the whole updating part for the 1st quadratic form
  
  ## c4 <- t(c3) %*% C_u_inv %*% c3  ## the whole updating part for the 1st quadratic form
  c4_1 <- t(c3_1) %*% C_u_inv_1 %*% c3_1  ## the whole updating part for the 1st quadratic form
  c4_2 <- t(c3_2) %*% C_u_inv_2 %*% c3_2  ## the whole updating part for the 1st quadratic form
  c4_3 <- t(c3_3) %*% C_u_inv_3 %*% c3_3  ## the whole updating part for the 1st quadratic form
  c4_4 <- t(c3_4) %*% C_u_inv_4 %*% c3_4  ## the whole updating part for the 1st quadratic form
  
  ## c7 <- rep(NA, 292) ## 
  
  ## c7 <- rep(NA, 292) ## 
  c7_1 <- rep(NA, 292) ## 
  c7_2 <- rep(NA, 292) ## 
  c7_3 <- rep(NA, 292) ## 
  c7_4 <- rep(NA, 292) ## 
  
  
  ## for (i in 1:292) {
  ##   c7[i] <- t(C_u_inv[, i]) %*% c6
  ## }
  
  for (i in 1:292) {
    ## c7[i] <- t(C_u_inv[, i]) %*% c6
    c7_1[i] <- t(C_u_inv_1[, i]) %*% c6_1
    c7_2[i] <- t(C_u_inv_2[, i]) %*% c6_2
    c7_3[i] <- t(C_u_inv_3[, i]) %*% c6_3
    c7_4[i] <- t(C_u_inv_4[, i]) %*% c6_4
  }
  
  ## c5 is the whole t(z) V21 y mixed z y quadratic part: 
  ## c5 <- - t(yf_pc1) %*% c7
  
  c5_1 <- - t(yf_pc1_1) %*% c7_1
  c5_2 <- - t(yf_pc1_2) %*% c7_2
  c5_3 <- - t(yf_pc1_3) %*% c7_3
  c5_4 <- - t(yf_pc1_4) %*% c7_4
  
  
  
  ## updated part: log(det(C_u))
  ## c10 <- determinant(C_u , logarithm=TRUE)$modulus
  
  ## updated part: log(det(C_u))
  ## c10 <- determinant(C_u, logarithm=TRUE)$modulus; c10
  c10_1 <- determinant(C_u_1, logarithm=TRUE)$modulus; c10_1
  c10_2 <- determinant(C_u_2, logarithm=TRUE)$modulus; c10_2
  c10_3 <- determinant(C_u_3, logarithm=TRUE)$modulus; c10_3
  c10_4 <- determinant(C_u_4, logarithm=TRUE)$modulus; c10_4
  
  
  
  
  # log likelihood: 
  ## log_like <- -1/2 *(c9 + c10) - (c1 + c2+ c4 + 2* c5) / 2
  log_like <- (-1/2 *(c9_1 + c10_1) - (c1_1 + c2_1 + c4_1 + 2* c5_1) / 2 
               -1/2 *(c9_2 + c10_2) - (c1_2 + c2_2 + c4_2 + 2* c5_2) / 2 
               -1/2 *(c9_3 + c10_3) - (c1_3 + c2_3 + c4_3 + 2* c5_3) / 2 
               -1/2 *(c9_4 + c10_4) - (c1_4 + c2_4 + c4_4 + 2* c5_4) / 2 ) 
  
  # log prior on theta: 
  logprior <-  sum(dbeta(u_cand, 2, 2, log = TRUE))  # beta(2, 2) prior
  #logprior <-  sum(dunif(theta0, 0, 1, log = TRUE))  # uniform(0, 1) prior
  # prior option 2: weak wrong prior (far away fro the truth)
  # logprior <-   (dbeta(u_cand[1], 1.5, 2.5, log = TRUE) + dbeta(u_cand[2], 1.5, 2.5, log = TRUE) + 
  #                  dbeta(u_cand[3], 2.5, 1.5, log = TRUE) + dbeta(u_cand[4], 2.5, 1.5, log = TRUE))  
  
  # prior option 3: prescient prior
  logprior <-   (dbeta(u_cand[1], 2.5, 1.5, log = TRUE) + dbeta(u_cand[2], 2.5, 1.5, log = TRUE) + 
                   dbeta(u_cand[3], 1.5, 2.5, log = TRUE) + dbeta(u_cand[4], 1.5, 2.5, log = TRUE))  
  
  
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
  
  par(mfrow = c(5, 1))
  plot(Us[1:iter, 1], ylim=c(0, 1), type="l", main="Trace plot of u1")
  plot(Us[1:iter, 2], ylim=c(0, 1), type="l", main="Trace plot of u2")
  plot(Us[1:iter, 3], ylim=c(0, 1), type="l", main="Trace plot of u3")
  plot(Us[1:iter, 4], ylim=c(0, 1), type="l", main="Trace plot of u4")
  plot(logposts[1:iter], type="l", main="Trace plot of log posterior")
  
  if (iter %% 5000) {
    save(list=c("mle_b_1", "mle_b_2","mle_b_3","mle_b_4",
                "Us",  "logposts", 
                "Tau2hat_1", "Tau2hat_2", "Tau2hat_3", "Tau2hat_4",
                "tau2hat_pc_1", "tau2hat_pc_2", "tau2hat_pc_3", "tau2hat_pc_4", 
                "urange", "steps"), file=paste("posterior_U_betaprior_all.RData", sep=""))
  }
}




###############################################################
############     Collect samples          #####################
###############################################################

## In scaled ranges: 

#x11()
par(mfrow = c(4, 1))
plot(Us[5000:iter, 1], ylim=c(0, 1), type="l", main="Trace plot of ff_Ns")
plot(Us[5000:iter, 2], ylim=c(0, 1), type="l", main="Trace plot of ff_Ms")
plot(Us[5000:iter, 3], ylim=c(0, 1), type="l", main="Trace plot of ff_Nr")
plot(Us[5000:iter, 4], ylim=c(0, 1), type="l", main="Trace plot of ff_Mr")
#dev.copy2pdf(file = "traceplot1.pdf")


#x11()
par(mfrow = c(4, 1))
plot(Us[, 1], ylim=c(0, 1), type="l", main="Trace plot of ff_Ns")
plot(Us[, 2], ylim=c(0, 1), type="l", main="Trace plot of ff_Ms")
plot(Us[, 3], ylim=c(0, 1), type="l", main="Trace plot of ff_Nr")
plot(Us[, 4], ylim=c(0, 1), type="l", main="Trace plot of ff_Mr")
#dev.copy2pdf(file = "traceplot2.pdf")


#x11()
par(mfrow = c(1, 4))
hist(Us[ , 1], xlim=c(0,1), main="Histogram of ff_Ns")
hist(Us[ , 2], xlim=c(0,1), main="Histogram of ff_Ms")
hist(Us[ , 3], xlim=c(0,1), main="Histogram of ff_Nr")
hist(Us[ , 4], xlim=c(0,1), main="Histogram of ff_Mr")
#dev.copy2pdf(file = "hist1.pdf")


par(mfrow = c(1, 1))
plot(logposts, type="l", main="Sample path of log posterior")

map <- c(Us[which.max(logposts),], max(logposts)); map

## Check the mixing:
library('coda')
mcmc.trace <- mcmc(Us)
summary(mcmc.trace)
acceptanceRate <- 1 - rejectionRate(mcmc.trace)
acceptanceRate
effectiveSize(mcmc.trace)

########################
## In original scales:  
########################

par(mfrow = c(1, 1))
plot(logposts, type="l", main="Sample path of log posterior")


# pairwise plots:
sols  <- cbind(Us, logposts)
unames <- c("ffactor_Ns", "ffactor_Ms","ffactor_Nr", "ffactor_Mr")
colnames(sols) <- c(unames, "logposts")

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

pairs(sols[,1:4], diag.panel = panel.hist, main="Estimated pairwise posterior distribution of U")


sols1 <- data.frame(sols[order(sols[, 5]),], rank=1:(nrow(sols))) #
cs <- terrain.colors(20)[21-(as.numeric(cut(rank(sols1[, 5]), breaks = 20)))]
cs[nrow(sols1)] <- "#000000" # plot in a reverse order of rows for not overlapping the optimal value
pairs(sols1[, 1:4], col=cs, pch=c(3, rep(20, nrow(sols1)-1)), main="Estimated pairwise posterior distribution of U")


## save(list=c("MLEs_Kdir_pc", "mle_b", "Us", "u0", "map", "logposts", "Tau2hat", "taub2hat", "urange", "steps"), file=paste("RData/pc_Kdir_cali_Bayes_betaprior.RData", sep=""))
save(list=c("MLEs_Kdir_pc", "MLEs_Kcross_pc", "MLEs_cdir_pc", "MLEs_ccross_pc", 
            "mle_b_1", "mle_b_2","mle_b_3","mle_b_4",
            "Us",  "logposts", 
            "Tau2hat_1", "Tau2hat_2", "Tau2hat_3", "Tau2hat_4",
            "tau2hat_pc_1", "tau2hat_pc_2", "tau2hat_pc_3", "tau2hat_pc_4", 
##             "urange", "steps"), file=paste("RData/all_outputs_posterior_U_weakwrong.RData", sep=""))
            "urange", "steps"), file=paste("RData/all_outputs_posterior_U_prescient.RData", sep=""))




