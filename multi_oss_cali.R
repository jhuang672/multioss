
########################################################################################################
## Optimization calibration using on-site surrogates on all 16 outputs together
## 16 -- PC --> 4 ----> 1
## Combining frequencies into PCs for 4 outputs
## On-site surrogates built from "multi_oss_emu.R"
## pointwise optimization calibration: using a zero-mean Gaussian process to model the discrepency
########################################################################################################


## Set up: 
library(laGP)        # for discrepancy term GP
library(plgp)        # for covariance matrix
library(nloptr)      # for optimization with optim
library(factoextra)  # for PCA
library(corrplot)    # for PCA visual 
library(ggpubr)      # grab PCA plots 
library(lhs)



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
Dfield$Cxx_9 <- Dfield$Cxx_9 * 10^3
Dfield$Cxx_11 <- Dfield$Cxx_11 * 10^3

Dfield$c_xy_2 <- Dfield$c_xy_2 * 10^3
Dfield$c_xy_5 <- Dfield$c_xy_5 * 10^3
Dfield$c_xy_9 <- Dfield$c_xy_9 * 10^3
Dfield$c_xy_11 <- Dfield$c_xy_11 * 10^3

summary(Dfield)


# load all isotseal data: 
# all inputs with 4 outputs at freq = 2: 
D_freq_2 <- read.csv("../../data/U_LHS_1000_2.csv") # 1000 maximinLHS design for all the same theta at each Xfield locations

# all outputs with 4 outputs at freqs 5, 7, 9, 11, 13: 
D_freq_5_13 <- read.csv("../../data/Multi_outputs_5_to_13_OSS_U_maximinLHS_1000.cvs") # 1000 maximinLHS design for all the same theta at each Xfield locations


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
  list_na_index[[i]] <- which(is.na(ym_1_2[1:1000 + 1000 * (i-1)]) == FALSE)
}   


# Scale both inputs X and U: 
XUori <- Disot[, 1:17]
Xori <- Dfield[, xnames]
Ulimits <- rbind(lower_limits, upper_limits)
colnames(Ulimits) <- unames
limits <- cbind(apply(Dfield[c(xnames)], 2, range), Ulimits)
XU <- matrix(NA, nrow = nrow(XUori), ncol = ncol(XUori))


# scale XU for isotseal: only 1000 runs need to be scaled 
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


# scale X for field data
X <- matrix(NA, nrow = nrow(Xori), ncol = ncol(Xori))
for (j in 1:ncol(Xori)) {
  X[ , j] <- (Xori[ , j] - limits[1, j]) / (limits[2, j] - limits[1, j]) 
}



## this part only needs to be pre-calculated once: 

# each of 292 locations have a different length due to a different rate of missing

## create CY() for all frequencies: kdir
CY_1_2 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1
CY_1_5 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1
CY_1_9 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1
CY_1_11 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1

## create CY() for all frequencies: kcross
CY_2_2 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1
CY_2_5 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1
CY_2_9 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1
CY_2_11 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1


## create CY() for all frequencies: cdir
CY_3_2 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1
CY_3_5 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1
CY_3_9 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1
CY_3_11 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1


## create CY() for all frequencies: 
CY_4_2 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1
CY_4_5 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1
CY_4_9 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1
CY_4_11 <- list()  # create a list of pre-calculated values, most of them are 1000 by 1


for (i in 1:292) {
  
  ## Fixed covariance part: 
  
  ## kdir
  Cm_1_2  <- covar.sep(U[list_na_index[[i]], ],  d = OSS_1$oss2[i, 1:4], g = OSS_1$oss2[i, 5]) # 1000 by 1000
  Cm_1_5  <- covar.sep(U[list_na_index[[i]], ],  d = OSS_1$oss5[i, 1:4], g = OSS_1$oss5[i, 5]) # 1000 by 1000
  Cm_1_9  <- covar.sep(U[list_na_index[[i]], ],  d = OSS_1$oss9[i, 1:4], g = OSS_1$oss9[i, 5]) # 1000 by 1000
  Cm_1_11 <- covar.sep(U[list_na_index[[i]], ],  d = OSS_1$oss11[i, 1:4], g = OSS_1$oss11[i, 5]) # 1000 by 1000
  
  ## kcross
  Cm_2_2  <- covar.sep(U[list_na_index[[i]], ],  d = OSS_2$oss2[i, 1:4], g = OSS_2$oss2[i, 5]) # 1000 by 1000
  Cm_2_5  <- covar.sep(U[list_na_index[[i]], ],  d = OSS_2$oss5[i, 1:4], g = OSS_2$oss5[i, 5]) # 1000 by 1000
  Cm_2_9  <- covar.sep(U[list_na_index[[i]], ],  d = OSS_2$oss9[i, 1:4], g = OSS_2$oss9[i, 5]) # 1000 by 1000
  Cm_2_11 <- covar.sep(U[list_na_index[[i]], ],  d = OSS_2$oss11[i, 1:4], g = OSS_2$oss11[i, 5]) # 1000 by 1000
  
  ## cdir
  Cm_3_2  <- covar.sep(U[list_na_index[[i]], ],  d = OSS_3$oss2[i, 1:4], g = OSS_3$oss2[i, 5]) # 1000 by 1000
  Cm_3_5  <- covar.sep(U[list_na_index[[i]], ],  d = OSS_3$oss5[i, 1:4], g = OSS_3$oss5[i, 5]) # 1000 by 1000
  Cm_3_9  <- covar.sep(U[list_na_index[[i]], ],  d = OSS_3$oss9[i, 1:4], g = OSS_3$oss9[i, 5]) # 1000 by 1000
  Cm_3_11 <- covar.sep(U[list_na_index[[i]], ],  d = OSS_3$oss11[i, 1:4], g = OSS_3$oss11[i, 5]) # 1000 by 1000
  
  
  ## Ccross
  Cm_4_2  <- covar.sep(U[list_na_index[[i]], ],  d = OSS_4$oss2[i, 1:4],  g = OSS_4$oss2[i, 5]) # 1000 by 1000
  Cm_4_5  <- covar.sep(U[list_na_index[[i]], ],  d = OSS_4$oss5[i, 1:4],  g = OSS_4$oss5[i, 5]) # 1000 by 1000
  Cm_4_9  <- covar.sep(U[list_na_index[[i]], ],  d = OSS_4$oss9[i, 1:4],  g = OSS_4$oss9[i, 5]) # 1000 by 1000
  Cm_4_11 <- covar.sep(U[list_na_index[[i]], ],  d = OSS_4$oss11[i, 1:4], g = OSS_4$oss11[i, 5]) # 1000 by 1000
  
  
  ## Inverse for 16 covariances: 
  
  ## Kdir 
  Cm_inv_1_2 <- solve(Cm_1_2)
  Cm_inv_1_5 <- solve(Cm_1_5)
  Cm_inv_1_9 <- solve(Cm_1_9)
  Cm_inv_1_11 <- solve(Cm_1_11)
  
  ## Kcross
  Cm_inv_2_2 <- solve(Cm_2_2)
  Cm_inv_2_5 <- solve(Cm_2_5)
  Cm_inv_2_9 <- solve(Cm_2_9)
  Cm_inv_2_11 <- solve(Cm_2_11)
  
  
  ## Cdir
  Cm_inv_3_2 <- solve(Cm_3_2)
  Cm_inv_3_5 <- solve(Cm_3_5)
  Cm_inv_3_9 <- solve(Cm_3_9)
  Cm_inv_3_11 <- solve(Cm_3_11)
  
  ## Ccross
  Cm_inv_4_2 <- solve(Cm_4_2)
  Cm_inv_4_5 <- solve(Cm_4_5)
  Cm_inv_4_9 <- solve(Cm_4_9)
  Cm_inv_4_11 <- solve(Cm_4_11)
  
  
  ## The fixed whole part: solve(Sigma(U, U)) %*% ym
  
  ## Kdir
  CY_1_2[[i]] <- Cm_inv_1_2 %*% ym_1_2[1:1000 + 1000 * (i-1)][list_na_index[[i]]]     # 1000 by 1 for all 292 runs
  CY_1_5[[i]] <- Cm_inv_1_5 %*% ym_1_5[1:1000 + 1000 * (i-1)][list_na_index[[i]]]     # 1000 by 1 for all 292 runs
  CY_1_9[[i]] <- Cm_inv_1_9 %*% ym_1_9[1:1000 + 1000 * (i-1)][list_na_index[[i]]]     # 1000 by 1 for all 292 runs
  CY_1_11[[i]] <- Cm_inv_1_11 %*% ym_1_11[1:1000 + 1000 * (i-1)][list_na_index[[i]]]     # 1000 by 1 for all 292 runs
  
  ## Kcross
  CY_2_2[[i]] <- Cm_inv_2_2 %*% ym_2_2[1:1000 + 1000 * (i-1)][list_na_index[[i]]]     # 1000 by 1 for all 292 runs
  CY_2_5[[i]] <- Cm_inv_2_5 %*% ym_2_5[1:1000 + 1000 * (i-1)][list_na_index[[i]]]     # 1000 by 1 for all 292 runs
  CY_2_9[[i]] <- Cm_inv_2_9 %*% ym_2_9[1:1000 + 1000 * (i-1)][list_na_index[[i]]]     # 1000 by 1 for all 292 runs
  CY_2_11[[i]] <- Cm_inv_2_11 %*% ym_2_11[1:1000 + 1000 * (i-1)][list_na_index[[i]]]     # 1000 by 1 for all 292 runs
  
  
  ## Cdir
  CY_3_2[[i]] <- Cm_inv_3_2 %*% ym_3_2[1:1000 + 1000 * (i-1)][list_na_index[[i]]]     # 1000 by 1 for all 292 runs
  CY_3_5[[i]] <- Cm_inv_3_5 %*% ym_3_5[1:1000 + 1000 * (i-1)][list_na_index[[i]]]     # 1000 by 1 for all 292 runs
  CY_3_9[[i]] <- Cm_inv_3_9 %*% ym_3_9[1:1000 + 1000 * (i-1)][list_na_index[[i]]]     # 1000 by 1 for all 292 runs
  CY_3_11[[i]] <- Cm_inv_3_11 %*% ym_3_11[1:1000 + 1000 * (i-1)][list_na_index[[i]]]     # 1000 by 1 for all 292 runs
  
  
  ## Ccross
  CY_4_2[[i]] <- Cm_inv_4_2 %*% ym_4_2[1:1000 + 1000 * (i-1)][list_na_index[[i]]]     # 1000 by 1 for all 292 runs
  CY_4_5[[i]] <- Cm_inv_4_5 %*% ym_4_5[1:1000 + 1000 * (i-1)][list_na_index[[i]]]     # 1000 by 1 for all 292 runs
  CY_4_9[[i]] <- Cm_inv_4_9 %*% ym_4_9[1:1000 + 1000 * (i-1)][list_na_index[[i]]]     # 1000 by 1 for all 292 runs
  CY_4_11[[i]] <- Cm_inv_4_11 %*% ym_4_11[1:1000 + 1000 * (i-1)][list_na_index[[i]]]     # 1000 by 1 for all 292 runs
  
}


############################################################################
## Option 3: fit one big principal component analysis on the 
##           discrepancies first, then explore different U values
############################################################################

## For a fair comparison on the likelihood: 
## Fit a big principal component on the observed discrepancy  ( YF - YM ) on all U values first! 

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


### fit PCA on observed biases: 
Bias_1 <- YF_1_rep - YM_1
Bias_2 <- YF_2_rep - YM_2
Bias_3 <- YF_3_rep - YM_3
Bias_4 <- YF_4_rep - YM_4

colnames(Bias_1) <- c("y2", "y5", "y9", "y11")
colnames(Bias_2) <- c("y2", "y5", "y9", "y11")
colnames(Bias_3) <- c("y2", "y5", "y9", "y11")
colnames(Bias_4) <- c("y2", "y5", "y9", "y11")

summary(Bias_1)
summary(Bias_2)
summary(Bias_3)
summary(Bias_4)


# Remove missing: 
Bias_1 <- na.omit(Bias_1)
Bias_2 <- na.omit(Bias_2)
Bias_3 <- na.omit(Bias_3)
Bias_4 <- na.omit(Bias_4)

cor(Bias_1)
cor(Bias_2)
cor(Bias_3)
cor(Bias_4)


# Perform PCA to combine frequencies: Kdir
Bias_pr_1 <- prcomp(Bias_1, scale = TRUE)
Bias_pr_1
Bias_pr_1$sdev^2
Bias_pr_1$sdev^2 /sum(Bias_pr_1$sdev^2)
summary(Bias_pr_1)
fviz_eig(Bias_pr_1)

# Perform PCA to combine frequencies: Kcross
Bias_pr_2 <- prcomp(Bias_2, scale = TRUE)
Bias_pr_2
Bias_pr_2$sdev^2
Bias_pr_2$sdev^2 /sum(Bias_pr_2$sdev^2)
summary(Bias_pr_2)
fviz_eig(Bias_pr_2)

# Perform PCA to combine frequencies: Cdir
Bias_pr_3 <- prcomp(Bias_3, scale = TRUE)
Bias_pr_3
Bias_pr_3$sdev^2
Bias_pr_3$sdev^2 /sum(Bias_pr_3$sdev^2)
summary(Bias_pr_3)
fviz_eig(Bias_pr_3)


# Perform PCA to combine frequencies: Ccross
Bias_pr_4 <- prcomp(Bias_4, scale = TRUE)
Bias_pr_4
Bias_pr_4$sdev^2
Bias_pr_4$sdev^2 /sum(Bias_pr_4$sdev^2)
summary(Bias_pr_4)
fviz_eig(Bias_pr_4)



## ------------------------------------------------------------------------
# generate priors for discrepancy term: 
da <- darg(list(mle=TRUE, max=5, min=eps), X, samp.size = nrow(X)); da  # also have stronger regularization on bias term


## prior for nugget effects: need to be careful in principal component space -----------------------------

ga_1 <- garg(list(mle=TRUE, max=1), Bias_pr_1$x[,1]); ga_1
ga_2 <- garg(list(mle=TRUE, max=1), Bias_pr_2$x[,1]); ga_2
ga_3 <- garg(list(mle=TRUE, max=1), Bias_pr_3$x[,1]); ga_3
ga_4 <- garg(list(mle=TRUE, max=1), Bias_pr_4$x[,1]); ga_4



## ------------------------------------------------------------------------
# Random initial values
n1 <- 500  # n1 = 500  LHS design in space-filling initial values for U
sol <- data.frame(matrix(NA, nrow=n1, ncol=10))
names(sol) <- c( c("Ns_ini", "Ms_ini", "Nr_ini", "Mr_ini"), c(unames),  "value", "evals")


##  Define bias function ------------------------------------------------------------------------
##  Use the PC from multiple levels of frequency to learn discrepancy

bhat.fit <- function(X, Y_1, Y_2, Y_3, Y_4, 
                     YM_1, YM_2, YM_3, YM_4, da, 
                     ga_1, ga_2, ga_3, ga_4, Uprior, clean=TRUE)  { 
  
  ## create a PCA of discrepancies: 

  ## use a fixed global PCA on every new U values for fair comparison: 
  bias_u_1 <-  Y_1 - YM_1    ##  observed bias given U value
  bias_u_2 <-  Y_2 - YM_2    ##  observed bias given U value
  bias_u_3 <-  Y_3 - YM_3    ##  observed bias given U value
  bias_u_4 <-  Y_4 - YM_4    ##  observed bias given U value
  
  colnames(bias_u_1) <- c("y2", "y5", "y9", "y11")  ## name it for PCA rotation
  colnames(bias_u_2) <- c("y2", "y5", "y9", "y11")  ## name it for PCA rotation
  colnames(bias_u_3) <- c("y2", "y5", "y9", "y11")  ## name it for PCA rotation
  colnames(bias_u_4) <- c("y2", "y5", "y9", "y11")  ## name it for PCA rotation
  
  bias_pc_1 <- predict(Bias_pr_1, bias_u_1)             ## rotate bias_u using fixed PCA
  bias_pc_2 <- predict(Bias_pr_2, bias_u_2)             ## rotate bias_u using fixed PCA
  bias_pc_3 <- predict(Bias_pr_3, bias_u_3)             ## rotate bias_u using fixed PCA
  bias_pc_4 <- predict(Bias_pr_4, bias_u_4)             ## rotate bias_u using fixed PCA
  
  
  
  ## Fit GP on the first principal component of all biases 
  bhat_1 <- newGPsep(X, bias_pc_1[, 1], d=rep(da$start, ncol(X)), g=ga_1$start, dK=TRUE)
  bhat_2 <- newGPsep(X, bias_pc_2[, 1], d=rep(da$start, ncol(X)), g=ga_2$start, dK=TRUE)
  bhat_3 <- newGPsep(X, bias_pc_3[, 1], d=rep(da$start, ncol(X)), g=ga_3$start, dK=TRUE)
  bhat_4 <- newGPsep(X, bias_pc_4[, 1], d=rep(da$start, ncol(X)), g=ga_4$start, dK=TRUE)
  
  
  
  if(ga_1$mle) {
    cmle_1 <- jmleGPsep(bhat_1, drange=c(da$min, da$max), grange=c(ga_1$min, ga_1$max), dab=da$ab, gab=ga_1$ab)
    cmle_2 <- jmleGPsep(bhat_2, drange=c(da$min, da$max), grange=c(ga_2$min, ga_2$max), dab=da$ab, gab=ga_2$ab)
    cmle_3 <- jmleGPsep(bhat_3, drange=c(da$min, da$max), grange=c(ga_3$min, ga_3$max), dab=da$ab, gab=ga_3$ab)
    cmle_4 <- jmleGPsep(bhat_4, drange=c(da$min, da$max), grange=c(ga_4$min, ga_4$max), dab=da$ab, gab=ga_4$ab)
  }
  else {
    cmle_1 <- mleGPsep(bhat_1, tmin=rep(da$min, ncol(X)), tmax=rep(da$max, ncol(X)), ab=da$ab)
    cmle_2 <- mleGPsep(bhat_2, tmin=rep(da$min, ncol(X)), tmax=rep(da$max, ncol(X)), ab=da$ab)
    cmle_3 <- mleGPsep(bhat_3, tmin=rep(da$min, ncol(X)), tmax=rep(da$max, ncol(X)), ab=da$ab)
    cmle_4 <- mleGPsep(bhat_4, tmin=rep(da$min, ncol(X)), tmax=rep(da$max, ncol(X)), ab=da$ab)
  }
    
    ## Put together as independent:
    ## cmle$nll <- (- llikGPsep(bhat_1, dab=da$ab, gab=ga_1$ab) - llikGPsep(bhat_2, dab=da$ab, gab=ga_2$ab) 
    u_post <- (- llikGPsep(bhat_1, dab=da$ab, gab=ga_1$ab) - llikGPsep(bhat_2, dab=da$ab, gab=ga_2$ab) 
    - llikGPsep(bhat_3, dab=da$ab, gab=ga_3$ab) - llikGPsep(bhat_4, dab=da$ab, gab=ga_4$ab)  +  Uprior)

    if(clean) {
    deleteGPsep(bhat_1)
    deleteGPsep(bhat_2)
    deleteGPsep(bhat_3)
    deleteGPsep(bhat_4)
  }
  return(u_post)
}

## Define the calibration function: 

calib <- function(u, X, 
                  Y_1, Y_2, Y_3, Y_4,  
                  da, ga_1, ga_2, ga_3, ga_4, 
                  U, 
                  CY_1_2, CY_1_5, CY_1_9, CY_1_11, 
                  CY_2_2, CY_2_5, CY_2_9, CY_2_11, 
                  CY_3_2, CY_3_5, CY_3_9, CY_3_11, 
                  CY_4_2, CY_4_5, CY_4_9, CY_4_11, 
                  MLEs)
  #calib <- function(u, X, Y, da, U, CY_2, CY_5, CY_9, CY_11, MLEs)
{
  # only one u at a time in this case: 
  Ym_1_2 <- rep(NA, 292)
  Ym_1_5 <- rep(NA, 292)
  Ym_1_9 <- rep(NA, 292)
  Ym_1_11 <- rep(NA, 292)
  
  Ym_2_2 <- rep(NA, 292)
  Ym_2_5 <- rep(NA, 292)
  Ym_2_9 <- rep(NA, 292)
  Ym_2_11 <- rep(NA, 292)
  
  Ym_3_2 <- rep(NA, 292)
  Ym_3_5 <- rep(NA, 292)
  Ym_3_9 <- rep(NA, 292)
  Ym_3_11 <- rep(NA, 292)
  
  Ym_4_2 <- rep(NA, 292)
  Ym_4_5 <- rep(NA, 292)
  Ym_4_9 <- rep(NA, 292)
  Ym_4_11 <- rep(NA, 292)
  
  
  # loop over 292 Xfield
  for (i in 1:292){
    
    ## Frequency 2: 
    CX_1_2 <- covar.sep(u, U[list_na_index[[i]], ], d = OSS_1$oss2[i, 1:4], g = OSS_1$oss2[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
    Ym_1_2[i] <- CX_1_2 %*% CY_1_2[[i]]  # predict from ym for means
    
    CX_2_2 <- covar.sep(u, U[list_na_index[[i]], ], d = OSS_2$oss2[i, 1:4], g = OSS_2$oss2[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
    Ym_2_2[i] <- CX_2_2 %*% CY_2_2[[i]]  # predict from ym for means
    
    CX_3_2 <- covar.sep(u, U[list_na_index[[i]], ], d = OSS_3$oss2[i, 1:4], g = OSS_3$oss2[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
    Ym_3_2[i] <- CX_3_2 %*% CY_3_2[[i]]  # predict from ym for means
    
    CX_4_2 <- covar.sep(u, U[list_na_index[[i]], ], d = OSS_4$oss2[i, 1:4], g = OSS_4$oss2[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
    Ym_4_2[i] <- CX_4_2 %*% CY_4_2[[i]]  # predict from ym for means
    

    ## Frequency 5: 
    CX_1_5 <- covar.sep(u, U[list_na_index[[i]], ], d = OSS_1$oss5[i, 1:4], g = OSS_1$oss5[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
    Ym_1_5[i] <- CX_1_5 %*% CY_1_5[[i]]  # predict from ym for means
    
    CX_2_5 <- covar.sep(u, U[list_na_index[[i]], ], d = OSS_2$oss5[i, 1:4], g = OSS_2$oss5[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
    Ym_2_5[i] <- CX_2_5 %*% CY_2_5[[i]]  # predict from ym for means
    
    CX_3_5 <- covar.sep(u, U[list_na_index[[i]], ], d = OSS_3$oss5[i, 1:4], g = OSS_3$oss5[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
    Ym_3_5[i] <- CX_3_5 %*% CY_3_5[[i]]  # predict from ym for means
    
    CX_4_5 <- covar.sep(u, U[list_na_index[[i]], ], d = OSS_4$oss5[i, 1:4], g = OSS_4$oss5[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
    Ym_4_5[i] <- CX_4_5 %*% CY_4_5[[i]]  # predict from ym for means
    
    
    
    ## Frequency 9: 
    CX_1_9 <- covar.sep(u, U[list_na_index[[i]], ], d = OSS_1$oss9[i, 1:4], g = OSS_1$oss9[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
    Ym_1_9[i] <- CX_1_9 %*% CY_1_9[[i]]  # predict from ym for means
    
    CX_2_9 <- covar.sep(u, U[list_na_index[[i]], ], d = OSS_2$oss9[i, 1:4], g = OSS_2$oss9[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
    Ym_2_9[i] <- CX_2_9 %*% CY_2_9[[i]]  # predict from ym for means
    
    CX_3_9 <- covar.sep(u, U[list_na_index[[i]], ], d = OSS_3$oss9[i, 1:4], g = OSS_3$oss9[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
    Ym_3_9[i] <- CX_3_9 %*% CY_3_9[[i]]  # predict from ym for means
    
    CX_4_9 <- covar.sep(u, U[list_na_index[[i]], ], d = OSS_4$oss9[i, 1:4], g = OSS_4$oss9[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
    Ym_4_9[i] <- CX_4_9 %*% CY_4_9[[i]]  # predict from ym for means
    
    
    ## Frequency 11: 
    CX_1_11 <- covar.sep(u, U[list_na_index[[i]], ], d = OSS_1$oss11[i, 1:4], g = OSS_1$oss11[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
    Ym_1_11[i] <- CX_1_11 %*% CY_1_11[[i]]  # predict from ym for means
    
    CX_2_11 <- covar.sep(u, U[list_na_index[[i]], ], d = OSS_2$oss11[i, 1:4], g = OSS_2$oss11[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
    Ym_2_11[i] <- CX_2_11 %*% CY_2_11[[i]]  # predict from ym for means
    
    CX_3_11 <- covar.sep(u, U[list_na_index[[i]], ], d = OSS_3$oss11[i, 1:4], g = OSS_3$oss11[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
    Ym_3_11[i] <- CX_3_11 %*% CY_3_11[[i]]  # predict from ym for means
    
    CX_4_11 <- covar.sep(u, U[list_na_index[[i]], ], d = OSS_4$oss11[i, 1:4], g = OSS_4$oss11[i, 5]) # here g=0 for out-of-sample prediction using Kronecker delta
    Ym_4_11[i] <- CX_4_11 %*% CY_4_11[[i]]  # predict from ym for means
  
  }
  
  
  YM_1 <- data.frame(Ym_1_2, Ym_1_5, Ym_1_9, Ym_1_11)
  YM_2 <- data.frame(Ym_2_2, Ym_2_5, Ym_2_9, Ym_2_11)
  YM_3 <- data.frame(Ym_3_2, Ym_3_5, Ym_3_9, Ym_3_11)
  YM_4 <- data.frame(Ym_4_2, Ym_4_5, Ym_4_9, Ym_4_11)
  
  ## Prior option 1: regularize with a marginal beta(2, 2) prior on U:
  ## Uprior <-  - sum(dbeta(u, 2, 2, log = TRUE))  # negative log likelihood to minimize
  
  ## Prior option 2: pick prior to be far away from the solutiion 
  ## Uprior <-  - sum(dbeta(u, 2, 2, log = TRUE))  # negative log likelihood to minimize
  Uprior <-  - (dbeta(u[1], 1.5, 2.5, log = TRUE) + dbeta(u[2], 1.5, 2.5, log = TRUE) + 
                  dbeta(u[3], 2.5, 1.5, log = TRUE) + dbeta(u[4], 2.5, 1.5, log = TRUE))  # negative log likelihood to minimize
  
  
  # Uniform [0, 1] prior on U:   
  #cmle <- bhat.fit(X, Y, Ym, da, ga, Uprior=0, clean=TRUE) # fit GP for bias and calculate likelihood for this u (takes time)
  
  # Beta prior on U:
  cmle <- bhat.fit(X,  Y_1, Y_2, Y_3, Y_4, 
                   YM_1, YM_2, YM_3, YM_4, da, 
                   ga_1, ga_2, ga_3, ga_4, Uprior, clean=TRUE) # fit GP for bias and calculate likelihood for this u (takes time)
  # cmle <- bhat.fit(X, Y, YM, da, Uprior, clean=TRUE) # fit GP for bias and calculate likelihood for this u (takes time)
  return(cmle)
}



# define objective function:
obj <- function(u, X,  Y_1, Y_2, Y_3, Y_4, 
                da, 
                ga_1, ga_2, ga_3, ga_4,  U, 
                CY_1_2, CY_1_5, CY_1_9, CY_1_11, 
                CY_2_2, CY_2_5, CY_2_9, CY_2_11, 
                CY_3_2, CY_3_5, CY_3_9, CY_3_11, 
                CY_4_2, CY_4_5, CY_4_9, CY_4_11, 
                MLEs) calib(u, X, Y_1, Y_2, Y_3, Y_4, da, 
        ga_1, ga_2, ga_3, ga_4, U, 
        CY_1_2, CY_1_5, CY_1_9, CY_1_11, 
        CY_2_2, CY_2_5, CY_2_9, CY_2_11, 
        CY_3_2, CY_3_5, CY_3_9, CY_3_11, 
        CY_4_2, CY_4_5, CY_4_9, CY_4_11,
        MLEs)
##        MLEs)$nll

## calibrate with different initial values for U: 

for(r in 1:n1) {
  ptm <- proc.time()
  opt <- nloptr(x0 = as.numeric(U[r, ]), eval_f = obj, lb=rep(eps, 4), ub=rep(1-eps, 4), # log(0) from beta prior gives -Inf
                opts = list("algorithm"="NLOPT_LN_COBYLA", maxeval = 1000, xtol_rel=eps),  # restart after 1000 runs if a starting value is bad
                X = X, 
                Y_1 = YF_1, Y_2 = YF_2, Y_3 = YF_3, Y_4 = YF_4, 
                da = da, 
                ga_1 = ga_1, ga_2 = ga_2, ga_3 = ga_3, ga_4 = ga_4, 
                U = U,
                CY_1_2 = CY_1_2, CY_1_5 = CY_1_5, CY_1_9 = CY_1_9,  CY_1_11 = CY_1_11, 
                CY_2_2 = CY_2_2, CY_2_5 = CY_2_5, CY_2_9 = CY_2_9,  CY_2_11 = CY_2_11, 
                CY_3_2 = CY_3_2, CY_3_5 = CY_3_5, CY_3_9 = CY_3_9,  CY_3_11 = CY_3_11, 
                CY_4_2 = CY_4_2, CY_4_5 = CY_4_5, CY_4_9 = CY_4_9,  CY_4_11 = CY_4_11, 
                MLEs = MLEs)    ## be careful about the priors

  
  ## save results: 
  sol[r,] <- c(as.numeric(XU[r, 14:17]), opt$solution, opt$objective, opt$iterations)
  
  ## print out current run: 
  print(r)
  print(proc.time() - ptm)
  print(sol[r, ])
  print(opt)
  save(list=c("sol"), file="sols_all_outputs.RData")  # save sols after each iteration
  pairs(sol[, 5:8])
  
  # Save a plot after each 5 iterations: 
  if (r %% 5){
  x11()
  sols <- na.omit(sol)
  sols <- sols[order(sols$value),]
  sols1 <- data.frame(sols[order(sols$value, decreasing=TRUE), ], rank=1:(nrow(sols))); #
  
  cs <- terrain.colors(100)[(as.numeric(cut(rank(sols1$value), breaks = 100)))]
  cs[r] <- "#000000" # plot in a reverse order of rows for not overlapping the optimal value
  pairs(sols1[, 5:8], col=cs, pch=c(rep(20, r-1), 3), main="Pairwise scatterplots of local optima(heat for likelihood)")
  
  dev.copy2pdf(file = paste("u_hat_optim_all_outputs.pdf", sep=""))
  dev.off()
  }
}


################################ End of loops #####################################
###################################################################################




head(sol, n=50)

solu <- sol[which.min(sol$value),]; solu
uhat <- solu[5:8]; uhat



# save 
save(list=c("sol"), file=paste("RData/all_outputs_cali_opti_betaprior", ".RData", sep=""))



# visualize in terminal:
sols <- sol
sols <- sols[order(sols$value),];head(sols, n=10)
sols1 <- data.frame(sols[order(sols$value, decreasing=TRUE), ], rank=1:(nrow(sol))); #

head(sols1, n=10)
cs <- heat.colors(100)[(as.numeric(cut(rank(sols1$value), breaks = 100)))]
##cs[nrow(sols)] <- "#000000" # plot in a reverse order of rows for not overlapping the optimal value
##pairs(sols1[, 5:8], col=cs, pch=c(rep(20, nrow(sols)-1), 3), main="Marginal scatterplots of local optima (heat for posterior)")
pairs(sols1[, 5:8], col=cs, pch=c(rep(20, nrow(sols))), main="Pairwise scatterplots of local optima(heat for likelihood)")

