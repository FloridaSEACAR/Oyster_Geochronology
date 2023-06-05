# 20th March 2023
# A test version to calibrate against a a bomb curve
# Code to implement NP Bayesian Calibration and Summaristion of Related 14C Determinations

set.seed(7)

# Read in postbomb curve
# jch made some edits
sampled_calcurve <- read.table("Curves/bombBahamasto10000calBP.14c", sep = "\t", header=TRUE)
sampled_calcurve <- sampled_calcurve[,1:5]
names(sampled_calcurve) <- c("Year", "D14C", "D14Csd", "F14C", "F14Csd")
PresentDay <- max(sampled_calcurve$Year) + 1 






# Interpolate curve onto regular 1/3 year intervals
calcurve <- data.frame(Year = seq(min(sampled_calcurve$Year), max(sampled_calcurve$Year), by = 1))

calcurve$calage <- PresentDay - calcurve$Year 
  
calcurve$c14age <- approx(sampled_calcurve$Year, sampled_calcurve$F14C, 
                          xout = seq(min(calcurve$Year), max(calcurve$Year), by = 1))$y
calcurve$c14sig <- approx(sampled_calcurve$Year, sampled_calcurve$F14Csd, 
                          xout = seq(min(calcurve$Year), max(calcurve$Year), by = 1))$y

# Test to check that interpolation gives the same results
# plot(sampled_calcurve$Year, sampled_calcurve$F14C, type = "l")
# lines(calcurve$Year, calcurve$c14age, col = "red", lty = 2)


# Read in the necessary functions
source('WalkerDirichletMixtureUpdateFunsFinal.R') # This also reads in the slice sampling SliceUpdateFuns.R
source("WalkerMasterFunctionFinal.R")
source("SimStudyFuncsFinal.R")

# Read in data
# x - c14ages
# xsig - corresponding 1 sigma 
Data <- read.csv("Curves/LeodiaAges.csv", header = TRUE, sep =  ",")
x <- Data[,2]
xsig <- Data[,3]

#############################################################################
# Now choose hyperparameters
############################################################
# Prior on the concentration parameter
# Place  a gamma prior on alpha
# alpha ~ Gamma(alphaprshape, alphaprrate) 
# A small alpha means more concentrated (i.e. few clusters)
# Large alpha not concentrated (many clusters)
cprshape <- alphaprshape <- 1 
cprrate <- alphaprrate <- 1 

#### Updated adaptive version
# Perform initial calibration of x and xsig
initprobs <- mapply(calibind, x, xsig, MoreArgs = list(calmu = calcurve$c14age, calsig = calcurve$c14sig))
# Find the SPD estimate (save as dataframe)
SPD <- data.frame(calage = calcurve$calage,
                  prob = apply(initprobs, 1, sum)/dim(initprobs)[2])
SPDcol <- grey(0.1, alpha = 0.5)



plot(calcurve$calage, calcurve$c14age, col = "blue", xlim = range(calcurve$calage),
     xlab = "Year (before present day)", ylab = expression(paste("F"^14, "C")),
     type = "l")
rug(x, side = 2)

par(new = TRUE, las = 1)
plot(calcurve$calage, SPD$prob, lty = 1, col = "purple", 
     ylim = c(0, 2*max(SPD$prob)), xlim = range(calcurve$calage),
     axes = FALSE, xlab = NA, ylab = NA, type = "n")
polygon(c(SPD$calage, rev(SPD$calage)), c(SPD$prob, rep(0, length(SPD$prob))), 
        border = NA, col = SPDcol)


inittheta <- calcurve$calage[apply(initprobs, 2, which.max)]
# Choose A and B from range of theta
A <- median(inittheta)
B <- 1 / (max(inittheta) - min(inittheta))^2
maxrange <- max(inittheta) - min(inittheta)

# Parameters for sigma2 (sigma^2 ~ InvGamma(nu1, nu2))
# E[tau] = (1/100)^2 Var[tau] = (1/100)^4
# Interval for sigma2 is approx 1/ c(nu2/nu1 - 2*nu2^2/nu1, nu2/nu1 + 2*nu2^2/nu1)
tempspread <- 0.1 * mad(inittheta)
tempprec <- 1/(tempspread)^2
nu1 <- 0.25
nu2 <- nu1 / tempprec

# Setup the NP method
lambda <- (100/maxrange)^2  # Each muclust ~ N(mutheta, sigma2/lambda)


# Choose number of iterations for sampler
niter <- 100000 
nthin <- 10 # Suggest not going higher than 5 if you want the histogram of posterior not to be sparse   

WalkerTemp <- WalkerBivarDirichlet(x = x, xsig = xsig, 
                                   lambda = lambda, nu1 = nu1, nu2 = nu2, 
                                   A = A, B = B, 
                                   cprshape = cprshape, cprrate = cprrate, 
                                   niter = niter, nthin = nthin, theta = inittheta,
                                   slicew = 70, m = 10, calcurve = calcurve, kstar = 8)

# To plot the predictive distribution then you run 
source("PostBombWalkerProcessingFinal.R")






# To access the posterior calendar age estimate for individual determination then you can look at:
# WalkerTemp$theta[,10] # MCMC chain for 10th determination (will need to remove burn in)

# If we want to plot e.g. the posterior calendar age density against the curve then we can run the below
# ident is the determination you want to calibrate
plotindpost(WalkerTemp, ident = , y = x, er = xsig, 
            calcurve = calcurve, PresentDay = PresentDay) 
