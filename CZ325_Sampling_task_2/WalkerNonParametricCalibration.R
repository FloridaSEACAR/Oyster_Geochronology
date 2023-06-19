# function to encapsulate the WalkerNP callibration

# DEP AGREEMENT NO. CZ325 Deliverable 2b (Radiocarbon sample size for age and time-averaging determination)
# Curve can only contain fraction maderns with columns
#  "Year", "D14C", "D14Csd", "F14C", "F14Csd"
WalkerNonParametricCalibration = function(Curve, Determinations){
require(dplyr)  
set.seed(7)
  sampled_calcurve = Curve
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
  source(here::here("CZ325_Sampling_task_2/WalkerDirichletMixtureUpdateFunsFinal.R")) # This also reads in the slice sampling SliceUpdateFuns.R
  source(here::here("CZ325_Sampling_task_2/WalkerMasterFunctionFinal.R"))
  source(here::here("CZ325_Sampling_task_2/SimStudyFuncsFinal.R"))
  
  # Read in data
  # x - c14ages
  # xsig - corresponding 1 sigma 
  Data <- Determinations
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
  
  # for us, this is TAV or the output of SUM in OxCal
  # Find the SPD estimate (save as dataframe)
  SPD <- data.frame(calage = calcurve$calage,
                    prob = apply(initprobs, 1, sum)/dim(initprobs)[2])
  
  calrange <- range(calcurve$calage)
  inittheta <- sample(x = calrange[1]:calrange[2], size = length(x))
  
  # Choose A and B from range of bombcurve rather than inittheta
  maxrange <- calrange[2] - calrange[1]
  A <- median(calcurve$calage)
  B <- 1 / (maxrange)^2
  
  # Parameters for sigma2 (sigma^2 ~ InvGamma(nu1, nu2))
  # E[tau] = (1/100)^2 Var[tau] = (1/100)^4
  # Interval for sigma2 is approx 1/ c(nu2/nu1 - 2*nu2^2/nu1, nu2/nu1 + 2*nu2^2/nu1)
  tempspread <- 0.1 * sd(calcurve$calage)
  tempprec <- 1/(tempspread)^2
  nu1 <- 0.25
  nu2 <- nu1 / tempprec
  
  # Setup the NP method
  lambda <- (100/maxrange)^2  # Each muclust ~ N(mutheta, sigma2/lambda)
  
  ####### END OF ALL PARAMETER CHANGES
  
  
  
  # Choose number of iterations for sampler
  niter <- 100000 
  nthin <- 10 # Suggest not going higher than 5 if you want the histogram of posterior not to be sparse   
  
  WalkerTemp <- WalkerBivarDirichlet(x = x, xsig = xsig, 
                                     lambda = lambda, nu1 = nu1, nu2 = nu2, 
                                     A = A, B = B, 
                                     cprshape = cprshape, cprrate = cprrate, 
                                     niter = niter, nthin = nthin, theta = inittheta,
                                     slicew = 70, m = 10, calcurve = calcurve, kstar = 8)
  
  
  Temp <- WalkerTemp
  
  npost <- dim(Temp$theta)[1]
  nburn <- floor(npost/2)
  npostsum <- 2000
  out <- Temp
  
  # Find predictive density
  npost <- dim(Temp$delta)[1]
  nburn <- floor(npost/2)
  npostsum <- 2000
  
  # Create the range over which to plot the density (note we will have to rescale the pdfs to account for calculation grid)
  ncalc <- length(calcurve$calage) 
  tempx <- seq(min(calcurve$calage), max(calcurve$calage), length = ncalc) 
  
  #Now choose the ids of the posterior sample
  sampid <- sample(x = nburn:npost, size = npostsum, replace = npostsum > (npost-nburn))
  
  # Now work out the actual posterior predictive density
  postDmat <- apply(as.row(sampid), 2, function(i, out, x, lambda, nu1, nu2) 
    WalkerFindpred(x, w = out$w[[i]], phi = out$phi[[i]], tau = out$tau[[i]], 
                   muphi = out$muphi[i], lambda = lambda, nu1 = nu1, nu2 = nu2), 
    out = Temp, x = tempx, lambda = lambda, nu1 = nu1, nu2 = nu2)
  
  # Rescale so it is a pdf
  postDmat <- postDmat  * (max(tempx) - min(tempx))/(ncalc) 
  # Find CI and mean
  postdenCI <- apply(postDmat, 1, quantile, probs = c(0.025, 0.975))
  postden <- apply(postDmat, 1, mean)
  
  
  
  return(data.frame(SPD) %>% mutate(NonParametric = rev(postden)))
  
}