
# DEP AGREEMENT NO. CZ325 Deliverable 1a:	Modeling approach to age and time-averaging estimation
# Helper function to make random draws from a posterior distribution


# start with a year, probability data frame

RandomEmpirical = function (n = 1, p)
{
  val = NULL
  p %>% mutate( cdf = cumsum(p$probability)) -> p
  for( ii in 1:n){
    y = runif(1)
    val = c(val, p$year[which.max(y <= p$cdf)])
  }
  return(val)
}

# assumes the distributions are subsetted
TAV = function(posterior_distribution){
  require(dplyr)
  specimen_sample = unique(posterior_distribution$name)
  n = length(specimen_sample)
  p =  posterior_distribution %>% group_by(year) %>% dplyr::summarize( probability = sum(probability)/n)
  return( p)
}

# doesn't assume the distributions are subsetted by sample
# TAV2 = function(posterior_distribution, samp_var = "igsn", year_var = "year", spec_var = "name", prob_var = "probability"){
#   require(data.table)
#   p = posterior_distribution[, .(probability = sum(eval(as.name(prob_var)))/length(unique(eval(as.name(spec_var))))), by = list(eval(as.name(samp_var)), eval(as.name(year_var)))]
#   return( p)
# }
TAV2 = function(posterior_distribution, samp_var = "igsn", year_var = "year", spec_var = "name", prob_var = "probability"){
  require(data.table)
  setDT(posterior_distribution)
  
  samp_var <- as.name(samp_var)
  year_var <- as.name(year_var)
  spec_var <- as.name(spec_var)
  prob_var <- as.name(prob_var)
  
  p = posterior_distribution[, .(probability = sum(eval(prob_var))/length(unique(eval(spec_var)))), by = list(eval(samp_var), eval(year_var))]
  setnames(p, c("samp_var", "year_var"), c(as.character(samp_var), as.character(year_var)))
  
  return( p)
}

Posterior_IQR = function(posterior_distribution, lo = 0.25, hi = 0.75)
{
  cummulative = cumsum(posterior_distribution$probability)
  year = posterior_distribution$year
  return(year[max(which(cummulative <= lo))] - year[max(which(cummulative <= hi))])
}

# 95% CR on TAV

CR_95_TAV = function(posterior_distribution){
  
  tmp = TAV(posterior_distribution)
  tav_iqr = Posterior_IQR(tmp, lo = 0.025, hi = 0.975)
  return(tav_iqr)
}






CPE = function(posterior_distribution){
  
  
  tmp = TAV2(posterior_distribution)
  tav_iqr = Posterior_IQR(tmp)
  
  specimen_sample = unique(posterior_distribution$name)
  # compute iqr's over 
  specimen_iqrs = numeric(0)
  
  for( gg in specimen_sample){
    specimen_iqrs = c(specimen_iqrs, Posterior_IQR(posterior_distribution %>% dplyr::filter(name == gg)))
  }
  
  aer = median(specimen_iqrs)
  return(list(tav_iqr = tav_iqr, cpe = tav_iqr - aer))
}


# length of highest density interval
HDI = function(posterior_distribution, credMass = 0.95){
  require(HDInterval)
  #hdi = hdi(density(RandomEmpirical(n = 1000, posterior_distribution)), allowSplit=TRUE, credMass = credMass)
  hdi = hdi(structure(list(x = posterior_distribution$year, y = posterior_distribution$probability, bw = 1, n = 1000, call = "posterior", data.name = "x", has.na = FALSE), class = "density"
            , names = c("x","y","bw", "n", "call","data.name","has.na")), allowSplit=TRUE, credMass = credMass)
  return(sum(hdi[,2]-hdi[,1])) # need the sum because there could be multiple intervals
}

# highest density interval
HDI_Interval = function(posterior_distribution, credMass = 0.95){
  require(HDInterval)
  #hdi = hdi(density(RandomEmpirical(n = 1000, posterior_distribution)), allowSplit=TRUE, credMass = credMass)
  hdi = hdi(structure(list(x = posterior_distribution$year, y = posterior_distribution$probability, bw = 1, n = 1000, call = "posterior", data.name = "x", has.na = FALSE), class = "density"
                      , names = c("x","y","bw", "n", "call","data.name","has.na")), allowSplit=TRUE, credMass = credMass)
  d = data.frame( )
  c = cumsum(posterior_distribution$probability)
  for( ii in 1:nrow(hdi))

    d = rbind(d, data.frame(begin = hdi[ii,1]
                            , end = hdi[ii,2]
                            , mass = c[posterior_distribution$year == hdi[ii,2]]-c[posterior_distribution$year == hdi[ii,1]]
                            , length = hdi[ii,2]-hdi[ii,1]+1
                            , density = (c[posterior_distribution$year == hdi[ii,2]]-c[posterior_distribution$year == hdi[ii,1]])/(hdi[ii,2]-hdi[ii,1]+1)))
  return(d) 
}

# zero out everything except highest density interval
HDI_Interval_Mask = function(posterior_distribution, credMass = 0.95){
  require(HDInterval)
  #hdi = hdi(density(RandomEmpirical(n = 1000, posterior_distribution)), allowSplit=TRUE, credMass = credMass)
  hdi = hdi(structure(list(x = posterior_distribution$year, y = posterior_distribution$probability, bw = 1, n = 1000, call = "posterior", data.name = "x", has.na = FALSE), class = "density"
                      , names = c("x","y","bw", "n", "call","data.name","has.na")), allowSplit=TRUE, credMass = credMass)
  d = data.frame( )
  p = posterior_distribution
  p$probability = 0
  pos = NULL
  for( ii in 1:nrow(hdi))
   {
    pos =c(pos, which(p$year >= hdi[ii,1] & p$year <= hdi[ii,2])) 
  }
  
   p$probability[pos] = posterior_distribution$probability[pos]
    
  return(p) 
}

# posterior_distribution contains a data frame of multiple posterior distributions
# this will produce an average of hdi's of each posterior
HDI_of_HDI = function(posterior_distribution,credMass = 0.5){
  n = length(unique(posterior_distribution$name))
  posterior_distribution %>% group_by(name) %>% reframe(HDI_Interval_Mask(data.frame(year, probability), credMass = credMass))  -> tmp
  tmp %>% group_by(year) %>% summarize( probability = sum(probability)) -> d
  d$probability = d$probability/n # average over names
  d$probability = d$probability/sum(d$probability)
  return(HDI(d, credMass = credMass)) # could use another credMass
}

#95% HDI on TAV

HDI_95_TAV = function(posterior_distribution, credMass = 0.95){
  
  tmp = TAV(posterior_distribution)
  return(HDI(tmp, credMass = 0.95))
}



# This is the "Olszewski" method, essentially taking the standard deviation of the posterior distributions averages
# This may not be meaningful for bimodal distributions lke the "post-bomb" posteriors found in Kowalewski et al 2018

ETA = function(posterior_distribution){
  group = unique(posterior_distribution$name)
  N = 10000
  sub = data.frame( x = numeric(), name = character())
  for( nn in group)
  {
    tmp = dplyr::filter(.data = posterior_distribution, name == nn)
    sub= rbind(sub,data.frame( x = RandomEmpirical(n = N,tmp), name = nn))
  }
  
  
  global_mean = mean(sub$x)
  SST = sum((sub$x-global_mean)^2)
  (sub %>% group_by(name) %>% summarize( SS = sum( (x - mean(x))^2)) %>% summarize( SSW = sum(SS)))$SSW -> SSW
  ( sub %>% group_by(name) %>% summarize( SS = N*sum( (mean(x) - global_mean)^2)) %>% summarize( SSB = sum(SS)))$SSB -> SSB
  # 2*sqrt(SSB/N) is the same thing
  return(2*sqrt(SST/N-SSW/N))
  
}

# Summary statistics that use population values

Bias = function( res, pop_val ){
  return( mean(res - pop_val))
}


StandardDeviation = function( res, pop_val){
  return( mean( (res - pop_val)*(res - pop_val)))
}

Median = function(p){
 return(p$year[min(which(cumsum(p$probability) >= 0.5))])
}



