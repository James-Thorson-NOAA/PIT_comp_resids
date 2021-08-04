

simulate_data <-
function(
         Z = 0.2,
         selex_50 = 3,
         selex_slope = 1,
         n_samp = 100,
         sigmaR = 0.6,
         max_age = 30,
         output = "C_a" ){

  # define dimensions
  age_set = 1:max_age

  # make abundance-at-age
  R_equil = 1e9
  N_a = R_equil * exp(-Z*age_set + rnorm(length(age_set), mean=-sigmaR^2/2, sd=sigmaR) )

  # make selectivity
  selex_a = plogis( selex_slope * (age_set-selex_50) )

  # simulate sampling
  prob_a = (N_a*selex_a) / sum(N_a*selex_a)
  C_a = rmultinom( n=1, size=n_samp, prob=prob_a )[,1]

  #
  if(output == "C_a") return(C_a)
  if(output == "prob_a") return(prob_a)
}

ddm <-
function(
         x,
         prob,
         size = sum(x),
         theta,
         beta = theta * size,
         log = FALSE,
         output = "loglike" ){

  loglike = lgamma(size+1) - sum(lgamma(x+1)) + lgamma(theta*size) - lgamma(size + beta)
  for( xI in seq_along(x) ){
    loglike = loglike + lgamma(x[xI] + beta*prob[xI]) - lgamma(beta * prob[xI]);
  }
  n_effective = 1/(1+theta) + size*(theta/(1+theta))
  if(output=="loglike") return(loglike);
  if(output=="n_effective") return(n_effective);
}

dmultinom <-
function(x, ..., theta, log, output="loglike"){
  if(output=="loglike") return(stats::dmultinom(x=x, ..., log=TRUE))
  if(output=="n_effective") return(sum(x))
}

#Z = 0.2
#selex_50 = 3
#selex_slope = 1
#n_samp = 100
#sigmaR = 0
#max_age = 20

fit_model <-
function(
         parlist,
         parvec = unlist(parlist),
         C_a,
         likelihood = dmultinom,
         output = "nll",
         sigmaR = 0,
         n_samp = 100 ){

  #
  params = relist( parvec, skeleton = parlist )

  #
  prob_a = simulate_data(
    selex_50 = params$selex_50,
    selex_slope = params$selex_slope,
    Z = params$Z,
    output = "prob_a",
    sigmaR = sigmaR    # Should be 0 during fitting
  )
  n_effective = likelihood(x=C_a, prob=prob_a, theta=exp(params$ln_theta), output="n_effective")

  # Replace this simulation for the DM with direct draws from a compound Dirichlet-multinomial process
  Csim_ar = matrix(NA, nrow=length(prob_a), ncol=100)
  for( rI in 1:100 ){
    n_sim = floor(n_effective) + ifelse(runif(1)<(n_effective%%1),1,0)
    Csim_ar[,rI] = simulate_data(
      selex_50 = params$selex_50,
      selex_slope = params$selex_slope,
      Z = params$Z,
      output = "C_a",
      sigmaR = sigmaR,  # Should be 0 during fitting
      n_samp = n_sim
    )
    Csim_ar[,rI] = Csim_ar[,rI] / sum(Csim_ar[,rI]) * n_samp
  }

  #
  nll = -1 * likelihood(x=C_a, prob=prob_a, theta=exp(params$ln_theta), log=TRUE)
  Chat_a = prob_a * sum(C_a)
  pearson_a = (C_a-Chat_a) / sqrt(n_effective * prob_a * (1-prob_a)) / (sum(C_a)/n_effective) # 3rd piece is a correction for sample sizes vs. DM variance
  PIT_a = runif( n=length(prob_a), min=rowSums(C_a%o%rep(1,100)<Csim_ar), max=rowSums(C_a%o%rep(1,100)<=Csim_ar) ) / 100

  #
  qPIT_a = qnorm( runif(n=length(PIT_a), min=PIT_a*99/100, max=1/100+PIT_a*99/100) )

  #
  if(output == "nll") return(nll)
  if(output == "residuals") return(cbind("pearson_a"=pearson_a,"C_a"=C_a,"Chat_a"=Chat_a,"prob_a"=prob_a,"PIT_a"=PIT_a,"qPIT_a"=qPIT_a))
  if(output == "Csim_ar") return(Csim_ar)
}

parlist = list(
  Z = 0.2,
  selex_50 = 3,
  selex_slope = 1,
  ln_theta = 0
)

if( FALSE ){
  # Explore
  C_a = simulate_data( sigmaR = 0.6 )
  parvec = unlist(parlist)

  # Run with initial value
  fit_model( parvec=unlist(parlist), parlist=parlist, C_a=C_a, likelihood=ddm )
  # Optimize
  parhat = nlminb(
    start = unlist(parlist),
    parlist=parlist,
    C_a=C_a,
    objective = fit_model,
    likelihood = ddm,
    control = list(trace=1)
  )
  # Get Pearson residuals
  fit_model( parvec=parhat$par, parlist=parlist, C_a=C_a, output="pearson_a" )
  fit_model( parvec=parhat$par, parlist=parlist, C_a=C_a, output="residuals" )
}

##############
# Simulation
# Problems with Pearson arise with n_samp <= 100, i.e., low bin observations
##############

sigmaR = 1.0  # If >0, then model is mis-specified with respective to disperion
n_samp = 100 # Affects power to detect mis-specification
likelihood = list( dmultinom, ddm )[[2]]

qPIT_ar = PIT_ar = pearson_ar = matrix(NA, nrow=30, ncol=100)

for( rI in 1:100 ){
  set.seed(rI + 101)
  # Simulate
  C_a = simulate_data(
    sigmaR = sigmaR,
    n_samp = n_samp
  )
  # Optimize
  parhat = nlminb(
    start = unlist(parlist),
    parlist=parlist,
    C_a=C_a,
    objective = fit_model,
    sigmaR = 0,
    n_samp = n_samp,
    likelihood = likelihood,
    control = list(trace=1)
  )
  # Get Pearson residuals
  resid_table = fit_model( parvec=parhat$par, parlist=parlist, C_a=C_a, n_samp=n_samp, likelihood=likelihood, output="residuals" )
  pearson_ar[,rI] = resid_table[,'pearson_a']
  PIT_ar[,rI] = resid_table[,'PIT_a']
  qPIT_ar[,rI] = resid_table[,'qPIT_a']
  # Csim_ar = fit_model( parvec=parhat$par, parlist=parlist, C_a=C_a, n_samp=n_samp, output="Csim_ar" )
}

Hist = function(x, threshold=5, ...){
  x = ifelse(x < -1*threshold, -1*threshold, x)
  x = ifelse(x > threshold, threshold, x)
  hist(x, ...)
}

setwd("C:/Users/James.Thorson/Desktop/Git/PIT_comp_resids")
ThorsonUtilities::save_fig( file=paste0(getwd(),"/demo"), width=6, height=6 )
  par(mfrow=c(2,1))
  xmax = 5.1
  xset = seq( -1*xmax, xmax, by=0.1)
  # Pearson
  Hist(pearson_ar, prob=TRUE, breaks=xset, main="Pearson" )
  lines( x=xset, y=dnorm(xset), col="blue", lwd=2 )
  # PIT
  #hist(PIT_ar, prob=TRUE, breaks=seq(0,1,by=0.01) )
  #abline(h=1, col="blue", lwd=2)
  # qPIT
  Hist(qPIT_ar, prob=TRUE, breaks=xset, main="PIT" )
  lines( x=xset, y=dnorm(xset), col="blue", lwd=2 )
dev.off()

