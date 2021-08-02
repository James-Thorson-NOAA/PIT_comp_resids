

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
    sigmaR = sigmaR
  )

  Csim_ar = matrix(NA, nrow=length(prob_a), ncol=100)
  for( rI in 1:100 ){
    Csim_ar[,rI] = simulate_data(
      selex_50 = params$selex_50,
      selex_slope = params$selex_slope,
      Z = params$Z,
      output = "C_a",
      sigmaR = sigmaR,
      n_samp = n_samp
    )
  }

  #
  nll = -1 * dmultinom(x=C_a, prob=prob_a, log=TRUE)
  Chat_a = prob_a * sum(C_a)
  pearson_a = (C_a-Chat_a) / sqrt(sum(C_a) * prob_a * (1-prob_a))
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
  selex_slope = 1
)

if( FALSE ){
  # Explore
  C_a = simulate_data( sigmaR = 0 )
  parvec = unlist(parlist)

  # Run with initial value
  fit_model( parvec=unlist(parlist), parlist=parlist, C_a=C_a )
  # Optimize
  parhat = nlminb(
    start = unlist(parlist),
    parlist=parlist,
    C_a=C_a,
    objective = fit_model
  )
  # Get Pearson residuals
  fit_model( parvec=parhat$par, parlist=parlist, C_a=C_a, output="pearson_a" )
  fit_model( parvec=parhat$par, parlist=parlist, C_a=C_a, output="residuals" )
}

##############
# Simulation
# Problems with Pearson arise with n_samp <= 100, i.e., low bin observations
##############

sigmaR = 0  # If >0, then model is mis-specified with respective to disperion
n_samp = 100 # Affects power to detect mis-specification

qPIT_ar = PIT_ar = pearson_ar = matrix(NA, nrow=30, ncol=100)

for( rI in 1:100 ){
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
    n_samp = n_samp
  )
  # Get Pearson residuals
  pearson_ar[,rI] = fit_model( parvec=parhat$par, parlist=parlist, C_a=C_a, n_samp=n_samp, output="residuals" )[,'pearson_a']
  PIT_ar[,rI] = fit_model( parvec=parhat$par, parlist=parlist, C_a=C_a, n_samp=n_samp, output="residuals" )[,'PIT_a']
  qPIT_ar[,rI] = fit_model( parvec=parhat$par, parlist=parlist, C_a=C_a, n_samp=n_samp, output="residuals" )[,'qPIT_a']
  # Csim_ar = fit_model( parvec=parhat$par, parlist=parlist, C_a=C_a, n_samp=n_samp, output="Csim_ar" )
}

Hist = function(x, threshold=5, ...){
  x = ifelse(x < -1*threshold, -1*threshold, x)
  x = ifelse(x > threshold, threshold, x)
  hist(x, ...)
}

ThorsonUtilities::save_fig( file=paste0("C:/Users/James.Thorson/Desktop/Work files/Collaborations/2021 -- Pearson residuals for comps/demo"), width=6, height=6 )
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

