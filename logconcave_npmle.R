library("rmutil")
library(lsei)
library(scales)
library(jmuOutlier)
library(parallel)
library(doSNOW)
source('~/Desktop/log_concave_compare_老电脑修改版_M1000/cnmlcd.R')

par(mfrow=c(1,1))
#################### sample from the true function ##########################
#seedd = 100
seedd = 111
set.seed(seedd)
n<-500 # number of observations 50 200 500 2500 10000///1000 200
m<- ceiling(log(n)*n^(1/5)) # the number of points, it should be taken to be C*n^1/5
intpoints<- 150 # number of points used to compute the integral (normalising term) 150
S<- 10000  #number of MCMC samples 500
N= 100#bootstrapping iteration numbers

### gamma distribution
alpha=2
#x<-rgamma(n, alpha, 1)
#x<-rnorm(n, 0, 0.1)
x<-rnorm(n, 100, 0.1)
#############################################
r0 <- cnmlcd(x)

M=ceiling(S/2)#bootstrapping sample sizes


x_ext=x

list_r = list(r0$lcd)


boots_log_concave <- function(M,r0,x_ext,N){
  r=r0
  for (j in 1:N){
    new_x = rlcd(r$lcd)
    x_ext = c(x_ext,new_x)
    r <- cnmlcd(x_ext)
  }
  return(r$lcd)
}

#res <- lapply(1:M, boots_log_concave,r0=r0,x_ext=x_ext,N=N)

#cl <- makeSOCKcluster(3)
#x <- pbSapply(cl, 1:100, function(i, j) {Sys.sleep(1); i + j}, 100)

pbSapply <- function(cl, X, FUN, ...) {
  registerDoSNOW(cl)
  pb <- txtProgressBar(max=length(X))
  on.exit(close(pb))
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  foreach(i=X,  .options.snow=opts) %dopar% {
    FUN(i, ...)
  }
}

clnum<-detectCores() 

cl <- makeCluster(getOption("cl.cores", clnum-1));

clusterExport(cl, c("rlcd","lcdmax","dlcd","indx","cnmlcd","x.weight","new.lcd",
                    "logLik.lcd","maxima.gradient","gradient.lcd","cpx","pnnls","line.lcd","simplify.lcd"), 
              envir=environment())

start_time1 <- Sys.time()
#list_r <- parLapply(cl, 1:M, boots_log_concave,r0=r0,x_ext=x_ext,N=N)
list_r <- pbSapply(cl, 1:M, boots_log_concave,r0=r0,x_ext=x_ext,N=N)
end_time1 <- Sys.time()
time1=end_time1 - start_time1
print(time1)


stopCluster(cl);

###################################### Functions ################################

### the concave, piecewise-linear function in the exponent of the log-concave density:
wprop <- function(theta, p, gamma1, gamma2,x)
{
  sumterms<-0
  for (j in 1:m){
    sumterms<-sumterms+min(theta[j],(x-an))/theta[j]*p[j]
  }	
  return(sumterms*gamma1-(x-an)*gamma2)
}

### the likelihood function: 
like <- function(theta, p, gamma1, gamma2){
  wlikeprop<-0
  for (i in 1:n){
    wlikeprop<- wlikeprop+wprop(theta,p,gamma1,gamma2,x[i])
  }
  return(exp(wlikeprop-n*log( normconst(theta,p,gamma1,gamma2) ) ) )
}


### likelihood ratio: important to work on the log-level if one divides small function values:
likeratio <- function(thetaprop,pprop,gamma1prop,gamma2prop,theta, p, gamma1, gamma2){
  wlikeprop<-0
  wlikepropprop<-0 
  for (i in 1:n){
    wlikeprop<- wlikeprop+wprop(theta,p,gamma1,gamma2,x[i])
    wlikepropprop<- wlikepropprop+wprop(thetaprop,pprop,gamma1prop,gamma2prop,x[i])
  }
  return(exp(wlikepropprop-wlikeprop +n*log( normconst(theta,p,gamma1,gamma2) ) -n*log( normconst(thetaprop,pprop,gamma1prop,gamma2prop) ) ) )
}


### the normalising constant corresponding to the prior:
normconst <- function(theta, p, gamma1, gamma2)
{
  int<-0
  for (i in 1:intpoints){
    int<-int+exp(wprop(theta,p,gamma1,gamma2, points[i]))
  }
  return(int/intpoints*(bn-an))
}


########################################## EB method #####################

an<-min(x)
bn<-max(x)

# deterministics choice of the end points:
#an<- -log(n)  
#bn<- log(n)  



points<-c(1:intpoints)/intpoints*(bn-an)+an  # design points for the numerical intergration


################# MCMC method #############################
######## Initialization
#Initial value
#theta<- runif(m,an,bn)
counter<-0
theta<-c(1:m)/m*(bn-an)


#p<- (1:m)/(m*(m+1)/2)
p<-rep(1/m,m)
stickbreakbeta<-rep(0,(m-1))
stickbreakbeta[1]<-p[1]
stickbreakbeta[2]<-(p[2])/(1-stickbreakbeta[1])
for(j in 3:(m-1)){
  stickbreakbeta[j]<-p[j]/p[j-1]*stickbreakbeta[j-1]/(1-stickbreakbeta[j-1])
}



######################## Initialisation: 
gamma1<-10
gamma2<-1/2
fvalue<-matrix(0,S+1,intpoints) # collection of the sampled log-concave function values
modes<-rep(0,S)

# normalizing constant:
constant<-normconst(theta,p,gamma1,gamma2) 

# starting function in the MCMC
for (i in 1: intpoints){
  fvalue[1,i]<-exp( wprop(theta, p, gamma1, gamma2,points[i]))/constant
}
# Acceptance probabilities in the Gibbs cycle (for each parameter):
accept1<-0
accept2<-0
accept3<-0
accept4<-0


##################################  MCMC sampling: Metropolis Hastings within Gibbs sampling #######################################################
# we sparsify the sample and keep only the functional parameter value in the end of the Gibbs cycle
start_time2 <- Sys.time()
for (i in 1:S){
  if (i%%50==0){print(i/S)}
  ### theta:
  thetaprop<-theta+rnorm(m)/30  # random walk MH: small, symmetric jumps in each iteration for proposal
  for (j in 1: m){ 
    if (thetaprop[j]<0) {thetaprop[j]<-(0.1)}
    if (thetaprop[j]>bn-an) {thetaprop[j]<-(bn-an-0.1)}
  } # if the number is out of boundary, then bounce back
  #### accept rejact step:
  u<- runif(1,0,1) # the randomisation step
  likefrac<- likeratio(thetaprop, p, gamma1, gamma2,theta, p, gamma1, gamma2) # posterior proportion = likelihood ratio, due to the uniform prior
  if (likefrac > u ){ 
    theta<-thetaprop
    constant<-normconst(theta,p,gamma1,gamma2)
    for (k in 1: intpoints){
      fvalue[i+1,k]<-exp( wprop(theta, p, gamma1, gamma2,points[k]))/constant
    }
    accept1<-accept1+1 # number of accepted iterations, should be around 23%
  }else {fvalue[i+1,]<-fvalue[i,]}
  
  ### p: Similar as above, again choose Dirichlet prior with parameters 1,1,1,..,1 hence in the posterior ratio only the likelihood appears, prior cancels.
  # proposal
  pprop<-abs(p+rnorm(m)/100)
  pprop<-pprop/sum(pprop)
  
  stickbreakbeta_prop<-rep(0,(m-1))
  stickbreakbeta_prop[1]<-pprop[1]
  stickbreakbeta_prop[2]<-(pprop[2])/(1-stickbreakbeta_prop[1])
  for(j in 3:(m-1)){
    stickbreakbeta_prop[j]<-pprop[j]/pprop[j-1]*stickbreakbeta_prop[j-1]/(1-stickbreakbeta_prop[j-1])
  }
  log_prior_frac<-sum(log(dbeta(stickbreakbeta_prop,1,m/2)))-sum(log(dbeta(stickbreakbeta,1,m/2)))
  
  # accept/reject step:
  u<- runif(1,0,1)
  likefrac<-likeratio(theta, pprop, gamma1, gamma2,theta, p, gamma1, gamma2)
  if ((likefrac*exp(log_prior_frac))>u ){
    p<-pprop
    stickbreakbeta<-stickbreakbeta_prop
    constant<-normconst(theta,p,gamma1,gamma2)
    for (k in 1: intpoints){
      fvalue[i+1,k]<-exp( wprop(theta, p, gamma1, gamma2,points[k]))/constant
    }
    accept2<-accept2+1
  }
  
  ### gamma1: For gamma1 we choose Cauchy prior
  #proposal
  gamma1prop<-abs(gamma1+rnorm(1))
  likefrac<- likeratio(theta, p, gamma1prop, gamma2,theta, p, gamma1, gamma2)
  postfrac<-likefrac*dcauchy(gamma1prop)/dcauchy(gamma1) 
  
  #accept/reject step:
  u<- runif(1,0,1)
  if (postfrac>u ){
    gamma1<-gamma1prop
    constant<-normconst(theta,p,gamma1,gamma2)
    for (k in 1: intpoints){
      fvalue[i+1,k]<-exp( wprop(theta, p, gamma1, gamma2,points[k]))/constant
    }
    accept3<-accept3+1
  }
  
  ### gamma2: For gamma2 we also choose Cauchy prior
  #proposal:
  gamma2prop<-abs(gamma2+rnorm(1)/10)
  #likefrac<- like(theta, p, gamma1, gamma2prop)/like(theta, p, gamma1, gamma2)
  likefrac<- likeratio(theta, p, gamma1, gamma2prop,theta, p, gamma1, gamma2)
  postfrac<-likefrac*dcauchy(gamma2prop)/dcauchy(gamma2)  
  
  #accept/reject step
  u<- runif(1,0,1)
  if (postfrac>u ){
    gamma2<-gamma2prop
    constant<-normconst(theta,p,gamma1,gamma2)
    for (k in 1: intpoints){
      #fvalue[(i-1)*4+5,k]<-exp( wprop(theta, p, gamma1, gamma2,points[k]))/constant
      fvalue[i+1,k]<-exp( wprop(theta, p, gamma1, gamma2,points[k]))/constant
    }
    accept4<-accept4+1
  }
  
  
  
  ### compute the mode
  LCmode<-1
  for (k in 2: intpoints){
    if (fvalue[i+1,LCmode]< fvalue[i+1,k]){LCmode<-k}
  }
  modes[i]<-LCmode/intpoints*(bn-an)+an
}
end_time2 <- Sys.time()
time2=end_time2 - start_time2
print(time2)



### We throw out the first S/2 sample (burn in period). Then for each design point we order the posterior draws to get the quantiles
for (i in 1:intpoints){
  vecthelp<-sort(fvalue[(S/2):(S+1),i], decreasing = FALSE)
  fvalue[(S/2):(S+1),i]<-vecthelp
}


#### Plotting the posterior mean and the 95% pointwise confidence band
name=paste("Density bootstrap for Gamma(2,1), n=", n, sep="")
#hist(x,freq=FALSE, xlim=c(an,bn),breaks=20, main=name)


#set.seed(1)


#name=paste("Bayesian Bootstrap Densities for Gamma(2,1), n=", n, sep="")
#name=paste("Bayesian Bootstrap Densities for Normal(0,sd=0.1), n=", n, sep="")
name=paste("Bayesian Bootstrap Densities for Normal(100,sd=0.1), n=", n, sep="")
s=plot(r0$lcd,x,log = F,xlim=c(r0$lcd$lower,r0$lcd$upper),col = alpha("red"),lwd=3, pch=2, main=name)

#plot(list_r[[1]],log = F,xlim=c(list_r[[1]]$lower,list_r[[1]]$upper),col = alpha("black",0.05),ylim=s)
for (i in 1:M){
  par(new=T)
  plot(list_r[[i]],log = F,xlim=c(list_r[[1]]$lower,list_r[[1]]$upper),col = alpha("black",0.05),ylim=s)
}

#par(new=T)
#plot(r0$lcd,log = F,xlim=c(r0$lcd$lower,r0$lcd$upper),col = alpha("red"),lwd=3, pch=2)
#temp=(1:100000)/10000
#gm21 = function(x){return(dgamma(x,2,1))}
#lines(temp,gm21(temp),col=alpha("purple"))


lines(points,fvalue[(S/2+S/2*0.975),],col=28, lty=2)
lines(points,fvalue[(S/2+S*0.025),],col=28,lty=2)
lines(points,colMeans(fvalue[(S/2):(S+1),]),col=28)
#lines(points, mle, col=3, lwd=2)
#lines(points,dnorm(points),col=2)
#lines(points,dgumbel(points,1,1),col=43)
#lines(points,dbeta(points,1,2),col=2)
#lines(points,dgamma(points,alpha,1),col=43,lwd=2)
lines(points,dnorm(points,100,0.1),col=43,lwd=2)
#lines(points,dnorm(points,0,0.1),col=43,lwd=2)


#hist(modes[(S/2+1):S],freq=FALSE,main = paste("n=", n),xlab="mode",ylab="relative frequency")
time2>=time1
time1
time2


#save.image(file='~/Desktop/environment/p???.RData')

