# This is a R code for Bayesian Statictics course project: 
# ***District level included***
# ***A Bayesian Heirarchial Model for Educational Infrastructure Qualityusing School Effectiveness, and other Spatial/Regional Information***
# Authors: Zheng Li, Chamari, Sulalitha, Miach 
# Southern Methodist University
# Dec, 2018
#--------------------------- conditional posterior distribution function ------------------ 
# conditonal posterior distribution for sigma
sigma.cond.post<- function(X, y, g, mu, beta, alpha, lambda, sigma)
{
  u <- beta-mu;
  sig<- solve(g*sigma^{2}*solve(t(X)%*%X));
  part1 <- 1/(det(sig))^{0.5};
  part2 <- exp(-0.5*t(u)%*%sig%*%u)*exp(-lambda*sigma)*sigma^{alpha-1}
  part1*part2
}

# conditonal posterior distribution for beta
beta.cond.post <- function(g, mu, beta, X, y, alpha, lambda, sigma){
  # returns a scalar
  # Note: beta is a vectors beta = [beta0, beta1, beta2,....]
  s1<- rep(0,dim(X)[1]);
  u <- beta-mu;
  XX <- t(X)%*%X
  #er <- matrix(0,dim(XX)[1],dim(XX)[2])
  #diag(er) <- 0.0001
  #XX <-XX + er
  sig<-g*sigma^{2}*solve(XX);
  part2<-exp(-0.5*t(u)%*%solve(sig)%*%u);

  XY <-cbind(X,y);
  s1<- apply(XY,1,bino_beta, beta.v = beta);
  part1<- sum(s1)
  part1*part2
}
# Compute binomial beta function
bino_beta<- function(xy,beta.v){
  # return a scalar
  y<-xy[length(xy)]
  x<-as.vector(xy[-1])
  p1<-(1/(1+exp(-t(beta.v)%*%x)))^{y};
  p2<-(1+exp(t(beta.v)%*%x))^{y-1};
  p1*p2
}

#------------------------------------ MH sampling function ------------------------------------
draw.sigma <- function(X, y, g, mu, beta, alpha, lambda, sigma) {
  T=1000;
  B=500;
  accept <- 0 ;
  sigma.store <- rep(0, T); 
  for (i in 1:T){
    sigma.star <- rgamma(1,alpha,lambda) # propose a random value from Gamma(1,1)
    logMH.part1 <- log(sigma.cond.post(X,
                                       y,
                                       g,
                                       mu,
                                       beta,
                                       alpha,
                                       lambda, 
                                       sigma.star)) - log(sigma.cond.post(X,
                                                                          y,
                                                                          g,
                                                                          mu,
                                                                          beta,
                                                                          alpha,
                                                                          lambda, 
                                                                          sigma))
    
    logMH.part2 <- dgamma(sigma,alpha,lambda,log=TRUE) - dgamma(sigma.star,alpha,lambda,log=TRUE)
    logMH <- logMH.part1 + logMH.part2
    # MH ratio
    MH = exp(logMH);
    u <- runif(1); 
    #wait for buring period to finish
    if (i<B){
      sigma.store[i] <- -1;
    }else{
      # After burn-in
      if (u < MH) {
        sigma <- sigma.star;
        sigma.store[i] <- sigma;
        accept <- accept +1 ;
      }else{
        sigma.store[i] <- -1;
      } 
    }
  }
  sigma.est = sigma.store[sigma.store>0]
  #t<-sort(table(alpha.est),decreasing=TRUE)
  #return(as.double(names(t)[1]))
  
  # only return the median of sampled sigma values
  cat(paste0(round(median(sigma.est),3)," (AR:",round(accept/(T-B), 3),")"));
  return(median(sigma.est))
}

draw.beta <- function(g, mu, beta, X, y, alpha, lambda, sigma) {
  # draw a vector of beta
  T=1000;
  B = 500;
  accept <-0;
  beta.store <- matrix(0, T, dim(X)[2]); # beta.store: T x dim(X)[2] matrix
  for (i in 1:T){
    beta.star <- as.vector(rzellner(1, g=g, sigma=sigma, X=X))
    logMH.part1 <-  log(beta.cond.post(g, 
                                       mu, 
                                       beta.star, 
                                       X, 
                                       y, 
                                       alpha, 
                                       lambda, 
                                       sigma)) - log(beta.cond.post(g, 
                                                                    mu, 
                                                                    beta, 
                                                                    X, 
                                                                    y,
                                                                    alpha, 
                                                                    lambda, 
                                                                    sigma))
    logMH.part2 <- dzellner(beta, g, sigma, X, log=TRUE) - dzellner(beta.star, g, sigma, X, log=TRUE)
    logMH <- logMH.part1 + logMH.part2
    # MH ratio
    MH = exp(logMH);
    u <- runif(1); 
    #burn-in period
    if (i<B){
      beta.store[i,] <- -1;
    }else{
      # After burn-in
      if (u < MH) {
        beta <- beta.star
        beta.store[i,] <- beta; # store vector
        accept <- accept+1;
      }else{
        beta.store[i,] <- -1;
      } 
    }
  }
  rowsum <- apply(beta.store, 1, sum)
  beta.est = beta.store[rowsum>0,]
  # return a 1xdim(X)[2] vector for all the betas
  cat(paste0("(AR:",round(accept/(T-B), 3),")"));
  #print(is.null(dim(beta.est)))
  if (is.null(dim(beta.est))){
      return(beta.est)
  }else{
      return(apply(beta.est, 2, median))
  }
}

#------------------------------------ MH sampler function ------------------------------------
GibbsMH <-function(seed=0xFACE,
                   iter=iter, 
                   mhiter = mhT, 
                   beta.init,
                   sigma.init,
                   X,
                   g,
                   mu,
                   alpha,
                   lambda){
  set.seed(seed);
  # The sampler will iterate T times and store estimated paramtesr into .tot variables
  
  # sigma.tot <- rep(0, iter); 
  beta.tot <- matrix(0,iter,8)
  
  for (t in 1:iter){
    cat(paste0("\nSimulation ", t,":"))
    # Now perform the metropolis steps for simulating beta and sigma using
    # a g prior distribution (multivariate normal)
    
    # Comment out using a fixed sigma instead
    # draw sigma
    # cat(' sigma-> ')
    # sigmac <- draw.sigma(X,
    #                      y,
    #                      g,
    #                      mu,
    #                      beta.init,
    #                      alpha,
    #                      lambda,
    #                      sigma.init)
    
    # draw beta
    cat(' beta->')
    betac <- draw.beta(g, 
                       mu, 
                       beta.init, 
                       X, 
                       y, 
                       alpha, 
                       lambda, 
                       sigma.init)
    
    
    
    # update initial values
    # if (length(sigmac)>0){
    #   sigma.init <- sigmac
    # }
    if (length(betac)>0){
      beta.init <- betac
    }
    # store the parameters being used
    # sigma.tot[t]<- sigma.init
    beta.tot[t,]<- beta.init
  }
  output<-list(beta.tot)
  return(output)
}

# -----------------------------------running code below-------------------------------------------#

##  MCMC: Gibbs Sampler with Metropolis Hastings steps
initial.estimates<- function(dataframe)
{
  # first fit a full logistic model to obtain initial estiamtes of coeffcients
  model <- glm(SEI~School.Type+Diversity+Choice.Magnet.School+X2018..Distinctions+AverageDistance+NumberOfRail,family=binomial(link='logit'),data=dataframe)
  X<- model.matrix(model)
  y<- dataframe$SEI
  coeff <- as.vector(model$coefficients)
  output<-list(X,y,coeff)
  return(output)
}

setwd("C:/Users/47513405/OneDrive - Southern Methodist University/SMU/2018Fall/Bayesian Statistics/Project")
## first fit a full logistic model to obtain initial estiamtes of coeffcients
df <- read.csv("data_final.csv")


#--------------------------- Initialize for the values of the subchains -----------------------------
iter <-2000; # num. of iterations of the chain
mhT <- 500; # num. of burn-in iterations
alpha <-9; lambda<-4;
sigma.init <-3.316843;

#------------------------------------- District 1 level parameters ----------------------------------------
df.district1 <- df[df$Trustee.District==1|df$Trustee.District==3|df$Trustee.District==2|df$Trustee.District==8,]
est.init <- initial.estimates(df.district1)
est.init[[3]][5]<-0
X <- est.init[[1]]
y <- est.init[[2]]
beta.init <-est.init[[3]]
D<- dim(X); 
g <- D[1];
mu <- rep(0,D[2]);
output.district1 <- GibbsMH(seed=0xFACE,
                            iter=iter, 
                            mhiter = mhT, 
                            beta.init= beta.init,
                            sigma.init = sigma.init,
                            X = X,
                            g = g,
                            mu = mu,
                            alpha = alpha,
                            lambda = lambda)

#------------------------------------- District 2 level parameters ----------------------------------------
df.district2 <- df[df$Trustee.District==4|df$Trustee.District==9,]
est.init <- initial.estimates(df.district2)
X <- est.init[[1]]
y <- est.init[[2]]
beta.init <-est.init[[3]]
D<- dim(X); 
g <- D[1];
mu <- rep(0,D[2]);
output.district2 <- GibbsMH(seed=0xFACE,
                            iter=iter, 
                            mhiter = mhT, 
                            beta.init= beta.init,
                            sigma.init = sigma.init,
                            X = X,
                            g = g,
                            mu = mu,
                            alpha = alpha,
                            lambda = lambda)

#------------------------------------- District 3 level parameters ----------------------------------------
df.district3 <- df[df$Trustee.District==5|df$Trustee.District==6|df$Trustee.District==7|df$Trustee.District==2,]
est.init <- initial.estimates(df.district3)
X <- est.init[[1]]
y <- est.init[[2]]
beta.init <-est.init[[3]]
D<- dim(X); 
g <- D[1];
mu <- rep(0,D[2]);
output.district3 <- GibbsMH(seed=0xFACE,
                            iter=iter, 
                            mhiter = mhT, 
                            beta.init= beta.init,
                            sigma.init = sigma.init,
                            X = X,
                            g = g,
                            mu = mu,
                            alpha = alpha,
                            lambda = lambda)

#------------------------------------- District 4 level parameters ----------------------------------------
df.district4 <- df[df$Trustee.District==4,]
est.init <- initial.estimates(df.district4)
X <- est.init[[1]]
y <- est.init[[2]]s
beta.init <-est.init[[3]]
D<- dim(X); 
g <- D[1];
mu <- rep(0,D[2]);
output.district4 <- GibbsMH(seed=0xFACE,
                            iter=iter, 
                            mhiter = mhT, 
                            beta.init= beta.init,
                            sigma.init = sigma.init,
                            X = X,
                            g = g,
                            mu = mu,
                            alpha = alpha,
                            lambda = lambda)

#------------------------------------- District 5 level parameters ----------------------------------------
df.district5 <- df[df$Trustee.District==5,]
est.init <- initial.estimates(df.district5)
X <- est.init[[1]]
y <- est.init[[2]]
beta.init <-est.init[[3]]
D<- dim(X); 
g <- D[1];
mu <- rep(0,D[2]);
output.district5 <- GibbsMH(seed=0xFACE,
                            iter=iter, 
                            mhiter = mhT, 
                            beta.init= beta.init,
                            sigma.init = sigma.init,
                            X = X,
                            g = g,
                            mu = mu,
                            alpha = alpha,
                            lambda = lambda)

#------------------------------------- District 6 level parameters ----------------------------------------
df.district6 <- df[df$Trustee.District==6,]
est.init <- initial.estimates(df.district6)
X <- est.init[[1]]
y <- est.init[[2]]
beta.init <-est.init[[3]]
D<- dim(X); 
g <- D[1];
mu <- rep(0,D[2]);
output.district6 <- GibbsMH(seed=0xFACE,
                            iter=iter, 
                            mhiter = mhT, 
                            beta.init= beta.init,
                            sigma.init = sigma.init,
                            X = X,
                            g = g,
                            mu = mu,
                            alpha = alpha,
                            lambda = lambda)

#------------------------------------- District 7 level parameters ----------------------------------------
df.district7 <- df[df$Trustee.District==7,]
est.init <- initial.estimates(df.district7)
X <- est.init[[1]]
y <- est.init[[2]]
beta.init <-est.init[[3]]
D<- dim(X); 
g <- D[1];
mu <- rep(0,D[2]);
output.district7 <- GibbsMH(seed=0xFACE,
                            iter=iter, 
                            mhiter = mhT, 
                            beta.init= beta.init,
                            sigma.init = sigma.init,
                            X = X,
                            g = g,
                            mu = mu,
                            alpha = alpha,
                            lambda = lambda)

#------------------------------------- District 8 level parameters ----------------------------------------
df.district8 <- df[df$Trustee.District==8,]
est.init <- initial.estimates(df.district8)
X <- est.init[[1]]
y <- est.init[[2]]
beta.init <-est.init[[3]]
D<- dim(X); 
g <- D[1];
mu <- rep(0,D[2]);
output.district8 <- GibbsMH(seed=0xFACE,
                            iter=iter, 
                            mhiter = mhT, 
                            beta.init= beta.init,
                            sigma.init = sigma.init,
                            X = X,
                            g = g,
                            mu = mu,
                            alpha = alpha,
                            lambda = lambda)

#------------------------------------- District 9 level parameters ----------------------------------------
df.district9 <- df[df$Trustee.District==9,]
est.init <- initial.estimates(df.district9)
X <- est.init[[1]]
y <- est.init[[2]]
beta.init <-est.init[[3]]
D<- dim(X); 
g <- D[1];
mu <- rep(0,D[2]);
output.district9 <- GibbsMH(seed=0xFACE,
                            iter=iter, 
                            mhiter = mhT, 
                            beta.init= beta.init,
                            sigma.init = sigma.init,
                            X = X,
                            g = g,
                            mu = mu,
                            alpha = alpha,
                            lambda = lambda)
# ---------------------------------------- Results plots -----------------------------------------

# get the data separately 
run1 <-output.district1
run2 <-output.district2
run3 <-output.district3

box.data1 <- data.frame('Sim1',
                        run1[[1]][,1],
                        run1[[1]][,2],
                        run1[[1]][,3],
                        run1[[1]][,4],
                        run1[[1]][,5],
                        run1[[1]][,6],
                        run1[[1]][,7],
                        run1[[1]][,8]); 
colnames(box.data1) <- c("Sim",
                         "beta0",
                         "beta1",
                         "beta2",
                         "beta3",
                         "beta4",
                         "beta5",
                         "beta6",
                         "beta7");

box.data2 <- data.frame('Sim2',
                        run2[[1]][,1],
                        run2[[1]][,2],
                        run2[[1]][,3],
                        run2[[1]][,4],
                        run2[[1]][,5],
                        run2[[1]][,6],
                        run2[[1]][,7],
                        run2[[1]][,8]); 
colnames(box.data2) <- c("Sim",
                         "beta0",
                         "beta1",
                         "beta2",
                         "beta3",
                         "beta4",
                         "beta5",
                         "beta6",
                         "beta7");
box.data3 <- data.frame('Sim3',
                        run3[[1]][,1],
                        run3[[1]][,2],
                        run3[[1]][,3],
                        run3[[1]][,4],
                        run3[[1]][,5],
                        run3[[1]][,6],
                        run3[[1]][,7],
                        run3[[1]][,8]); 
colnames(box.data3) <- c("Sim",
                         "beta0",
                         "beta1",
                         "beta2",
                         "beta3",
                         "beta4",
                         "beta5",
                         "beta6",
                         "beta7");
box.data <-rbind(box.data1,box.data2,box.data3)
pdf("boxplot_districts.pdf",width=8, height=5,family="Times")
layout(matrix(c(1:8), nrow=2, byrow=T))
par(mar=c(3,2,3,1))
boxplot(beta0~Sim,data = box.data,notch = TRUE,outline=FALSE,names = c('DG1','DG2','DG3'))
title(expression(beta[0]),line=0.7)

boxplot(beta1~Sim,data = box.data,notch = TRUE,outline=FALSE,names = c('DG1','DG2','DG3'))
title(expression(beta[1]),line=0.7)

boxplot(beta2~Sim,data = box.data,notch = TRUE,outline=FALSE,names = c('DG1','DG2','DG3'))
title(expression(beta[2]),line=0.7)

boxplot(beta3~Sim,data = box.data,notch = TRUE,outline=FALSE,names = c('DG1','DG2','DG3'))
title(expression(beta[3]),line=0.7)

boxplot(beta4~Sim,data = box.data,notch = TRUE,outline=FALSE,names = c('DG1','DG2','DG3'))
title(expression(beta[4]),line=0.7)

boxplot(beta5~Sim,data = box.data,notch = TRUE,outline=FALSE,names = c('DG1','DG2','DG3'))
title(expression(beta[5]),line=0.7)

boxplot(beta6~Sim,data = box.data,notch = TRUE,outline=FALSE,names = c('DG1','DG2','DG3'))
title(expression(beta[6]),line=0.7)
boxplot(beta7~Sim,data = box.data,notch = TRUE,outline=FALSE,names = c('DG1','DG2','DG3'))
title(expression(beta[7]),line=0.7)


dev.off()

## compute CPI for beta
cpi.data <-data.frame(t(as.vector(quantile(output.district1[[1]][,1],probs = c(2.5, 50, 97.5)/100))))
colnames(cpi.data)<- c('5%DG1','50%DG1','95%DG1')
for (i in 2:8){
  row<- data.frame(t(as.vector(quantile(output.district1[[1]][,i],probs = c(2.5, 50, 97.5)/100))))
  colnames(row)<- c('5%DG1','50%DG1','95%DG1')
  cpi.data<-rbind(cpi.data,row)
  }
cpi.DG1 <- cpi.data

cpi.data <-data.frame(t(as.vector(quantile(output.district2[[1]][,1],probs = c(2.5, 50, 97.5)/100))))
colnames(cpi.data)<- c('5%DG2','50%DG2','95%DG2')
for (i in 2:8){
  row<- data.frame(t(as.vector(quantile(output.district2[[1]][,i],probs = c(2.5, 50, 97.5)/100))))
  colnames(row)<- c('5%DG2','50%DG2','95%DG2')
  cpi.data<-rbind(cpi.data,row)
}
cpi.DG2 <- cpi.data

cpi.data <-data.frame(t(as.vector(quantile(output.district3[[1]][,1],probs = c(2.5, 50, 97.5)/100))))
colnames(cpi.data)<- c('5%DG3','50%DG3','95%DG3')
for (i in 2:8){
  row<- data.frame(t(as.vector(quantile(output.district3[[1]][,i],probs = c(2.5, 50, 97.5)/100))))
  colnames(row)<- c('5%DG3','50%DG3','95%DG3')
  cpi.data<-rbind(cpi.data,row)
}
cpi.DG3 <- cpi.data

cpi.data<-cbind(cpi.DG1,cpi.DG2,cpi.DG3)
# write to csv
write.csv(cpi.data,file='beta_CPI.csv')
