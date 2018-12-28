# This is a R code for Bayesian Statictics course project: 
#***District level NOT included**
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
  sig<-g*sigma^{2}*solve(t(X)%*%X);
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
  T=1500;
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
  return(apply(beta.est, 2, median))
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
  
  sigma.tot <- rep(0, iter); 
  beta.tot <- matrix(0,iter,8)
  
  for (t in 1:iter){
    cat(paste0("\nSimulation ", t,":"))
    # Now perform the metropolis steps for simulating beta and sigma using
    # a g prior distribution (multivariate normal)
    
    # draw sigma
    cat(' sigma-> ')
    sigmac <- draw.sigma(X,
                         y,
                         g,
                         mu,
                         beta.init,
                         alpha,
                         lambda,
                         sigma.init)
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
    if (length(sigmac)>0){
      sigma.init <- sigmac
    }
    if (length(betac)>0){
      beta.init <- betac
    }
    # store the parameters being used
    sigma.tot[t]<- sigma.init
    beta.tot[t,]<- beta.init
  }
  output<-list(sigma.tot,beta.tot)
  return(output)
}

# -----------------------------------running code below-------------------------------------------#

##  MCMC: Gibbs Sampler with Metropolis Hastings steps

## first fit a full logistic model to obtain initial estiamtes of coeffcients
df <- read.csv("data_final.csv")
attach(df)
mod1 <- glm(SEI~School.Type+Diversity+Choice.Magnet.School+X2018..Distinctions+AverageDistance+NumberOfRail,family=binomial(link='logit'),data=df)
# summary(mod1)
X<- model.matrix(mod1)
y<- df$SEI
coeff <- as.vector(mod1$coefficients)

#--------------------------- Initialize for the values of the subchains -----------------------------
iter <- 2000; # num. of iterations of the chain
mhT <- 500; # num. of burn-in iterations

D<- dim(X); g <- D[1];
alpha <-5; lambda<-4;
mu <- rep(0,D[2]); 
beta.init<-coeff; sigma.init <-2;

#------------------------------------- Run MCMC algorithm --------------------------------------------
output <- GibbsMH(seed=0xFACE,
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

# plot(1:1000,output[[1]],type='l',xlab = 'Iterations',ylab = expression(sigma),main=expression(paste('Trace plot of ',sigma)))
# hist(output[[1]], main = expression(paste('histogram of ',sigma)), xlab = expression(sigma));
# quantile(output[[1]],probs = c(2.5, 50, 97.5)/100)
# 
# plot(1:1000,output[[2]][,1],type='l',xlab = 'Iterations',ylab = expression(beta[0]),main=expression(paste('Trace plot of ',beta[0])))
# hist(output[[2]][,1], main = expression(paste('histogram of ',beta[0])), xlab = expression(beta[0]));
# quantile(output[[2]][,1],probs = c(2.5, 50, 97.5)/100)
# 
# plot(1:1000,output[[2]][,2],type='l',xlab = 'Iterations',ylab = expression(beta[1]),main=expression(paste('Trace plot of ',beta[1])))
# hist(output[[2]][,2], main = expression(paste('histogram of ',beta[1])), xlab = expression(beta[1]));
# quantile(output[[2]][,2],probs = c(2.5, 50, 97.5)/100)
# 
# plot(1:1000,output[[2]][,3],type='l',xlab = 'Iterations',ylab = expression(beta[2]),main=expression(paste('Trace plot of ',beta[2])))
# hist(output[[2]][,3], main = expression(paste('histogram of ',beta[2])), xlab = expression(beta[2]));
# quantile(output[[2]][,3],probs = c(2.5, 50, 97.5)/100)
# 
# plot(1:1000,output[[2]][,4],type='l',xlab = 'Iterations',ylab = expression(beta[3]),main=expression(paste('Trace plot of ',beta[3])))
# hist(output[[2]][,4], main = expression(paste('histogram of ',beta[3])), xlab = expression(beta[3]));
# quantile(output[[2]][,4],probs = c(2.5, 50, 97.5)/100)
# 
# plot(1:1000,output[[2]][,5],type='l',xlab = 'Iterations',ylab = expression(beta[4]),main=expression(paste('Trace plot of ',beta[4])))
# hist(output[[2]][,5], main = expression(paste('histogram of ',beta[4])), xlab = expression(beta[4]));
# quantile(output[[2]][,5],probs = c(2.5, 50, 97.5)/100)
# 
# plot(1:1000,output[[2]][,6],type='l',xlab = 'Iterations',ylab = expression(beta[5]),main=expression(paste('Trace plot of ',beta[5])))
# hist(output[[2]][,6], main = expression(paste('histogram of ',beta[5])), xlab = expression(beta[5]));
# quantile(output[[2]][,6],probs = c(2.5, 50, 97.5)/100)
# 
# plot(1:1000,output[[2]][,7],type='l',xlab = 'Iterations',ylab = expression(beta[6]),main=expression(paste('Trace plot of ',beta[6])))
# hist(output[[2]][,7], main = expression(paste('histogram of ',beta[6])), xlab = expression(beta[6]));
# quantile(output[[2]][,7],probs = c(2.5, 50, 97.5)/100)
# 
# plot(1:1000,output[[2]][,8],type='l',xlab = 'Iterations',ylab = expression(beta[7]),main=expression(paste('Trace plot of ',beta[7])))
# hist(output[[2]][,8], main = expression(paste('histogram of ',beta[7])), xlab = expression(beta[7]));
# quantile(output[[2]][,8],probs = c(2.5, 50, 97.5)/100)
# 
# acf(output[[1]],main=expression(sigma))
# acf(output[[2]][,1],main=expression(beta[0]))
# acf(output[[2]][,2],main=expression(beta[1]))
# acf(output[[2]][,3],main=expression(beta[2]))
# acf(output[[2]][,4],main=expression(beta[3]))
# acf(output[[2]][,5],main=expression(beta[4]))
# acf(output[[2]][,6],main=expression(beta[5]))
# acf(output[[2]][,7],main=expression(beta[6]))
# acf(output[[2]][,8],main=expression(beta[7]))

setwd("/Users/ZhengLi/SMU_OneDrive/OneDrive - Southern Methodist University/SMU/2018Fall/Bayesian Statistics/Project")
pdf("beta_boxplot_9_4.pdf",width=8, height=5,family="Times")
layout(matrix(c(1:9), 3,3, byrow = T))
par(mar=c(4,5,3,3))
hist(output[[1]], main = expression(paste('histogram of ',sigma)), xlab = expression(sigma));

hist(output[[2]][,1], main = expression(paste('histogram of ',beta[0])), xlab = expression(beta[0]));

hist(output[[2]][,2], main = expression(paste('histogram of ',beta[1])), xlab = expression(beta[1]));

hist(output[[2]][,3], main = expression(paste('histogram of ',beta[2])), xlab = expression(beta[2]));

hist(output[[2]][,4], main = expression(paste('histogram of ',beta[3])), xlab = expression(beta[3]));


hist(output[[2]][,5], main = expression(paste('histogram of ',beta[4])), xlab = expression(beta[4]));


hist(output[[2]][,6], main = expression(paste('histogram of ',beta[5])), xlab = expression(beta[5]));


hist(output[[2]][,7], main = expression(paste('histogram of ',beta[6])), xlab = expression(beta[6]));


hist(output[[2]][,8], main = expression(paste('histogram of ',beta[7])), xlab = expression(beta[7]));

dev.off()

pdf("ACF_9_4.pdf",width=8, height=5,family="Times")
layout(matrix(c(1:9), 3,3, byrow = T))
par(mar=c(4,5,3,3))
acf(output[[1]],main=expression(sigma))
acf(output[[2]][,1],main=expression(beta[0]))
acf(output[[2]][,2],main=expression(beta[1]))
acf(output[[2]][,3],main=expression(beta[2]))
acf(output[[2]][,4],main=expression(beta[3]))
acf(output[[2]][,5],main=expression(beta[4]))
acf(output[[2]][,6],main=expression(beta[5]))
acf(output[[2]][,7],main=expression(beta[6]))
acf(output[[2]][,8],main=expression(beta[7]))
dev.off()

pdf("Trace_9_4.pdf",width=8, height=5,family="Times")
layout(matrix(c(1:9), 3,3, byrow = T))
par(mar=c(4,5,3,3))
iterations= 2000
plot(1:iterations,output[[1]],type='l',xlab = 'Iterations',ylab = expression(sigma),main=expression(paste('Trace plot of ',sigma)))

plot(1:iterations,output[[2]][,1],type='l',xlab = 'Iterations',ylab = expression(beta[0]),main=expression(paste('Trace plot of ',beta[0])))

plot(1:iterations,output[[2]][,2],type='l',xlab = 'Iterations',ylab = expression(beta[1]),main=expression(paste('Trace plot of ',beta[1])))

plot(1:iterations,output[[2]][,3],type='l',xlab = 'Iterations',ylab = expression(beta[2]),main=expression(paste('Trace plot of ',beta[2])))

plot(1:iterations,output[[2]][,4],type='l',xlab = 'Iterations',ylab = expression(beta[3]),main=expression(paste('Trace plot of ',beta[3])))

plot(1:iterations,output[[2]][,5],type='l',xlab = 'Iterations',ylab = expression(beta[4]),main=expression(paste('Trace plot of ',beta[4])))

plot(1:iterations,output[[2]][,6],type='l',xlab = 'Iterations',ylab = expression(beta[5]),main=expression(paste('Trace plot of ',beta[5])))

plot(1:iterations,output[[2]][,7],type='l',xlab = 'Iterations',ylab = expression(beta[6]),main=expression(paste('Trace plot of ',beta[6])))

plot(1:iterations,output[[2]][,8],type='l',xlab = 'Iterations',ylab = expression(beta[7]),main=expression(paste('Trace plot of ',beta[7])))

dev.off()

# get the data separately 
run1 <-output
run2 <-output
run3 <-output

box.data1 <- data.frame('Sim1',
                        run1[[1]],
                        run1[[2]][,1],
                        run1[[2]][,2],
                        run1[[2]][,3],
                        run1[[2]][,4],
                        run1[[2]][,5],
                        run1[[2]][,6],
                        run1[[2]][,7],
                        run1[[2]][,8]); 
colnames(box.data1) <- c("Sim",
                         'sigma',
                         "beta0",
                         "beta1",
                         "beta2",
                         "beta3",
                         "beta4",
                         "beta5",
                         "beta6",
                         "beta7");

box.data2 <- data.frame('Sim2',
                        run2[[1]],
                        run2[[2]][,1],
                        run2[[2]][,2],
                        run2[[2]][,3],
                        run2[[2]][,4],
                        run2[[2]][,5],
                        run2[[2]][,6],
                        run2[[2]][,7],
                        run2[[2]][,8]); 
colnames(box.data2) <- c("Sim",
                         'sigma',
                         "beta0",
                         "beta1",
                         "beta2",
                         "beta3",
                         "beta4",
                         "beta5",
                         "beta6",
                         "beta7");
box.data3 <- data.frame('Sim3',
                        run3[[1]],
                        run3[[2]][,1],
                        run3[[2]][,2],
                        run3[[2]][,3],
                        run3[[2]][,4],
                        run3[[2]][,5],
                        run3[[2]][,6],
                        run3[[2]][,7],
                        run3[[2]][,8]); 
colnames(box.data3) <- c("Sim",
                         'sigma',
                         "beta0",
                         "beta1",
                         "beta2",
                         "beta3",
                         "beta4",
                         "beta5",
                         "beta6",
                         "beta7");
box.data <-rbind(box.data1,box.data2,box.data3)
pdf("boxplot_sensitivity.pdf",width=8, height=5,family="Times")
layout(matrix(c(1:9), nrow=3, byrow=T))
par(mar=c(3,4,3,3))
boxplot(sigma~Sim,data = box.data,notch = TRUE,outline=FALSE,names = c(expression(paste(alpha,'=8, ',beta,'=2')),
                                                                       expression(paste(alpha,'=9, ',beta,'=3')),
                                                                       expression(paste(alpha,'=9, ',beta,'=4'))))
title(expression(sigma),line=0.5)

boxplot(beta0~Sim,data = box.data,notch = TRUE,outline=FALSE,names = c(expression(paste(alpha,'=8, ',beta,'=2')),
                                                                       expression(paste(alpha,'=9, ',beta,'=3')),
                                                                       expression(paste(alpha,'=9, ',beta,'=4'))))
title(expression(beta[0]),line=0.5)

boxplot(beta1~Sim,data = box.data,notch = TRUE,outline=FALSE,names = c(expression(paste(alpha,'=8, ',beta,'=2')),
                                                                       expression(paste(alpha,'=9, ',beta,'=3')),
                                                                       expression(paste(alpha,'=9, ',beta,'=4'))))
title(expression(beta[1]),line=0.5)

boxplot(beta2~Sim,data = box.data,notch = TRUE,outline=FALSE,names = c(expression(paste(alpha,'=8, ',beta,'=2')),
                                                                       expression(paste(alpha,'=9, ',beta,'=3')),
                                                                       expression(paste(alpha,'=9, ',beta,'=4'))))
title(expression(beta[2]),line=0.5)

boxplot(beta3~Sim,data = box.data,notch = TRUE,outline=FALSE,names = c(expression(paste(alpha,'=8, ',beta,'=2')),
                                                                       expression(paste(alpha,'=9, ',beta,'=3')),
                                                                       expression(paste(alpha,'=9, ',beta,'=4'))))
title(expression(beta[3]),line=0.5)

boxplot(beta4~Sim,data = box.data,notch = TRUE,outline=FALSE,names = c(expression(paste(alpha,'=8, ',beta,'=2')),
                                                                       expression(paste(alpha,'=9, ',beta,'=3')),
                                                                       expression(paste(alpha,'=9, ',beta,'=4'))))
title(expression(beta[4]),line=0.5)

boxplot(beta5~Sim,data = box.data,notch = TRUE,outline=FALSE,names = c(expression(paste(alpha,'=8, ',beta,'=2')),
                                                                       expression(paste(alpha,'=9, ',beta,'=3')),
                                                                       expression(paste(alpha,'=9, ',beta,'=4'))))
title(expression(beta[5]),line=0.5)

boxplot(beta6~Sim,data = box.data,notch = TRUE,outline=FALSE,names = c(expression(paste(alpha,'=8, ',beta,'=2')),
                                                                       expression(paste(alpha,'=9, ',beta,'=3')),
                                                                       expression(paste(alpha,'=9, ',beta,'=4'))))
title(expression(beta[6]),line=0.5)
boxplot(beta7~Sim,data = box.data,notch = TRUE,outline=FALSE,names = c(expression(paste(alpha,'=8, ',beta,'=2')),
                                                                       expression(paste(alpha,'=9, ',beta,'=3')),
                                                                       expression(paste(alpha,'=9, ',beta,'=4'))))
title(expression(beta[7]),line=0.5)


dev.off()
