df <- read.csv("clean.csv")
attach(df)
mod1 <- glm(SEI~School.Type+Diversity+Choice.Magnet.School+X2018..Distinctions+AverageDistance+NumberOfRail,family=binomial(link='logit'),data=df)
summary(mod1)
x<-model.matrix(mod1)

# basic function to fit the g-prior described in the paper
# x=design matrix without an intercept; intercept added automatically
# y=vector of Bernoulli responses, b and g as in paper
#gprior=function(x,y,b,g){
#  start=4*lm(y~x)$coef # crude least-squares starting values
#  n=length(y); x=cbind(rep(1,n),x); p=length(x[1,]); xtx=t(x)%*%x
#  bm=rep(0,p); bm[1]=b
#  ll=function(beta){
#    p=exp(x%*%beta)/(1+exp(x%*%beta))
#    -sum(x%*%beta*y-log(1+exp(x%*%beta)))+0.5*(beta-bm)%*%xtx%*%(beta-bm)/(g*n)}
#  fit=optim(start,ll,hessian=T)
#  cov=solve(fit$hessian)
# cat("Logistic Regression with informative g-prior \n")
# cat("b =",b," g =",g,"\n")
#  cat("Parameter PostMode PostSD \n")
#  cat("------------------------- \n")
#  for(i in 1:p){cat("beta[",i-1,"]",fit$par[i]," ",sqrt(cov[i,i]),"\n")}
#  cat("\n")
  
#}


x<-model.matrix(mod1)   # Gives design matrix of logistic model
cof<-as.vector(mod1$coefficients)  # Coefficient vector
g <- 149   # Set number of observations as g
sigma <- 2
dhyperg(g, alpha=3)
dzellner(cof, g, sigma, x)
rzellner(1, g, sigma, x)  


