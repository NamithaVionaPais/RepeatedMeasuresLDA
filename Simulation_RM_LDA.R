setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/Thesis_NVP/Chapter2")
library(rstan)
library(gdata)
library(bayesplot) 
library(gtools)
library(MASS)
library(matrixcalc)
library(rBeta2009)
library(cmdstanr)
library(invgamma)
library(trialr)


set.seed(221196)
options(scipen = 999) # disabling scientific notation
J <-3  #number of topics/motivation
V <-6 #number of words/products
M <-100 #number of documents
n_d<-10
n_d1<-2
N<-300 #words in each sub doc n 

ytrue<-matrix(,nrow=M,ncol=n_d)
ytest<-matrix(,nrow=M,ncol=n_d1)
bp <- rep(200/V,V) 
beta <- array(NA,c(J,V)) 
beta[1,]<-c(0.1,0.6,0.1,0.1,0.05,0.05)
beta[2,]<-c(0.01,0.05,0.1,0.6,0.01,0.24)
beta[3,]<-c(0.2,0.23,0.5,0.05,0.01,0.01)

w <- rep(NA,N) #represents a particular word n
doc <- rep(NA,N)#represents a particular document id for word n 


#Parameters
sigma<-list()
mu<-list()
for(d in 1:M)
{
  sigma[[d]]<-rlkjcorr(1, J, eta = 2)
  mu[[d]]<-rep(0,J)
}


#glmer---rho
alpha<-matrix(,nrow=M,ncol=n_d)
w<-array(, dim = c(M, n_d, N))
wtest<-array(,dim = c(M, n_d1, N))

for(d in 1:M)
{
  la1<- mvrnorm(n = 1,mu =mu[[d]],sigma[[d]])
  for(r in 1:n_d)
  {
    theta<-rdirichlet(1, exp(la1))
    zbar<-c(0,0,0)
    for(n in 1:N)
    {
      z <- which(rmultinom(1,1,theta) == 1) 
      zbar[z]<-zbar[z]+1
      w[d,r,n] <- which(rmultinom(1,1,beta[z,]) == 1)
    }
    zbar<-zbar/N
    ytrue[d,r]<-which.max(zbar)
  }
  for(r in 1:n_d1)
  {
    theta<-rdirichlet(1, exp(la1))
    zbar<-c(0,0,0)
    for(n in 1:N)
    {
      z <- which(rmultinom(1,1,theta) == 1) 
      zbar[z]<-zbar[z]+1
      wtest[d,r,n] <- which(rmultinom(1,1,beta[z,]) == 1)
    }
    zbar<-zbar/N
    ytest[d,r]<-which.max(zbar)
  }
}

data <- list(J=J, V=V, M=M, N=N, n_d=n_d,w=w,mu=rep(0,J))
data


model<- cmdstan_model(stan_file = "RM_LDA.stan")
model$print()
fit <- model$variational(data = data,
                         seed=1,output_samples = 1000,
                         eval_elbo=20,grad_samples =20,
                         elbo_samples = 20,
                         algorithm = "meanfield",
                         output_dir = NULL,
                         iter = 1000,
                         adapt_iter = 500,
                         save_latent_dynamics=TRUE,
                         tol_rel_obj = 10^-3)


library(rstan)


fitmcmc<-model$sample(data = data,
                      seed = 1,
                      chains = 1,
                      iter_warmup = 10,
                      iter_sampling = 100)

b1<-matrix(fit$summary("beta")$mean,nrow=3,byrow=FALSE);b1
beta

b1se<-matrix(fit$summary("beta")$sd,nrow=3,byrow=FALSE);b1se

databeta<-cbind(fit$summary("beta")$mean,fit$summary("beta")$sd)
colnames(databeta)<-c("mean","sd")
plot <- ggplot(databeta, aes( y = mean)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5), size = 3) +
  labs(x = "Parameter", y = "Estimate", title = paste("Matrix", matrix_idx, "Parameter Estimates with SE")) +
  theme_minimal()



beta1<-b1
beta1[1,]<-b1[3,]
beta1[2,]<-b1[2,]
beta1[3,]<-b1[1,]
library(SMFilter)
round(FDist2(beta,beta1),4)

sigmapred<-array(fit$summary("sigma")$mean,dim=c(M,J,J))
sigma1<-list()
for(i in 1:M)
{
  sigma1[[i]]<-sigmapred[i,,]
  sigma1[[i]][1,2]=sigma1[[i]][2,1]=sigmapred[i,,][2,3]
  sigma1[[i]][2,3]=sigma1[[i]][3,2]=sigmapred[i,,][1,2]
  sigma1[[i]][1,3]=sigma1[[i]][3,1]=sigmapred[i,,][1,3]
  
}

fnorm<-c()
for(i in 1:M)
{
  fnorm[i]<-round(FDist2(sigma[[i]],sigma1[[i]]),4)
}
mean(fnorm)

vb_diag <- utils::read.csv(fit$latent_dynamics_files()[1], comment.char = "#")
ELBO <- data.frame(Iteration = vb_diag[,1], ELBO = vb_diag[,3])
ggplot(data = ELBO, aes(x = Iteration, y = ELBO)) + geom_line(lwd=1.5) + 
  theme(text = element_text(size = 20),
        panel.background = element_rect(fill = "transparent", color = "lightgrey"),
        panel.grid.major = element_line(colour = "lightgrey")) + xlim(0,450)
b1<-matrix(fit$summary("beta")$mean,nrow=3,byrow=FALSE);b1


ypred<-matrix(,nrow=M,ncol=n_d)
predtheta<-array(fit$summary("theta")$mean,dim=c(M,n_d,J))
for(d in 1:M)
{
  for( r in 1:n_d)
  {
    ypred[d,r]<-which.max(predtheta[d,r,])
  }
}



d1<-as.factor(as.vector(ytrue))
d2<-as.factor(as.vector(ypred))
d3<-c()

for(i in 1:length(ypred))
{
  if(d2[i]==1)
  {
    d3[i]=3
  }else if(d2[i]==2){
    d3[i]=2
  }else if(d2[i]==3){
    d3[i]=1
  }
}

y_true<-d1
y_pred<-d3
tab<-table(y_pred,y_true);tab
sum(diag(tab))/sum(tab)



beta
beta1
topics <- c("Topic 1", "Topic 2", "Topic 3")

# Convert data to long format
df <- data.frame(topic = rep(topics, each = ncol(beta1)),
                 Word = rep(c("Word 1", "Word 2", "Word 3", "Word 4", "Word 5", "Word 6"), times = nrow(beta1)),
                 proportion = as.vector(t(beta1)))

# Barplot with proportions
ggplot(df, aes(x = topic, y = proportion, fill = Word)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +ylim(0,1)+
  geom_text(aes(label = round(proportion, 2)), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  labs(x = "Topics", y = "Proportions", title = "Distribution of words over J topics") +
  scale_fill_manual(values = c("#FFD9D9", "#D9FFD9", "#D9D9FF", "#FFFFCC", "#FFCCFF", "#CCFFFF")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))



###TEST setup1
data1<- list(J=J, V=V, M=M, N=N, n_d=n_d1,w=wtest,mu=rep(0,J),beta=round(beta1,5),sigma=sigma1)
data1

model1<- cmdstan_model(stan_file = "RM_LDA_inf1.stan")
set.seed(3)
fit1 <- model1$variational(data = data1,
                           seed=3,output_samples = 1000,
                           eval_elbo=20,grad_samples =20,
                           elbo_samples = 20,
                           algorithm = "meanfield",
                           output_dir = NULL,
                           iter = 1000,
                           adapt_iter = 100,
                           save_latent_dynamics=TRUE,
                           tol_rel_obj = 10^-3)


ypred1<-matrix(,nrow=M,ncol=n_d1)
predtheta1<-array(fit1$summary("theta")$mean,dim=c(M,n_d1,J))
for(d in 1:M)
{
  for( r in 1:n_d1)
  {
    ypred1[d,r]<-which.max(predtheta1[d,r,])
  }
}

d1<-as.factor(as.vector(ytest))
d2<-as.factor(as.vector(ypred1))

y_true<-d1
y_pred<-d2
tab<-table(y_pred,y_true);tab
sum(diag(tab))/sum(tab)



###TEST setup 2

set.seed(65)
J <-3  #number of topics/motivation
V <-6 #number of words/products
M <-2 #number of documents
n_d<-50
N<-300 #words in each sub doc n 

ytrue<-matrix(,nrow=M,ncol=n_d)

w <- rep(NA,N) #represents a particular word n
doc <- rep(NA,N)#represents a particular document id for word n 


#Parameters
sigma<-list()
mu<-list()
for(d in 1:M)
{
  sigma[[d]]<-rlkjcorr(1, J, eta = 1)
  mu[[d]]<-rep(0,J)
}


#glmer---rho
alpha<-matrix(,nrow=M,ncol=n_d)
w<-array(, dim = c(M, n_d, N))
wtest<-array(,dim = c(M, n_d1, N))

for(d in 1:M)
{
  la1<- mvrnorm(n = 1,mu =mu[[d]],sigma[[d]])
  for(r in 1:n_d)
  {
    theta<-rdirichlet(1, exp(la1))
    zbar<-c(0,0,0)
    for(n in 1:N)
    {
      z <- which(rmultinom(1,1,theta) == 1) 
      zbar[z]<-zbar[z]+1
      w[d,r,n] <- which(rmultinom(1,1,beta1[z,]) == 1)
    }
    zbar<-zbar/N
    ytrue[d,r]<-which.max(zbar)
  }
}
summary(as.factor(ytrue))

data2 <- list(J=J, V=V, M=M, N=N, n_d=n_d,w=w,mu=rep(0,J),beta=round(beta1,5))
data2

model2<- cmdstan_model(stan_file = "RM_LDA_inf2.stan")

set.seed(3)
fit2 <- model2$variational(data = data2,
                           seed=1,output_samples = 1000,
                           eval_elbo=20,grad_samples =20,
                           elbo_samples = 20,
                           algorithm = "meanfield",
                           output_dir = NULL,
                           iter = 1000,
                           adapt_iter = 20,
                           save_latent_dynamics=TRUE,
                           tol_rel_obj = 10^-3)


ypred2<-matrix(,nrow=M,ncol=n_d)
predtheta<-array(fit2$summary("theta")$mean,dim=c(M,n_d,J))
for(d in 1:M)
{
  for( r in 1:n_d)
  {
    ypred2[d,r]<-which.max(predtheta[d,r,])
  }
}

d1<-as.factor(as.vector(ytrue))
d2<-as.factor(as.vector(ypred2))
y_true<-d1
y_pred<-d2
tab<-table(y_pred,y_true);tab
sum(diag(tab))/sum(tab)


fit2$summary("sigma")


