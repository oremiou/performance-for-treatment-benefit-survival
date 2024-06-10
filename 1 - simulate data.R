remove(list=ls())
set.seed(42) # the answer to life the universe and everything
options(digits=2)

# simulate data ----
library(MASS)
N <- 2000
Sigma <- outer(1:4, 1:4, function(x,y) 0.5^abs(x-y)) #variance covariance matrix for covariates

x <- mvrnorm(n = N, rep(0, 4), Sigma)
x[,3] <- cut(x[,3], breaks=c(-Inf, -1, 1, 2, Inf)) #categorical with 4 categories
x[,4] <- ifelse(x[,4] > 0, 1, 0) #binary predictor

tr <- rbinom(size=1,n=N, p=0.5)
survdat <- data.frame(x,tr)
survdat[,3:5] <- lapply(survdat[,3:5], factor)
colnames(survdat)[1:4] <- paste0("x", 1:4)

survdat$rate0 <- with(survdat, 
                      exp(-1+0.3*x1+0.3*x2+0.1*(x3=="2")
                          -0.2*(x3=="3")+0.1*(x3=="4")
                          -0.1*(x4=="1")))
survdat$benefit <- with(survdat, 
                        exp(-0.1-0.1*x1-0.2*x2+0.2*(x3=="1")-0.2*(x4=="1"))) # treatment benefit
survdat$rate1 <- with(survdat, rate0*benefit)

survdat$survivaltime0 <- rexp(N, rate = survdat$rate0)
survdat$survivaltime1 <- rexp(N, rate = survdat$rate1)

survdat$survival.with.treat.taken=with(survdat, survivaltime0*(tr=="0")+survivaltime1*(tr=="1"))

censtimes <- 1 + 2*runif(N)
mean(censtimes)
survdat$time <- pmin(survdat$survival.with.treat.taken, censtimes)
survdat$status <- as.numeric(censtimes > survdat$survival.with.treat.taken)
table(survdat$status)
head(survdat)
median(survdat$time)

sum(survdat$benefit<1)
# visualize data ----
library(ggplot2)
ggplot(survdat, aes(x=survival.with.treat.taken)) + geom_histogram(color="black", fill="white", bins=50)
median(survdat$survivaltime0)
median(survdat$survivaltime1)
sum(survdat$status[survdat$tr==0]); sum(survdat$tr==0); sum(survdat$status[survdat$tr==0])/ sum(survdat$tr==0);median(survdat$time[survdat$tr==0])
sum(survdat$status[survdat$tr==1]); sum(survdat$tr==1); sum(survdat$status[survdat$tr==1])/ sum(survdat$tr==1);median(survdat$time[survdat$tr==1])

# KM curves
library(ggfortify)
library(survival)
model_fit <- survfit(Surv(time, status) ~ 1, data=survdat) 

autoplot(model_fit) + 
  labs(x = "\n Survival Time (years) ", y = "Survival Probabilities \n", 
       title = "Survival Times") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", size = 12),
        axis.title.y = element_text(face="bold", size = 12),
        legend.title = element_text(face="bold", size = 10))


### plot KM by treatment group
model_fit2 <- survfit(Surv(time, status) ~ tr, data=survdat) 
autoplot(model_fit2) + 
  labs(x = "\n Survival Time (Years) ", y = "Survival Probabilities \n", 
       title = "Survival Times") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", size = 12),
        axis.title.y = element_text(face="bold", size = 12),
        legend.title = element_text(face="bold", size = 10))+ xlim(0,3)

