
# mean bias in survival, two years 
km0 <- survfit(Surv(time, status)~tr, data=dat.CV[dat.CV$tr==0,])
km1 <- survfit(Surv(time, status)~tr, data=dat.CV[dat.CV$tr==1,])

summary(km1, times=t.0)[6][[1]] - summary(km0, times=t.0)[6][[1]]-
  mean(dat.CV$benefit.survival.at.fixed.time)


# mean bias in RMST, two years
library(survRM2 )
t.0=2
rmst.fit <- rmst2(dat.CV$time,dat.CV$status,dat.CV$tr, tau=t.0)

observed <- rmst.fit$unadjusted.result[1,1]
predicted <- mean(dat.CV$predicted.benefit.RMSE3)
observed-predicted
  

# estimate c for benefit 
group0 <- dat.CV[dat.CV$tr == 0,c("time", "status", "benefit.survival.at.fixed.time")]
group1 <- dat.CV[dat.CV$tr == 1,c("time", "status", "benefit.survival.at.fixed.time")]

n0<-dim(group0)[1]
n1<-dim(group1)[1]

c.benefit=c()
c.benefit.se=c()
repeats=200

for (rep in 1:repeats){
  if(n0>n1){
    group00=group0[sample(1:n0,n1, replace = F),]
    group11=group1
  }  else  {
    group00=group0
    group11=group1[sample(1:n1,n0, replace = F),]}
  
  group00 <- group00[order(group00$benefit.survival.at.fixed.time),]
  group11 <- group11[order(group11$benefit.survival.at.fixed.time),]
  
  colnames(group11)=c("time1", "status1", "benefit.survival.at.fixed.time1")
  
  dat.all=cbind(group00,group11)
  
  dat.all$pred.ben.avg <- with(dat.all, 
                               (benefit.survival.at.fixed.time+benefit.survival.at.fixed.time1)/2)>0
  
  dat.all$observed.ben <-0
  dat.all$status.at.smallest.time<-with(dat.all, status*(time<time1)+status1*(time1<time))
  dat.all$observed.ben<- dat.all$status.at.smallest.time*(2*(dat.all$time1>dat.all$time)-1)
  
  ## removing 0 benefit
  dat.all1=dat.all[dat.all$observed.ben!=0,]
# Benefit c-statistic
  library(Hmisc)
  cindex <- rcorr.cens(dat.all1$pred.ben.avg, dat.all1$observed.ben)
  c.benefit <- c(c.benefit, cindex["C Index"][[1]])
  c.benefit.se <- c(c.benefit.se, cindex["S.D."][[1]]/2	)# Half the sd of Dxy
}

mean(c.benefit)
mean(c.benefit - 1.96*c.benefit.se)
mean(c.benefit + 1.96*c.benefit.se)		




# match by X
library(Matching)
c.by.covariates=c(); se.c.by.covariates=c()
repeats=100
X=dat.CV[,c("x1", "x2","x3", "x4")]
c.by.covariates=c(); se.c.by.covariates=c()
n.0 <- sum(survdat$tr==0)
n.1 <- sum(survdat$tr==1)
n.10=min(n.0, n.1)
dat.CV3=dat.CV
for (i in 1:repeats){
  rows.in<-c(which(dat.CV$tr==1)[sample(n.1, n.10)], which(dat.CV$tr==0)[sample(n.0, n.10)])
  dat.CV4=dat.CV3[rows.in[order(rows.in)],]
  X2=X[rows.in[order(rows.in)],]
  
  X2$x3.2=X2$x3=="2"
  X2$x3.3=X2$x3=="3"
  X2$x3.4=X2$x3=="4"
  X2$x4=as.numeric(X2$x4)-1
  X2=X2[,c("x1","x2", "x3.2", "x3.3", "x3.4")]
  
  rr1 <- Match(Tr=as.numeric(dat.CV4$tr)-1 ,X=X2, M=1,ties=F,replace=FALSE)
  ind.0 <- rr1$index.control
  ind.1 <- rr1$index.treated
  ### Calculation of predicted and observed benefit in matched pairs
  dat.all=cbind(dat.CV4[ind.0,c("time", "status","benefit.survival.at.fixed.time")],
                dat.CV4[ind.1,c("time", "status","benefit.survival.at.fixed.time")])
  colnames(dat.all)=c("time", "status", "benefit.survival.at.fixed.time", "time1", 
                      "status1", "benefit.survival.at.fixed.time1")
  dat.all$pred.ben.avg <- with(dat.all, 
                               (benefit.survival.at.fixed.time+benefit.survival.at.fixed.time1)/2)>0
  
  dat.all$observed.ben <-0
  dat.all$status.at.smallest.time<-with(dat.all, status*(time<time1)+status1*(time1<time))
  dat.all$observed.ben<- dat.all$status.at.smallest.time*(2*(dat.all$time1>dat.all$time)-1)
  
    ## removing 0 benefit
  dat.all1=dat.all[dat.all$observed.ben!=0,]
  # Benefit c-statistic
  library(Hmisc)
  cindex <- rcorr.cens(dat.all1$pred.ben.avg, dat.all1$observed.ben)
  c.by.covariates <- c(c.by.covariates, cindex["C Index"][[1]])
  se.c.by.covariates <- c(se.c.by.covariates, cindex["S.D."][[1]]/2	)# Half the sd of Dxy
}

mean(c.by.covariates)
mean(c.by.covariates - 1.96*se.c.by.covariates)
mean(c.by.covariates + 1.96*se.c.by.covariates)		



