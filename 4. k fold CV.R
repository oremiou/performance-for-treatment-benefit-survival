
N.folds<-10
# Randomly shuffle the data
survdat1<-survdat[sample(nrow(survdat)),]
# Create random folds
survdat1$folds <- cut(seq(1,nrow(survdat1)),breaks=N.folds,labels=FALSE) # 10-fold CV
t.0 <- 2 # assess measures at 2 years
t.star<-2 # for RMS
pct <- 1:99/100 # percentiles


prediction.weibull <-
  function(new.patient, fit.model = NULL){
    pct1 <- 1:300/100
    mu_hat <- predict(fit.model, newdata = new.patient, type = "lp")
    predicted.exp <- 1 - pweibull(pct1, shape = 1/fit.model$scale, 
                                  scale = exp(mu_hat))
    results <- data.frame(time = pct1, means = predicted.exp)
    return(results)
  }  
predictions.both.groups=list()
predictions.group0 <- list() 
predictions.group1 <- list() 

# Obtain out-of-sample predictions 
dat.out.CV<-list()
for (i in 1:N.folds){
  dat.in.CV=survdat1[survdat1$folds!=i,]
  dat.out.CV[[i]]=survdat1[survdat1$folds==i,]
  dat1<-dat.out.CV[[i]]; dat1$t=1
  dat0<-dat.out.CV[[i]]; dat0$t=0
  fit.model.1 <- survreg(Surv(time, status)~x1+x2+x3+x4+tr+x1*tr+x2*tr+x3*tr+x4*tr, 
                         dat.in.CV, dist="weibul")
  
  # benefit in survival at fixed time point
  for( k in 1:length(dat.out.CV[[i]]$x1)) # cycle through all patients in left out fold
  {
    new.patient0 <- dat.out.CV[[i]][k,]; new.patient0$tr="0"
    new.patient1 <- dat.out.CV[[i]][k,]; new.patient1$tr="1"
    mu_hat0 <- predict(fit.model.1, newdata=new.patient0, type="link")
    mu_hat1 <- predict(fit.model.1, newdata=new.patient1, type="link")
    
    dat.out.CV[[i]]$predicted.survival.tr1[k] <- # predicted survival at t.0 under tr=1
      1 - pweibull(t.0,shape=1/fit.model.1$scale,scale=exp(mu_hat1))
    
    dat.out.CV[[i]]$predicted.survival.tr0[k] <- # predicted survival at t.0 under tr=0
      1 - pweibull(t.0,shape=1/fit.model.1$scale,scale=exp(mu_hat0))
    
    dat.out.CV[[i]]$benefit.survival.at.fixed.time[k] <- # benefit survival at t.0 
      dat.out.CV[[i]]$predicted.survival.tr1[k]- 
      dat.out.CV[[i]]$predicted.survival.tr0[k]

    dat.out.CV[[i]]$predicted.outcome.tr0[k]<-1- # predicted event probability at t.0 under tr=0
      dat.out.CV[[i]]$predicted.survival.tr0[k]
    dat.out.CV[[i]]$predicted.outcome.tr0.cll[k]<-
      log(-log(1-dat.out.CV[[i]]$predicted.outcome.tr0[k]))# log-log transformed prob. at t.0 under tr=0
  
    dat.out.CV[[i]]$benefit.survival.at.fixed.time.cll[k] <- # log-log transformed benefit at t.0
      log(-log(  dat.out.CV[[i]]$predicted.survival.tr1[k]))-
      log(-log(  dat.out.CV[[i]]$predicted.survival.tr0[k]))
    
    
    ptime.mean0 <- predict(fit.model.1, new.patient0, 
                           type='quantile', p=pct, se=TRUE)$fit
    ptime.mean1 <- predict(fit.model.1, new.patient1, 
                           type='quantile', p=pct, se=TRUE)$fit
    
    ptime.mean1[ptime.mean1>t.star]=t.star
    ptime.mean0[ptime.mean0>t.star]=t.star
    dat.out.CV[[i]]$predicted.benefit.RMS[k]=sum(ptime.mean1 *0.01)-sum(ptime.mean0 *0.01) #   predicted benefit in RMS
    
  }
  
  group0 <- dat.out.CV[[i]][dat.out.CV[[i]]$tr == 0,] 
  group1 <- dat.out.CV[[i]][dat.out.CV[[i]]$tr == 1,]
  
 # predicted probabilities under control
  for (j in 1:dim(group0)[1]){
    predictions.group0=append( predictions.group0, 
      list(c(prediction.weibull(group0[j,],fit.model.1)$means)))
  }
 
  
# predicted probabilities under treatment
  for (j in 1:dim(group1)[1]){
    predictions.group1=append( predictions.group1, 
            list(c(prediction.weibull(group1[j,],fit.model.1)$means)))
  }
  
  }

dat.CV<-dat.out.CV[[1]]
for (i in 2:N.folds){  dat.CV<-rbind(dat.CV,dat.out.CV[[i]])}

predictions.mean.group1 <- apply(do.call(cbind, predictions.group1), 1, mean)
predictions.mean.group0 <- apply(do.call(cbind, predictions.group0), 1, mean)


predictions.both.groups.CV<- data.frame(means= c(predictions.mean.group0,predictions.mean.group1), 
                                           group= c(rep("Z=0, fitted", 300), rep("Z=1, fitted", 300)), percentiles.surv = 1:300/100)
