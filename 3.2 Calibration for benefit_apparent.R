## calibration -----------------------------

# graphical calibration with KM curves
group0 <- survdat[survdat$tr == 0,] 
group1 <- survdat[survdat$tr == 1,]

prediction.weibull <-
  function(new.patient, fit.model = NULL){
    pct <- 1:300/100
    mu_hat <- predict(fit.model, newdata = new.patient, type = "lp")
    predicted.exp <- 1 - pweibull(pct, shape = 1/fit.model$scale, 
                                  scale = exp(mu_hat))
    results <- data.frame(time = pct, means = predicted.exp)
    return(results)
  }  

predictions.group0 <- list()  # predicted probabilities under control
for (i in 1:dim(group0)[1]){
  predictions.group0[[i]] <- 
    prediction.weibull(group0[i,],fit.model)$means
}

predictions.group1 <- list() # predicted probabilities under treatment
for (i in 1:dim(group1)[1]){
  predictions.group1[[i]] <- prediction.weibull(group1[i,],fit.model)$means
}


predictions.mean.group1 <- apply(do.call(cbind, predictions.group1), 1, mean)
predictions.mean.group0 <- apply(do.call(cbind, predictions.group0), 1, mean)


predictions.both.groups <- data.frame(means= c(predictions.mean.group0,predictions.mean.group1), 
                                      group= c(rep("Z=0, fitted", 300), rep("Z=1, fitted", 300)), percentiles.surv = 1:300/100)

m2 <- survfit(Surv(time, status) ~ tr, data=survdat) 
res <- summary(m2, censored = T)
dt1 <- with(res, data.frame(time = time, surv = surv, upper = upper,
                            lower = lower, strata=strata))

dt2 <- data.frame(time = c(1:300/100, 1:300/100, dt1$time), 
                  surv = c(predictions.both.groups$means, dt1$surv), 
                  group = c(predictions.both.groups$group, dt1$strata))
dt2$group[dt2$group == "1"] <- "Z=0, observed"
dt2$group[dt2$group == "2"] <- "Z=1, observed"

ggplot(data = dt2) + 
  geom_line(aes(x = time, y = surv, group = group, colour=group, linetype=group,
                size = group)) +
  labs(x="Percent survival", y="Time (years)")+
  scale_size_manual(values=c(0.8,0.1,0.8,0.1))+ 
  scale_color_manual(values=c("red", "red", "black", "black"))+
  scale_linetype_manual(values=c("solid", "dashed","solid", "dashed"))+
  xlim(0,3)+ theme(axis.text.y = element_text(size = 12))+ theme(axis.text.x = element_text(size = 12))+
  theme(axis.title = element_text(size = 14))


# model-based calibration
## set time
t.0=2
# assess calibration
group0 <- survdat; group0$tr="0"
for(i in 1: N){  group0$wb[i]=1-prediction.weibull(group0[i,], fit.model)[t.0*100,2] }
survdat$prob.t0.cll=log(-log(1-group0$wb)) # complimentary log-log transformation

for( i in 1:N){
  new.patient0 <- survdat[i,]; new.patient0$tr="0" 
  new.patient1 <- survdat[i,]; new.patient1$tr="1"
  mu_hat0 <- predict(fit.model, newdata=new.patient0, type="link")
  mu_hat1 <- predict(fit.model, newdata=new.patient1, type="link")
  
  survdat$benefit.survival.at.fixed.time[i] <- 
    1 - pweibull(t.0,shape=1/fit.model$scale,scale=exp(mu_hat1))-(
      1 - pweibull(t.0,shape=1/fit.model$scale,scale=exp(mu_hat0)))
  
  survdat$benefit.survival.at.fixed.time.cll[i] <- 
    log(-log( 1 - pweibull(t.0,shape=1/fit.model$scale,scale=exp(mu_hat1))))-
    log(-log(
      1 - pweibull(t.0,shape=1/fit.model$scale,scale=exp(mu_hat0))) )
}
survdat2=survdat
survdat2$tr=as.numeric(survdat2$tr)-1
cox.aux=coxph(Surv(time, status)~ prob.t0.cll+tr+tr:benefit.survival.at.fixed.time.cll,x=T,
      data=survdat2)
confint(cox.aux)
library(Hmisc)
library(rms)
library(pec)
calibrate.cox <- coxph(Surv(time, status)~ rcs(prob.t0.cll,3)+tr+tr:rcs(benefit.survival.at.fixed.time.cll,3),x=T,
                       data=survdat2)

predict.grid.wb <- seq(quantile(survdat$benefit.survival.at.fixed.time.cll,probs=0.01),
                       quantile(survdat$benefit.survival.at.fixed.time.cll,probs=0.99),length=100)

predict.calibrate.cox <- 1 - predictSurvProb(calibrate.cox,newdata=data.frame(prob.t0.cll=survdat2$prob.t0.cll, 
                                                                              benefit.survival.at.fixed.time.cll=survdat2$benefit.survival.at.fixed.time.cll,
                                                                              tr=0), times=t.0)-
  (1 - predictSurvProb(calibrate.cox,newdata=data.frame(prob.t0.cll=survdat2$prob.t0.cll, 
                                                        benefit.survival.at.fixed.time.cll=survdat2$benefit.survival.at.fixed.time.cll,
                                                        tr=1), times=t.0))

dd=data.frame(predicted=survdat2$benefit.survival.at.fixed.time, observed=predict.calibrate.cox) 



ggplot(dd, aes(predicted, observed)) + geom_smooth(method = "loess")+ geom_abline(intercept = 0, slope = 1)





# predict benefit in survival probability
t.0=2
Ngroups=5
survdat$benefit.survival.at.fixed.time<- -100
for( i in 1:N)
{
  new.patient0 <- survdat[i,]; new.patient0$tr="0" 
  new.patient1 <- survdat[i,]; new.patient1$tr="1"
  mu_hat0 <- predict(fit.model, newdata=new.patient0, type="link")
  mu_hat1 <- predict(fit.model, newdata=new.patient1, type="link")
  survdat$benefit.survival.at.fixed.time[i] <- 
    1 - pweibull(t.0,shape=1/fit.model$scale,scale=exp(mu_hat1))-(
      1 - pweibull(t.0,shape=1/fit.model$scale,scale=exp(mu_hat0))) }

d1 <- quantile(survdat$benefit.survival.at.fixed.time, probs = seq(0, 1, 1/Ngroups))
g1 <- list()
for (i in 1:Ngroups) {
  g1[[i]] <- survdat[survdat$benefit.survival.at.fixed.time >= d1[i] & 
                       survdat$benefit.survival.at.fixed.time < d1[i + 1], ]}

hist(survdat$benefit.survival.at.fixed.time)
predicted <- c()
observed <- c()
SE.observed <- c()
lower.predicted<-c()
upper.predicted<-c()

library(survRM2 )
for (i in 1:Ngroups){
  km0 <- survfit(Surv(time, status)~tr, data=g1[[i]][g1[[i]]$tr==0,])
  km1 <- survfit(Surv(time, status)~tr, data=g1[[i]][g1[[i]]$tr==1,])
  observed <- c(observed,summary(km1, times=t.0)[6][[1]] - summary(km0, times=t.0)[6][[1]])
  SE.observed <- c(SE.observed, sqrt(summary(km0, times=t.0)[7][[1]]^2+summary(km1, times=t.0)[7][[1]]^2))
  predicted <- c(predicted, mean(g1[[i]]$benefit.survival.at.fixed.time))
  lower.predicted <- c(lower.predicted,min(g1[[i]]$benefit.survival.at.fixed.time)  )
  upper.predicted <- c(upper.predicted,max(g1[[i]]$benefit.survival.at.fixed.time)  )
}


dat1 <- data.frame(pred = predicted, obs = observed,
                   SE.observed ,  lower.predicted , upper.predicted)

library(ggplot2)
suppressMessages(ggplot(dat1,aes(x=pred,y=obs))+
                   geom_point(size=3,shape=20)+
                   labs(x="Predicted benefit (survival probability at fixed time)", y="Observed benefit (survival probability at fixed time)")+
                   geom_abline(intercept=0,slope=1,color="black",linetype="dashed",
                               size=0.5)+
                   geom_errorbar(aes(ymin=obs-1.96*SE.observed,ymax=obs+1.96*SE.observed,),width=0.01)+
                   geom_errorbarh(aes(xmin=lower.predicted,xmax=upper.predicted), height=0.01)+
                   geom_smooth(method="lm",
                               colour="blue",size=0.2)+theme(aspect.ratio=1))

#slope
lm(dat1$obs~dat1$pred)
confint(lm(dat1$obs~dat1$pred))
# MSE 
sqrt(mean((dat1$pred-dat1$obs)^2))
# mean bias 
(mean((dat1$pred-dat1$obs)))

# predict benefit in restricted mean survival time
Ngroups=5
t.star=2
pct=1:99/100
survdat$predicted.benefit.RMSE3=-100
for (i in 1:N){
  new.patient0 <- survdat[i,]; new.patient0$tr="0" 
  new.patient1 <- survdat[i,]; new.patient1$tr="1"
  ptime.mean0 <- predict(fit.model, new.patient0, 
                         type='quantile', p=pct, se=TRUE)$fit
  ptime.mean1 <- predict(fit.model, new.patient1, 
                         type='quantile', p=pct, se=TRUE)$fit
  
  ptime.mean1[ptime.mean1>t.star]=t.star
  ptime.mean0[ptime.mean0>t.star]=t.star
  survdat$predicted.benefit.RMSE3[i]=sum(ptime.mean1 *0.01)-sum(ptime.mean0 *0.01)
}
summary(survdat$predicted.benefit.RMSE3)

d1 <- quantile(survdat$predicted.benefit.RMSE3, probs = seq(0, 1, 1/Ngroups))
g1 <- list()
for (i in 1:Ngroups) {
  g1[[i]] <- survdat[survdat$predicted.benefit.RMSE3 >= d1[i] & 
                       survdat$predicted.benefit.RMSE3 < d1[i + 1], ]
}
predicted <- c()
observed <- c()
lower.observed <- c()
upper.observed <- c()
lower.predicted<-c()
upper.predicted<-c()

library(survRM2 )
for (i in 1:Ngroups){
  rmst2.fit <- rmst2(g1[[i]]$time, g1[[i]]$status, g1[[i]]$tr, tau=t.star)
  
  observed <- c(observed,rmst2.fit$unadjusted.result[1,1])
  lower.observed <- c(lower.observed, rmst2.fit$unadjusted.result[1,2])
  upper.observed <- c(upper.observed, rmst2.fit$unadjusted.result[1,3])
  predicted <- c(predicted, mean(g1[[i]]$predicted.benefit.RMSE3))
  lower.predicted <- c(lower.predicted,min(g1[[i]]$predicted.benefit.RMSE3)  )
  upper.predicted <- c(upper.predicted,max(g1[[i]]$predicted.benefit.RMSE3)  )
}


dat1 <- data.frame(pred = predicted, obs = observed, 
                   lower.predicted , upper.predicted,
                   lower.observed , upper.observed)

library(ggplot2)
suppressMessages(ggplot(dat1,aes(x=pred,y=obs))+
                   geom_point(size=3,shape=20)+
                   labs(x="Predicted benefit (RMST)", y="Observed benefit (RMST)")+
                   geom_abline(intercept=0,slope=1,color="black",linetype="dashed",
                               size=0.5)+
                   geom_errorbar(aes(ymin=lower.observed,
                                     ymax=upper.observed),width=0.01)+
                   geom_errorbarh(aes(xmin=lower.predicted,xmax=upper.predicted), height=0.01)+
                   geom_smooth(method="lm",
                               colour="blue",size=0.2)+theme(aspect.ratio=1))

#slope
lm(dat1$obs~dat1$pred)
confint(lm(dat1$obs~dat1$pred))
# RMSE 
sqrt(mean((dat1$pred-dat1$obs)^2))
# mean bias 
(mean((dat1$pred-dat1$obs)))

