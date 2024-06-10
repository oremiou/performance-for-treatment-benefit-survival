# graphical calibration with KM curves ---------------------------
m2 <- survfit(Surv(time, status) ~ tr, data=survdat) 
res <- summary(m2, censored = T)
dt1 <- with(res, data.frame(time = time, surv = surv, upper = upper,
                            lower = lower, strata=strata))

dt2 <- data.frame(time = c(1:300/100, 1:300/100, dt1$time), 
                  surv = c(predictions.both.groups.CV$means, dt1$surv), 
                  group = c(predictions.both.groups.CV$group, dt1$strata))
dt2$group[dt2$group == "1"] <- "Z=0, observed"
dt2$group[dt2$group == "2"] <- "Z=1, observed"

ggplot(data = dt2) + 
  geom_line(aes(x = time, y = surv, group = group, colour=group, linetype=group,
                size = group)) +
  labs(y="Survival probability", x="Time (years)")+
  scale_size_manual(values=c(0.8,0.1,0.8,0.1))+ 
  scale_color_manual(values=c("red", "red", "black", "black"))+
  scale_linetype_manual(values=c("solid", "dashed","solid", "dashed"))+
  xlim(0,3)+ theme(axis.text.y = element_text(size = 12))+ theme(axis.text.x = element_text(size = 12))+
  theme(axis.title = element_text(size = 14))



# model-based calibration ----------------------------
dat.cv.copy<-dat.CV
dat.cv.copy$tr<-as.numeric(dat.cv.copy$tr)-1
coxph(Surv(time, status)~ predicted.outcome.tr0.cll+tr+tr:benefit.survival.at.fixed.time.cll,x=T,
      data=dat.cv.copy)
confint(coxph(Surv(time, status)~ predicted.outcome.tr0.cll+tr+tr:benefit.survival.at.fixed.time.cll,x=T,
              data=dat.cv.copy))
library(Hmisc)
library(rms)
library(pec)
dat.cv.copy$tr<-as.factor(dat.cv.copy$tr)
calibrate.cox <- coxph(Surv(time, status)~ predicted.outcome.tr0.cll+#rcs(predicted.outcome.tr0.cll,3)+
                         tr+tr:rcs(benefit.survival.at.fixed.time.cll,5),x=T,
                       data=dat.cv.copy)


predict.calibrate.cox <- 1 - predictSurvProb(calibrate.cox,newdata=data.frame(predicted.outcome.tr0.cll=dat.cv.copy$predicted.outcome.tr0.cll, 
                                                                              benefit.survival.at.fixed.time.cll=dat.cv.copy$benefit.survival.at.fixed.time.cll,
                                                                              tr="0"), times=t.0)-
  (1 - predictSurvProb(calibrate.cox,newdata=data.frame(predicted.outcome.tr0.cll=dat.cv.copy$predicted.outcome.tr0.cll, 
                                                        benefit.survival.at.fixed.time.cll=dat.cv.copy$benefit.survival.at.fixed.time.cll,
                                                        tr="1"), times=t.0))

dd=data.frame(predicted=dat.CV$benefit.survival.at.fixed.time, observed=predict.calibrate.cox) 



ggplot(dd, aes(predicted, observed)) + 
  geom_smooth(method = "loess")+ geom_abline(intercept = 0, slope = 1)


# calibrate benefit in survival probability -----------
t.0=2
Ngroups=5
d1 <- quantile(dat.CV$benefit.survival.at.fixed.time, probs = seq(0, 1, 1/Ngroups))
g1 <- list()
for (i in 1:Ngroups) {
  g1[[i]] <- dat.CV[dat.CV$benefit.survival.at.fixed.time >= d1[i] & 
                      dat.CV$benefit.survival.at.fixed.time < d1[i + 1], ]}

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
p1=suppressMessages(ggplot(dat1,aes(x=pred,y=obs))+
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

# predict benefit in restricted mean survival time -------
Ngroups=5
t.star=2
pct=1:99/100

d1 <- quantile(dat.CV$predicted.benefit.RMS, probs = seq(0, 1, 1/Ngroups))
g1 <- list()
for (i in 1:Ngroups) {
  g1[[i]] <- dat.CV[dat.CV$predicted.benefit.RMS >= d1[i] & 
                      dat.CV$predicted.benefit.RMS < d1[i + 1], ]
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
  predicted <- c(predicted, mean(g1[[i]]$predicted.benefit.RMS))
  lower.predicted <- c(lower.predicted,min(g1[[i]]$predicted.benefit.RMS)  )
  upper.predicted <- c(upper.predicted,max(g1[[i]]$predicted.benefit.RMS)  )
}


dat1 <- data.frame(pred = predicted, obs = observed, 
                   lower.predicted , upper.predicted,
                   lower.observed , upper.observed)

library(ggplot2)
p2=suppressMessages(ggplot(dat1,aes(x=pred,y=obs))+
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

gridExtra::grid.arrange(p1,p2, nrow=1)
