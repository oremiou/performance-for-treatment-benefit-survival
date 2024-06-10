# time of interest
t.0=2

# split by predicted survival
survdat$agreed <- as.numeric((survdat$tr==1)==(survdat$benefit.survival.at.fixed.time>0))
table(survdat$agreed)
table(survdat$tr==0, survdat$benefit.survival.at.fixed.time>0)
km23 <- survfit(Surv(time, status)~1, data=survdat[survdat$agreed==0,]) ## group 2U3, see main paper
km14 <- survfit(Surv(time, status)~1, data=survdat[survdat$agreed==1,]) ## group 1U4, see main paper
km24 <- survfit(Surv(time, status)~1, data=survdat[survdat$tr==0,]) ## group 2U4, see main paper
km13 <- survfit(Surv(time, status)~1, data=survdat[survdat$tr==1,]) ## group 1U3, see main paper
km1 <- survfit(Surv(time, status)~1, data=survdat[survdat$tr==1&survdat$benefit.survival.at.fixed.time>0,]) ## group 1, see main paper

### benefit for survival probability --------
## PB
PB <- summary(km14, times=t.0)[6][[1]]-summary(km23, times=t.0)[6][[1]]
SE.PB <- sqrt(summary(km23, times=t.0)[7][[1]]^2+summary(km14, times=t.0)[7][[1]]^2)
lower.PB <- PB-1.96*SE.PB
upper.PB <- PB+1.96*SE.PB
paste(round(PB, 3), " [", round(lower.PB,3)," ," ,round(upper.PB,3),"]", sep="")

# PB of treat all
PB.treat.all <- summary(km13, times=t.0)[6][[1]]-summary(km24, times=t.0)[6][[1]]
SE.PB.treat.all <- sqrt(summary(km13, times=t.0)[7][[1]]^2+summary(km24, times=t.0)[7][[1]]^2)
lower.PB.treat.all <- PB.treat.all-1.96*SE.PB.treat.all
upper.PB.treat.all <- PB.treat.all+1.96*SE.PB.treat.all
paste(round(PB.treat.all, 3), " [", round(lower.PB.treat.all,3)," ," ,round(upper.PB.treat.all,3),"]", sep="")

# estimate PB^(1) benefit of model vs. treating everybody
p.neg=mean(survdat$benefit.survival.at.fixed.time<0)
v.p.neg=p.neg*(1-p.neg)/N

km3 <- survfit(Surv(time, status)~1, data=survdat[survdat$tr==1&survdat$benefit.survival.at.fixed.time<0,]) ## group 3, see main paper
km4 <- survfit(Surv(time, status)~1, data=survdat[survdat$tr==0&survdat$benefit.survival.at.fixed.time<0,]) ## group 4, see main paper

PB.1=p.neg*(summary(km4, times=t.0)[6][[1]]-summary(km3, times=t.0)[6][[1]])
v43=summary(km4, times=t.0)[7][[1]]^2+summary(km3, times=t.0)[7][[1]]^2
se.PB.1=sqrt(
  v.p.neg*v43+
    v43*(p.neg)^2+
    v.p.neg*(summary(km4, times=t.0)[6][[1]]-summary(km3, times=t.0)[6][[1]])^2
)
lower.PB.1 <- PB.1-1.96*se.PB.1
upper.PB.1 <- PB.1+1.96*se.PB.1
paste(round(PB.1, 3), " [", round(lower.PB.1,3)," ," ,round(upper.PB.1,3),"]", sep="")


# estimate PB^(0) benefit of model vs. treating no one
p.pos=1-p.neg
v.p.pos=v.p.neg
km1 <- survfit(Surv(time, status)~1, data=survdat[survdat$tr==1&survdat$benefit.survival.at.fixed.time>0,]) ## group 1, see main paper
km2 <- survfit(Surv(time, status)~1, data=survdat[survdat$tr==0&survdat$benefit.survival.at.fixed.time>0,]) ## group 4, see main paper
PB.0=p.pos*(summary(km1, times=t.0)[6][[1]]-summary(km2, times=t.0)[6][[1]])
v12=summary(km1, times=t.0)[7][[1]]^2+summary(km2, times=t.0)[7][[1]]^2
se.PB.0=sqrt(
  v.p.pos*v12+
    v12*(p.pos)^2+
    v.p.pos*(summary(km1, times=t.0)[6][[1]]-summary(km2, times=t.0)[6][[1]])^2
)
lower.PB.0 <- PB.0-1.96*se.PB.0
upper.PB.0 <- PB.0+1.96*se.PB.0
paste(round(PB.0, 3), " [", round(lower.PB.0,3)," ," ,round(upper.PB.0,3),"]", sep="")



# estimate treatment effect in subpopulation of patients with negative predicted benefit
km3 <- survfit(Surv(time, status)~1, data=survdat[survdat$tr==1&survdat$benefit.survival.at.fixed.time<0,]) ## group 3, see main paper
km4 <- survfit(Surv(time, status)~1, data=survdat[survdat$tr==0&survdat$benefit.survival.at.fixed.time<0,]) ## group 4, see main paper
PB.2 <- summary(km3, times=t.0)[6][[1]]-summary(km4, times=t.0)[6][[1]]
SE.PB.2 <- sqrt(summary(km3, times=t.0)[7][[1]]^2+summary(km4, times=t.0)[7][[1]]^2)
lower.PB.2 <- PB.2-1.96*SE.PB.2
upper.PB.2 <- PB.2+1.96*SE.PB.2
paste(round(PB.2, 3), " [", round(lower.PB.2,3)," ," ,round(upper.PB.2,3),"]", sep="")



#### benefit for RMST -------
library(survRM2 )
t.star=1

## PB and treatment effect
rmst_treatment <- with(survdat,rmst2(time, status, tr,tau=t.star))
rmst_model <- with(survdat,rmst2(time, status, agreed,tau=t.star))

# PB 1
PB.RMST=rmst_model$RMST.arm1$result[1]-rmst_treatment$RMST.arm1$result[1]
SE.PB.RMST=sqrt(rmst_model$RMST.arm1$result[1,2]^2+rmst_treatment$RMST.arm1$result[1,2]^2)
lower.PB.RMST<- PB.RMST-1.96*SE.PB.RMST
upper.PB.RMST <- PB.RMST+1.96*SE.PB.RMST
paste(round(PB.RMST, 3), " [", round(lower.PB.RMST,3)," ," ,round(upper.PB.RMST,3),"]", sep="")

# estimate treatment effect in subpopulation of patients with negative predicted benefit
with(survdat[survdat$benefit.survival.at.fixed.time<0,],rmst2(time, status, tr,tau=t.star))



## population benefit for median survival -----
mtr=survfit(Surv(time, status) ~ tr, data = survdat, type = "kaplan-meier")
m.model=survfit(Surv(time, status) ~ agreed, data = survdat, type = "kaplan-meier")

#PB 1
m14=survfit(Surv(time, status) ~ 1, data = survdat[survdat$agreed==1,], type = "kaplan-meier") # Group 1U4
m13=survfit(Surv(time, status) ~ 1, data = survdat[survdat$tr==1,], type = "kaplan-meier") # Group 1U3
summary(m14)$table[7]-summary(m13)$table[7]
   
#Treatment effect in patients with negative benefit 
m34=survfit(Surv(time, status) ~ tr, data = survdat[survdat$benefit.survival.at.fixed.time<0,], type = "kaplan-meier") # Group 1U4


      
