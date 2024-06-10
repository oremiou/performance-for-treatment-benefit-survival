## fit prediction model
library(survival)
fit.model <- survreg(Surv(time, status)~x1+x2+x3+x4+tr+x1*tr+x2*tr+x3*tr+x4*tr, 
                     survdat, dist="weibul")

#  predict survival curve for new patient ---- 
new.patient0 <- survdat[1,]; new.patient0$tr="0" 
new.patient1 <- survdat[1,]; new.patient1$tr="1" 

pct <- 1:99/100  
ptime0 <- predict(fit.model, newdata=new.patient0, type='quantile', p=pct, se=TRUE)
ptime1 <- predict(fit.model, newdata=new.patient1, type='quantile', p=pct, se=TRUE)
predicted <- rbind(data.frame(Z="0","percentiles.surv"=pct,"pred"=ptime0$fit), 
                   data.frame(Z="1", "percentiles.surv"=pct,"pred"=ptime1$fit))

library(ggplot2)
colnames(predicted)[1]="Intervention_group"
ggplot(predicted, aes(x = pred, y =1 - percentiles.surv , color=Intervention_group)) + 
  geom_line( linewidth=.6)+ xlab("Years of follow-up")+
  ylab("Predicted survival probability")+
  scale_color_manual(values = c("black", "red")) +
  xlim(0,3) 



# for patient 2
#  predict survival curve for new patient ---- 
new.patient0 <- survdat[2,]; new.patient0$tr="0" 
new.patient1 <- survdat[2,]; new.patient1$tr="1" 

pct <- 1:99/100  
ptime0 <- predict(fit.model, newdata=new.patient0, type='quantile', p=pct, se=TRUE)
ptime1 <- predict(fit.model, newdata=new.patient1, type='quantile', p=pct, se=TRUE)
predicted <- rbind(data.frame(Z="0","percentiles.surv"=pct,"pred"=ptime0$fit), 
                   data.frame(Z="1", "percentiles.surv"=pct,"pred"=ptime1$fit))

library(ggplot2)
colnames(predicted)[1]="Intervention_group"
ggplot(predicted, aes(x = pred, y =1 - percentiles.surv , color=Intervention_group)) + 
  geom_line( linewidth=.6)+ xlab("Years of follow-up")+
  ylab("Predicted survival probability")+
  scale_color_manual(values = c("black", "red")) +
  xlim(0,3) 


# benefit in survival probability
prediction.weibull <-
  function(new.patient, model = NULL){
    pct <- 1:300/100
    mu_hat <- predict(model, newdata = new.patient, type = "lp")
    predicted.exp <- 1 - pweibull(pct, shape = 1/model$scale, 
                                  scale = exp(mu_hat))
    
    results <- data.frame(time = pct, means = predicted.exp)
    return(results)
  }  

prediction.weibull(new.patient0, model = fit.model)[200,2]-
prediction.weibull(new.patient1, model = fit.model)[200,2] # benefit at 2 years

# predicted benefit expected survival time
pct=1:999/1000
ptime.mean0 <- predict(fit.model, new.patient0, 
                       type='quantile', p=pct, se=TRUE)$fit
sum(ptime.mean0 *0.001)

ptime.mean1 <- predict(fit.model, new.patient1, 
                       type='quantile', p=pct, se=TRUE)$fit

sum(ptime.mean0 *0.001)-sum(ptime.mean1 *0.001)

# predicted benefit restricted mean survival time
ptime.mean1[ptime.mean1>2]=2
ptime.mean0[ptime.mean0>2]=2
sum(ptime.mean0 *0.001)-sum(ptime.mean1 *0.001)


# predicted benefit median survival time
predict(fit.model, new.patient0, 
        type='quantile', p=0.5, se=TRUE)$fit-
  predict(fit.model, new.patient1, 
          type='quantile', p=0.5, se=TRUE)$fit






