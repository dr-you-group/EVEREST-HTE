# herral's c-index for x-learner model
c_indx<-function(data, gp){
  
  Surv.obj<-Surv(data[,c('Time')], data[,c('Y')])
  cox.obj<-coxph(Surv.obj ~data[,c(gp)], data=data)
  summary(cox.obj)
  tidy(cox.obj)
  
  data.frame(tidy(cox.obj, exponentiate = T, conf.int = T)[, c("term", "estimate", "conf.low", "conf.high", "p.value")])
  conco_ben=concordance(cox.obj)$concordance
  
  
  test.benefit<-subset(data, get(gp)=="benefit")
  Surv.obj2<-Surv(test.benefit$Time, test.benefit$Y)
  cox.obj2<-coxph(Surv.obj2 ~ DAPT , data=test.benefit)
  summary(cox.obj2);summary(cox.obj2)
  tidy(cox.obj2)
  
  data.frame(tidy(cox.obj2, exponentiate = T, conf.int = T)[, c("term", "estimate", "conf.low", "conf.high", "p.value")])
  conco_dapt=concordance(cox.obj2)$concordance
  
  return(list(conco_ben=conco_ben, conco_dapt=conco_dapt))
  
}




# calibration plot(for ARR with x-learner model)
Calib_arr_x <- function(data, n.tile){
  
  test.ex <- mutate(data, quantile_rank = ntile(data$predict_x_learner_fit, n.tile))
  
  pred_rr1 = data.frame()  
  pred_rr2 = data.frame()  
  obs_rr1 = data.frame()
  obs_rr2 = data.frame()
  
  for (i in (1:n.tile)) {
    
    
    dat = test.ex[test.ex$quantile_rank==i,]
    summary(dat$predict_x_learner_fit)
    
    #Predicted ARR for X-learner
    #cal1
    pred_rr1[i,1] = mean(dat$predict_x_learner_fit)#predict_x_learner_fit
    
    #cal2
    pred_rr2[i,1] = ifelse(is.na(mean(dat[dat$predict_x_learner_fit<=0 & dat$W2==0,"predict_x_learner_fit"]) - 
                                   mean(dat[dat$predict_x_learner_fit<=0 & dat$W2==1,"predict_x_learner_fit"]))==FALSE, 
                           mean(dat[dat$predict_x_learner_fit<=0 & dat$W2==0,"predict_x_learner_fit"]) + 
                             mean(dat[dat$predict_x_learner_fit<=0 & dat$W2==1,"predict_x_learner_fit"]),
                           mean(dat[dat$predict_x_learner_fit>0 & dat$W2==0, "predict_x_learner_fit"]) + 
                             mean(dat[dat$predict_x_learner_fit>0 & dat$W2==1, "predict_x_learner_fit"]))
    
    
    #Observed ARR for Kaplan-Meier Estimation
    surv.val.w0 = survfit(Surv(Time, as.factor(Y))~1, data = subset(dat, W2==0))
    surv.val.w1 = survfit(Surv(Time, as.factor(Y))~1, data = subset(dat, W2==1))
    
    surv.val.000 = survfit(Surv(Time, as.factor(Y))~ W2, data = dat)
    summary(surv.val.000)
    
    inx_w0 = nrow(surv.val.w0$cumhaz);inx_w1 = nrow(surv.val.w1$cumhaz);min_n<-min(inx_w0, inx_w1)
    obs_rr1[i,1] = -((surv.val.w0$cumhaz[inx_w0-1,]-surv.val.w1$cumhaz[inx_w1-1,]))
    
    
    #Observed ARR for Logistic model
    m_ps = glm(as.factor(Y)~W2, family=binomial(), data=dat)
    
    prs_df = data.frame(y_hat = predict(m_ps, type = 'response'), W2=m_ps$model$W2)
    
    obs_rr2[i,1] = (mean(prs_df[prs_df$W2==0,"y_hat"])-mean(prs_df[prs_df$W2==1,"y_hat"]))
    
    i+1
  }
  
  paired_rr = cbind(pred_rr1, pred_rr2, obs_rr1, obs_rr2)
  names(paired_rr) = c("pred_rr1", "pred_rr2", "obs_rr1", "obs_rr2")
  
  return(paired_rr)
  
  
}

# calibration plot(for Risk-Probabbility with x-learner model)
Calib_risk_x <- function(data, n.tile){
  
  m_ps = glm(as.factor(test$Y)~predict_x_learner_fit2+W2, family=binomial(), data=data)
  prs_df = data.frame(pr_score = predict(m_ps, type = 'response'))
  
  test2<-cbind(data, prs_df)
  
  test.ex <- mutate(test2, quantile_rank = ntile(test2$pr_score, n.tile))
  
  pred_prr = data.frame()  
  obs_prr = data.frame()  

  
  for (i in (1:n.tile)) {
    
    pdat = test.ex[test.ex$quantile_rank==i,]
    
    pred_prr[i,1] =mean(pdat$pr_score)
    
    #Observed ARR for Propensity score
    
    obs_prr[i,1] = mean(pdat$Y)
    
  }
  paired_rr = cbind(pred_prr, obs_prr)
  names(paired_rr) = c("pred_prr", "obs_prr")
  
  return(paired_rr)
}

# calibration plot(for Risk-Probabbility with PRECISE-DAPT model)
Calib_risk_dapt <- function(data, n.tile){
  
  m_ps = glm(as.factor(dapt.test$Y)~precise_DAPT+W2, family=binomial(), data=dapt.test)
  prs_df = data.frame(pr_score = predict(m_ps, type = 'response'))
  
  dapt.test2<-cbind(dapt.test, prs_df)
  
  test.ex <- mutate(dapt.test2, quantile_rank = ntile(dapt.test2$pr_score, n.tile))
  
  pred_prr = data.frame()  
  obs_prr = data.frame()  
  
  for (i in (1:n.tile)) {
    
    pdat = test.ex[test.ex$quantile_rank==i,]
    
    pred_prr[i,1] =mean(pdat$pr_score)
    
    #Observed ARR for Propensity score
    m_ps2 = glm(as.factor(Y)~as.factor(W2), family=binomial(), data=pdat)
    summary(m_ps2)
    
    prs_df2 = data.frame(y_hat = predict(m_ps2, type = 'response'), W2=pdat$W2)
    
    #obs_prr[i,1] = mean(prs_df2$y_hat)
    obs_prr[i,1] = mean(pdat$Y)
  }
  paired_rr = cbind(pred_prr, obs_prr)
  names(paired_rr) = c("pred_prr", "obs_prr")
  
  return(paired_rr)
}




# calibration plot(for ARR with PRECISE-DAPT model)
Calib_ARR_dapt <- function(data, n.tile){
  
  test.ex <- mutate(data, quantile_rank = ntile(data$predict_x_learner_fit, n.tile))
  
  pred_rr1 = data.frame()  
  pred_rr2 = data.frame()  
  obs_rr1 = data.frame()
  obs_rr2 = data.frame()
  
  
  for (i in (1:n.tile)) {
    
    
    dat = test.ex[test.ex$quantile_rank==i,]
    summary(dat$predict_x_learner_fit)
    #Predicted ARR for X-learner
    #cal1
    pred_rr1[i,1] = mean(dat$predict_x_learner_fit)#predict_x_learner_fit
    
    #cal2
    pred_rr2[i,1] = ifelse(is.na(mean(dat[dat$predict_x_learner_fit<=0 & dat$W2==0,"predict_x_learner_fit"]) - mean(dat[dat$predict_x_learner_fit<=0 & dat$W2==1,"predict_x_learner_fit"]))==FALSE, 
                           mean(dat[dat$predict_x_learner_fit<=0 & dat$W2==0,"predict_x_learner_fit"]) + mean(dat[dat$predict_x_learner_fit<=0 & dat$W2==1,"predict_x_learner_fit"]),
                           mean(dat[dat$predict_x_learner_fit>0 & dat$W2==0, "predict_x_learner_fit"]) + mean(dat[dat$predict_x_learner_fit>0 & dat$W2==1, "predict_x_learner_fit"]))
    
    
    #Observed ARR for Kaplan-Meier Estimation

    if (sum(dat[dat$W2==0,"Y"])==0 & sum(dat[dat$W2==1,"Y"])==0) {
      obs_rr1[i,1]= 0.0001
      
      
    } else if (sum(dat[dat$W2==0,"Y"])!=0 & sum(dat[dat$W2==1,"Y"])!=0) {
      surv.val.w0 = survfit(Surv(Time, as.factor(Y))~1, data = subset(dat, W2==0))
      surv.val.w1 = survfit(Surv(Time, as.factor(Y))~1, data = subset(dat, W2==1))
      
      inx_w0 = nrow(surv.val.w0$cumhaz);inx_w1 = nrow(surv.val.w1$cumhaz);min_n<-min(inx_w0, inx_w1)
      
      obs_rr1[i,1] = -((surv.val.w0$cumhaz[inx_w0-1,]-surv.val.w1$cumhaz[inx_w1-1,]))
      
    } else if (sum(dat[dat$W2==0,"Y"])==0 & sum(dat[dat$W2==1,"Y"])!=0){
      surv.val.w1 = survfit(Surv(Time, as.factor(Y))~1, data = subset(dat, W2==1))
      
      inx_w1 = nrow(surv.val.w1$cumhaz)
      
      obs_rr1[i,1] = -((0-surv.val.w1$cumhaz[inx_w1-1,]))
      
    } else if (sum(dat[dat$W2==0,"Y"])!=0 & sum(dat[dat$W2==1,"Y"])==0){ 
      surv.val.w0 = survfit(Surv(Time, as.factor(Y))~1, data = subset(dat, W2==0))
      
      inx_w0 = nrow(surv.val.w0$cumhaz)
      
      obs_rr1[i,1] = -((surv.val.w0$cumhaz[inx_w0-1,]-0))
    }
    

    #Observed ARR for logistic
    m_ps = glm(as.factor(Y)~W2, family=binomial(), data=dat)
    
    prs_df = data.frame(y_hat = predict(m_ps, type = 'response'), W2=m_ps$model$W2)
    
    obs_rr2[i,1] = (mean(prs_df[prs_df$W2==0,"y_hat"])-mean(prs_df[prs_df$W2==1,"y_hat"]))
    
    i+1
  }
  
  paired_rr = cbind(pred_rr1, pred_rr2, obs_rr1, obs_rr2)
  names(paired_rr) = c("pred_rr1", "pred_rr2", "obs_rr1", "obs_rr2")
  
  return(paired_rr)
  
  
}
