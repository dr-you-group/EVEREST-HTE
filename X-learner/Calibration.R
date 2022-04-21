
##Form a calibration plot of predicted risk reduction, return the
##slope, intercept, predicted and observed ARRs.

test$Time<-ifelse(test$DAYS.TIMI.MAJOR.BLEEDING<360, test$DAYS.TIMI.MAJOR.BLEEDING, 360)

Calibration <- function(data, n.tile){
  sort(data$predict_x_learner_fit)
  test.ex <- mutate(data, quantile_rank = ntile(data$predict_x_learner_fit, n.tile))
  
  table(test.ex$quantile_rank)
  
    pred_rr1 = data.frame()  
    pred_rr2 = data.frame()  
    pred_rr3 = data.frame() 
    obs_rr1 = data.frame()
    obs_rr2 = data.frame()
    
    for (i in (1:n.tile)) {
  
  
      dat = test.ex[test.ex$quantile_rank==i,]
      table(test$Y)
      
      #Predicted ARR for X-learner
      #cal1
      pred_rr1[i,1] = mean(dat[dat$W2==0,]$predict_x_learner_fit)-mean(dat[dat$W2==1,]$predict_x_learner_fit)
      
      #cal2
      pred_rr2[i,1] = ifelse(is.na(mean(dat[dat$predict_x_learner_fit<=0 & dat$W2==0,"Y"]) - mean(dat[dat$predict_x_learner_fit<=0 & dat$W2==1,"Y"]))==FALSE, 
                                   mean(dat[dat$predict_x_learner_fit<=0 & dat$W2==0,"Y"]) - mean(dat[dat$predict_x_learner_fit<=0 & dat$W2==1,"Y"]),
                                   mean(dat[dat$predict_x_learner_fit>0 & dat$W2==0, "Y"]) - mean(dat[dat$predict_x_learner_fit>0 & dat$W2==1, "Y"]))

      pred_rr3[i,1]=pred_rr2[i,1]-pred_rr1[i,1]
      
      
      #Observed ARR for Kaplan-Meier Estimation
      surv.val = survfit(Surv(Time, as.factor(Y))~W2, data = dat)
      inx_w0 = surv.val$strata[[1]];inx_w1=surv.val$strata[[2]]
      
      obs_rr1[i,1] = -(surv.val$pstate[inx_w0+inx_w1,2]-surv.val$pstate[inx_w0,2])
      
      #Observed ARR for Propensity score
      #fmla_ct <- paste("as.factor(Y) ~", paste(covariate_names, collapse = " + "))
      m_ps = glm(as.factor(Y)~W2, family=binomial(), data=dat)
      
      prs_df = data.frame(pr_score = predict(m_ps, type = 'response'), W2=m_ps$model$W2)
      head(prs_df)
      
      obs_rr2[i,1] = mean(prs_df[prs_df$W2==0,"pr_score"])-mean(prs_df[prs_df$W2==1,"pr_score"])
      
    }
  
  paired_rr = cbind(pred_rr1, pred_rr2, pred_rr3, obs_rr1, obs_rr2)
  names(paired_rr) = c("pred_rr1", "pred_rr2", "pred_rr3","obs_rr1", "obs_rr2")
  
  return(paired_rr)
  
  
}

cb_x <- Calibration(test, 5)

model = lm(cb[,4]~cb[,3], data=cb)
slope<-model$coefficients[[2]]#ideally 1
intercept<-model$coefficients[[1]]#ideally 0
summary(model)


plot(cb[,3],cb[,4], pch=19, col="darkgreen", xlab="Predicted ARR", ylab="Observed ARR", xlim=c(0, 0.05), ylim=c(0, 0.05), main="Calibration plot")
abline(a=0,b=1,col="gray",lty=6)
abline(a=intercept, b=slope, col="darkgreen", lty=5, lwd=1.8)
text(x=0.015, y= 0.046, labels = c("Slope: 0.99, Intercept: -0.00"))


##C-statistics
##Return concordance for benefit, the proportion of all matched pairs with 
##unequal observed benefit, in which the patient pair receiving greater
##treatment benefit was predicted to do so.
##

#ensure results are reproducible
seed(123)

