
# Evaluation of model bucketing for "Benefit", "no-benefit" according to ARR 
# Benefit: predicted risk<=0; no-benefit: predicted risk>0

##W=1 then treat is short dapt(<6mo)




Xbucket_arr <- function(df_x_pred, predict_x_learner_fit, Y, W) {
  
  Xarr_ben = mean(df_x_pred[predict_x_learner_fit<=0 & df_x_pred$W==0,"Y"]) - mean(df_x_pred[predict_x_learner_fit<=0 & df_x_pred$W==1,"Y"])
  Xarr_noben = mean(df_x_pred[predict_x_learner_fit>0 & df_x_pred$W==0, "Y"]) - mean(df_x_pred[predict_x_learner_fit>0 & df_x_pred$W==1, "Y"])
  
  # emprical 95% range of a bucket of benefit/no-benefit
  #Xarr_benL = -quantile(df_x_pred[df_x_pred$predict_x_learner_fi<=0,"predict_x_learner_fit"], 0.025)
  #Xarr_benU = -quantile(df_x_pred[df_x_pred$predict_x_learner_fi<=0,"predict_x_learner_fit"], 0.975)
  
  # 95% CI for normal dist
  #Xarr_ben_sd = stdev(df_x_pred[df_x_pred$predict_x_learner_fit<=0,"predict_x_learner_fit"])
  #Xarr_benL = Xarr_ben-Xarr_ben_sd*1.96
  #Xarr_benU = Xarr_ben+Xarr_ben_sd*1.96
  #Xarr_benL;Xarr_ben;Xarr_benU
  
  #Xarr_noben_sd = stdev(df_x_pred[df_x_pred$predict_x_learner_fit>0,"predict_x_learner_fit"])
  #Xarr_nobenL = Xarr_noben-Xarr_noben_sd*1.96
  #Xarr_nobenU = Xarr_noben+Xarr_noben_sd*1.96
  #Xarr_nobenL;Xarr_noben;Xarr_nobenU
  
  
  # 95% CI for binomial dist
  Xarr_ben_sd = sqrt(Xarr_ben*(1-Xarr_ben)/nrow(df_x_pred[predict_x_learner_fit<=0,]))
  Xarr_benL = Xarr_ben-1.96*Xarr_ben_sd
  Xarr_benU = Xarr_ben+1.96*Xarr_ben_sd
  Xarr_benL;Xarr_ben;Xarr_benU
  
  Xarr_noben_sd = sqrt(Xarr_noben*(1-Xarr_noben)/nrow(df_x_pred[predict_x_learner_fit>0,]))
  Xarr_nobenL = -(Xarr_noben-1.96*Xarr_noben_sd)
  Xarr_nobenU = -(Xarr_noben+1.96*Xarr_noben_sd)
  Xarr_nobenL;-Xarr_noben;Xarr_nobenU
  
  # Wald test to return p-value for the event rate in each of the buckets
  for (i in 1:2) {
    
    if(i==1) { # Benefit
      
      #trt = df_x_pred[predict_x_learner_fit<=0 & df_x_pred$W==1,"Y"]
      #ctr = df_x_pred[predict_x_learner_fit<=0 & df_x_pred$W==0, "Y"]
      
      benefit.dat = df_x_pred[predict_x_learner_fit<=0,]
      
      #trt_y1 = sum(trt);trt_y0 = length(trt) - trt_y1;
      #ctr_y1 = sum(ctr);ctr_y0 = length(ctr) - ctr_y1;
      
      #benefit = matrix(c(trt_y1, trt_y0, ctr_y1, ctr_y0), nrow=2, ncol=2)  
      
      benefit_chisq <- chisq.test(df_x_pred[predict_x_learner_fit<=0,]$Y, df_x_pred[predict_x_learner_fit<=0,]$W)
      benefit_p.val = benefit_chisq$p.value
      
    } else { #No-benefit
      
      #trt <- df_x_pred[predict_x_learner_fit>0 & df_x_pred$W==1,"Y"] 
      #ctr <- df_x_pred[predict_x_learner_fit>0 & df_x_pred$W==0,"Y"]
      
      nobenefit.dat = df_x_pred[predict_x_learner_fit>0,]
      
      #trt_y1 = sum(trt);trt_y0 = length(trt) - trt_y1;
      #ctr_y1 = sum(ctr);ctr_y0 = length(ctr) - ctr_y1;
      
      #nobenefit = matrix(c(trt_y1, trt_y0, ctr_y1, ctr_y0), nrow=2, ncol=2)  
      
      nobenefit_chisq = chisq.test(df_x_pred[predict_x_learner_fit>0,]$Y, df_x_pred[predict_x_learner_fit>0,]$W)
      nobenefit_p.val = nobenefit_chisq$p.value
      
    }
  }
  
  
  return(list(Xarr_ben=Xarr_ben, Xarr_noben=Xarr_noben, benefit=benefit, nobenefit=nobenefit))
  
}



#Forestplot for Average Risk Reduction for Each group(Benefit/No-benefit)
library(ggplot2)

mytheme <- function(base_size = 12, base_family = "sans"){
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      panel.grid.major = element_line(color = "grey"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "aliceblue"),
      strip.background = element_rect(fill = "darkgrey", color = "grey", size = 1),
      strip.text = element_text(face = "bold", size = 12, color = "white"),
      legend.position = " ",
      legend.justification = "top", 
      panel.border = element_rect(color = "grey", fill = NA, size = 0.5)
    )
}

Parameter  = c("Benefit", "Nobenefit")
Model = c(rep("X-learn", 2))

RR = c(Xarr_ben, -Xarr_noben)
Lower = c(Xarr_benL, Xarr_nobenL)
Upper = c(Xarr_benU, Xarr_nobenU)


dat = data.frame(Parameter, RR, Lower, Upper)

dat$Parameter = factor(dat$Parameter, levels=c("Benefit", "Nobenefit")) 


g = ggplot(data=dat, aes(x=Model, y=RR, ymin=Lower, ymax=Upper, color=Parameter))#, color=Parameter
g = g + geom_pointrange(aes(col=c("red", "blue")), lwd=0.5)  
g = g + geom_hline(aes(fill=Parameter), yintercept = 0, linetype=2)
g = g + xlab("") + ylab("Average Risk Reduction (95% Confidence Interval)")
g = g + geom_errorbar(aes(ymin=Lower, ymax=Upper, col=Model), width=0.5, cex=1)
g = g + facet_wrap(~Parameter, strip.position="left", nrow=9, scales = "free_y")
g = g + theme(plot.title=element_text(size=16,face="bold"),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.x=element_text(face="bold"),
              axis.title=element_text(size=12,face="bold"),
              strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))

g + coord_flip() + theme_bw() + theme(legend.position="none")




  

  