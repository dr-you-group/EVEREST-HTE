precise.dapt <- read.csv(file.path(dataFolder,"precise_dapt.csv")) 

# merge
dapt.test<-inner_join(precise.dapt, df_test, by='patient_id')

nrow(dapt.test)

# survival plot for DAPT score
dapt.test$DAPT<-ifelse(dapt.test$W2==1,"Short", "Long")
dapt.test$Time<-ifelse(dapt.test$DAYS.TIMI.MAJOR.BLEEDING<360, 
                       dapt.test$DAYS.TIMI.MAJOR.BLEEDING, 360)
dapt.test$class<-ifelse(dapt.test$precise_DAPT<25, "benefit", "no-benefit")



surv_func = function(test){
  
  
  fit1<-survfit(Surv(Time, Y) ~ class, data = test)
  fit2<-survfit(Surv(Time, Y) ~ DAPT, data = subset(test, class=="benefit"))
  fit3<-survfit(Surv(Time, Y) ~ DAPT, data = subset(test, class=="no-benefit"))
  
  
  
  summary(fit1);summary(fit2);summary(fit3)
  
  return(list(test=test, fit1=fit1, fit2=fit2, fit3=fit3, test=test))
  
}      

surv_d<-surv_func(dapt.test)

survplot = c()

form<-c()
dat<-c()

form[[1]]<-formula("Surv(Time, Y) ~ class")
dat[[1]]<-surv_d$test

form[[2]]<-formula("Surv(Time, Y) ~ DAPT")
dat[[2]]<-subset(surv_d$test, class=="benefit")

form[[3]]<-formula("Surv(Time, Y) ~ DAPT")
dat[[3]]<-subset(surv_d$test, class=="no-benefit")


for (i in 1:3) {
  
  title.nm<-c("TIMI.MAJOR.BLEEDING", 
              "TIMI.MAJOR.BLEEDING for Benefit", 
              "TIMI.MAJOR.BLEEDING for No-benefit")
  
  y.max<-c(0.2, 0.1, 0.2)
  
  legend.lab<-c()
  legend.lab[[1]]<-c("Benefit(Short)", "No-benefit")
  legend.lab[[2]]<-c("DAPT 12 month", "DAPT 6 month")
  legend.lab[[3]]<-c("DAPT 12 month", "DAPT 6 month")
  
  
  palette<-c()
  palette[[1]]<-c("#7CAE00", "#FF6600")
  palette[[2]]<-"circulation"
  palette[[3]]<-"circulation"
  
  #plotting
  survplot[[i]]<-ggsurvplot(survfit(form[[i]], data =dat[[i]]),
                            conf.int = T,
                            pval = F,
                            censor=F,
                            break.time.by=90, 
                            break.y.by=0.03,
                            pval.method=TRUE,
                            risk.table.col = "strata", # Change risk table color by groups
                            ggtheme = theme_bw(), # Change ggplot2 theme
                            palette = palette[[i]],
                            fun = "cumhaz",
                            xlim = c(0, 360), ylim = c(0, y.max[i]),
                            risk.table = TRUE,
                            title = title.nm[i],#Cardiac death, MI, stent thrombosis, CVA
                            fontsize=4.0,     
                            risk.table.pos="out",
                            risk.table.y.text.col=FALSE,
                            #tables.theme = theme_cleantable(),
                            #ggtheme=theme_classic2(base_size=12, base_family="Arial"),
                            #font.family="Arial",
                            legend.labs = legend.lab[[i]],
                            legend.title=""
                            
  ) 
  ggpar(survplot[[i]],
        font.main = c(13, "bold.italic", "darkblue"),
        font.x = c(11, "bold", "black"),
        font.y = c(11, "bold", "black"),
        font.tickslab = c(11, "plain", "black"), 
        font.legend = c(11, "bold") 
  )
}

survplot


# Log-Lank Test
km_diff<-c()
km_pval<-c()
tmp<-c()
cox_out<-c()
cox<-c()

dat[[1]]$class <- relevel(as.factor(dat[[1]]$class), ref = "no-benefit")
for (i in 1 : 3){
  
  km_diff[[i]] <-survdiff(form[[i]], data =dat[[i]])
  km_pval[i] <- 1 - pchisq(km_diff[[i]]$chisq, length(km_diff[[i]]$n) - 1)
  
  
  cox_out[[i]] <- coxph(form[[i]], data =dat[[i]]);
  cox[[i]]<-as.matrix(summary(eval(cox_out[[i]]))$conf.int)
  
  options(digits=3)
  
  print(cox_out[[i]])
  
  i=i+1
  
}

km_pval_all<-data.frame(rbind(km_pval[[1]], km_pval[[2]], km_pval[[3]]))
names(km_pval_all)<-c("bind.km_pval..1....km_pval..2....km_pval..3..."="log-lank")
cox_all<-data.frame(rbind(cox[[1]], cox[[2]], cox[[3]]))
Result<-cbind(km_pval_all, cox_all)
Result


# Apply Model for other outcome(NACE, MACE)
source("C:/git/EVEREST-HTE/Final/Apply.R")

result_ex<-aply("AD.MI.ST.CVA.MAJOR.BLEEDING")

form<-c();dat<-c()

form[[1]]<-formula("Surv(Time, Y) ~ class_ex")
dat[[1]]<-result_ex$test_ex

form[[2]]<-formula("Surv(Time, Y) ~ DAPT")
dat[[2]]<-subset(result_ex$test_ex, class_ex=="benefit")

form[[3]]<-formula("Surv(Time, Y) ~ DAPT")
dat[[3]]<-subset(result_ex$test_ex, class_ex=="no-benefit")


for (i in 1:3) {
  
  title.nm<-c("AD.MI.ST.CVA.MAJOR.BLEEDING", 
              "AD.MI.ST.CVA.MAJOR.BLEEDING for Benefit", 
              "AD.MI.ST.CVA.MAJOR.BLEEDING for No-benefit")
  
  y.max<-c(0.2, 0.1, 0.2)
  
  legend.lab<-c()
  legend.lab[[1]]<-c("Benefit(Short)", "No-benefit")
  legend.lab[[2]]<-c("DAPT 12 month", "DAPT 6 month")
  legend.lab[[3]]<-c("DAPT 12 month", "DAPT 6 month")
  
  
  palette<-c()
  palette[[1]]<-c("#7CAE00", "#FF6600")
  palette[[2]]<-"circulation"
  palette[[3]]<-"circulation"
  
  # plotting
  survplot[[i]]<-ggsurvplot(survfit(form[[i]], data =dat[[i]]),
                            conf.int = T,
                            pval = F,
                            censor=F,
                            break.time.by=90, 
                            break.y.by=0.05,
                            pval.method=TRUE,
                            risk.table.col = "strata", # Change risk table color by groups
                            ggtheme = theme_bw(), # Change ggplot2 theme
                            palette = palette[[i]],
                            fun = "cumhaz",
                            xlim = c(0, 360), ylim = c(0, y.max[i]),
                            risk.table = TRUE,
                            title = title.nm[i],#Cardiac death, MI, stent thrombosis, CVA
                            fontsize=4.0,     
                            risk.table.pos="out",
                            risk.table.y.text.col=FALSE,
                            #tables.theme = theme_cleantable(),
                            #ggtheme=theme_classic2(base_size=12, base_family="Arial"),
                            #font.family="Arial",
                            legend.labs = legend.lab[[i]],
                            legend.title=""
                            
  ) 
  ggpar(survplot[[i]],
        font.main = c(13, "bold.italic", "darkblue"),
        font.x = c(11, "bold", "black"),
        font.y = c(11, "bold", "black"),
        font.tickslab = c(11, "plain", "black"), 
        font.legend = c(11, "bold")
  )
}

survplot

# Log-Lank Test
km_diff<-c();km_pval<-c()
cox_out<-c();cox<-c()

dat[[1]]$class_ex <- relevel(as.factor(dat[[1]]$class_ex), ref = "no-benefit")

for (i in 1 : 3){
  
  km_diff[[i]] <-survdiff(form[[i]], data =dat[[i]])
  km_pval[i] <- 1 - pchisq(km_diff[[i]]$chisq, length(km_diff[[i]]$n) - 1)
  
  cox_out[[i]] <- coxph(form[[i]], data =dat[[i]]);
  cox[[i]]<-as.matrix(summary(eval(cox_out[[i]]))$conf.int)
  
  options(digits=3)
  
  print(cox_out[[i]])
  
  i=i+1
  
}

km_pval_all<-data.frame(rbind(km_pval[[1]], km_pval[[2]], km_pval[[3]]))
names(km_pval_all)<-c("bind.km_pval..1....km_pval..2....km_pval..3..."="log-lank")
cox_all<-data.frame(rbind(cox[[1]], cox[[2]], cox[[3]]))
Result<-cbind(km_pval_all, cox_all)
Result


result_ex<-aply("CD.MI.ST.CVA")

form<-c()
dat<-c()

form[[1]]<-formula("Surv(Time, Y) ~ class_ex")
dat[[1]]<-result_ex$test_ex
form[[2]]<-formula("Surv(Time, Y) ~ DAPT")
dat[[2]]<-subset(result_ex$test_ex, class_ex=="benefit")
form[[3]]<-formula("Surv(Time, Y) ~ DAPT")
dat[[3]]<-subset(result_ex$test_ex, class_ex=="no-benefit")


for (i in 1:3) {
  
  title.nm<-c("CD.MI.ST.CVA", 
              "CD.MI.ST.CVA for Benefit", 
              "CD.MI.ST.CVA for No-benefit")
  
  y.max<-c(0.1, 0.1, 0.1)
  
  legend.lab<-c()
  legend.lab[[1]]<-c("Benefit(Short)", "No-benefit")
  legend.lab[[2]]<-c("DAPT 12 month", "DAPT 6 month")
  legend.lab[[3]]<-c("DAPT 12 month", "DAPT 6 month")
  
  
  palette<-c()
  palette[[1]]<-c("#7CAE00", "#FF6600")
  palette[[2]]<-"circulation"
  palette[[3]]<-"circulation"
  
  # plotting
  survplot[[i]]<-ggsurvplot(survfit(form[[i]], data =dat[[i]]),
                            conf.int = T,
                            pval = T,
                            censor=F,
                            break.time.by=90, 
                            break.y.by=0.01,
                            pval.method=TRUE,
                            risk.table.col = "strata", # Change risk table color by groups
                            ggtheme = theme_bw(), # Change ggplot2 theme
                            palette = palette[[i]],
                            fun = "cumhaz",
                            xlim = c(0, 360), ylim = c(0, y.max[i]),
                            risk.table = TRUE,
                            title = title.nm[i],#Cardiac death, MI, stent thrombosis, CVA
                            fontsize=4.0,     
                            risk.table.pos="out",
                            risk.table.y.text.col=FALSE,
                            #tables.theme = theme_cleantable(),
                            #ggtheme=theme_classic2(base_size=12, base_family="Arial"),
                            #font.family="Arial",
                            legend.labs = legend.lab[[i]],
                            legend.title=""
                            
  ) 
  ggpar(survplot[[i]],
        font.main = c(13, "bold.italic", "darkblue"),
        font.x = c(11, "bold", "black"),
        font.y = c(11, "bold", "black"),
        font.tickslab = c(11, "plain", "black"), 
        font.legend = c(11, "bold")
  )
}

survplot

# Log-Lank Test
km_diff<-c();km_pval<-c()
cox_out<-c();cox<-c()

dat[[1]]$class_ex <- relevel(as.factor(dat[[1]]$class_ex), ref = "no-benefit")

for (i in 1 : 3){
  
  km_diff[[i]] <-survdiff(form[[i]], data =dat[[i]])
  km_pval[i] <- 1 - pchisq(km_diff[[i]]$chisq, length(km_diff[[i]]$n) - 1)
  
  cox_out[[i]] <- coxph(form[[i]], data =dat[[i]]);
  cox[[i]]<-as.matrix(summary(eval(cox_out[[i]]))$conf.int)
  
  options(digits=3)
  
  print(cox_out[[i]])
  
  i=i+1
  
}

km_pval_all<-data.frame(rbind(km_pval[[1]], km_pval[[2]], km_pval[[3]]))
names(km_pval_all)<-c("bind.km_pval..1....km_pval..2....km_pval..3..."="log-lank")
cox_all<-data.frame(rbind(cox[[1]], cox[[2]], cox[[3]]))
Result<-cbind(km_pval_all, cox_all)
Result

# herral's c-index for PRECISE-DAPT model
source("C:/git/EVEREST-HTE/Final/Calibration.R")
c_indx(surv_d$test, "class")


# C-for-benefit optimism
mv.cph <- cph(Surv(as.numeric(Time), as.numeric(Y))  ~ 
                DAPT,x=T,y = T, data = subset(dapt.test, class=="benefit")) 
bootcov(mv.cph)

set.seed(2222)
mv.cph.val <- validate(mv.cph, method = "boot", B = 1000, )
0.5*(1+mv.cph.val[1,5]) # Corrected C for benefit

# RMST
rms_ben<-rmst2(dapt.test[dapt.test$class=="benefit",]$Time, dapt.test[dapt.test$class=="benefit",]$Y, dapt.test[dapt.test$class=="benefit",]$W2, tau=360)
rms_noben<-rmst2(dapt.test[dapt.test$class!="benefit",]$Time, dapt.test[dapt.test$class!="benefit",]$Y, dapt.test[dapt.test$class!="benefit",]$W2, tau=360)

rms_ben$RMST.arm1
rms_noben$RMST.arm0

(rms_ben$RMST.arm1$rmtl[[2]]+rms_noben$RMST.arm0$rmtl[[2]])/2

(rms_ben$RMST.arm1$rmtl[[1]] * rms_ben$RMST.arm1$fit[[1]] + rms_noben$RMST.arm0$rmtl[[1]] * rms_noben$RMST.arm0$fit[[1]])/(sum(rms_noben$RMST.arm0$fit$n.event)+sum(rms_ben$RMST.arm1$fit$n.event))
(rms_ben$RMST.arm1$rmst[[1]] + rms_noben$RMST.arm0$rmst[[1]]) /(sum(rms_noben$RMST.arm0$fit$n.event)+sum(rms_ben$RMST.arm1$fit$n.event))



fit1<-lrm(Y ~ DAPT,
          data = subset(dapt.test, class=="benefit"), x = T, y = T)
cal1 <- calibrate(fit1, method = 'boot', B = 1000)
plot(cal1,xlim = c(0,0.1),ylim = c(0,0.1))

# calibration plot(for Risk-Probability with PRECISE-DAPT model)
dapt.test$pred_Y<-ifelse(dapt.test$precise_DAPT<25,0,1)

cb_precise <- Calib_risk_dapt(dapt.test, 5)
cb_precise<-as.data.frame(cb_precise)

model = glm(obs_prr~pred_prr, data=cb_precise)
slope<-model$coefficients[[2]]#ideally 1
intercept<-model$coefficients[[1]]#ideally 0
summary(model)

Predicted.ARR<-cb_precise[,1]
Observed.ARR<-cb_precise[,2]

cb_precise<-data.frame(Predicted.ARR, Observed.ARR)

#cb_precise<-cbind(cb_precise[,1], cb_precise[,2])

plot(cb_precise[,1],cb_precise[,2], 
     pch=19, col="darkblue", xlab="Predicted Risk-Probability", 
     ylab="Observed Risk-Probability", xlim=c(0, 0.1), ylim=c(0, 0.1), 
     main="Precise Dapt Score")

abline(a=0,b=1,col="gray",lty=6)
abline(a=intercept, b=slope, col="darkblue", lty=5, lwd=1.8)
text(x=0.033, y= 0.095, 
     labels = c(paste0("Slope: ", round(slope, digit=2), ",  Intercept: ", 
                                          round(intercept, digit=3))))




plt_risk_dapt<-ggplot(data=cb_precise, aes(x=Predicted.ARR, y=Observed.ARR)) + 
  geom_point( size=2, colour="darkblue") + # shape 15: solid square
  ggtitle("Calibration for PRECISE-DAPT SCORE")+
  theme_test()+
  theme(axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11))+
  geom_abline(intercept= intercept, slope=slope, color='darkblue', size = 0.5, linetype="dashed")+
  geom_abline(intercept= 0, slope=1, color='darkgray', size = 0.5, linetype="solid")



plt_risk_dapt<-plt_risk_dapt+
               coord_cartesian(xlim=c(0, 0.1), ylim=c(0, 0.1))#+
              #annotate("text", x=0.03, y=0.08, label="Versicolor", size=7)




# fitting PRECISE-DAPT score into X-learner model
covariate_names<-c("precise_DAPT")

source("C:/git/EVEREST-HTE/Final/X-learner.R")
set.seed(1111)
train_x_learner_fit<-train_x_learner_dapt(dapt.test[,covariate_names, drop=FALSE], 
                                          dapt.test$W2, dapt.test$Y, ipcw=0.5)
predict_x_learner_fit<-predict_x_learner(dapt.test[,covariate_names, drop=FALSE], 
                                         dapt.test$W2, FALSE, TRUE)

dapt.test.pred<-cbind(dapt.test, predict_x_learner_fit)

surv_d<-surv_func(dapt.test.pred)

survplot = c()

form<-c()
dat<-c()

form[[1]]<-formula("Surv(Time, Y) ~ class")
dat[[1]]<-surv_d$test

form[[2]]<-formula("Surv(Time, Y) ~ DAPT")
dat[[2]]<-subset(surv_d$test, class=="benefit")

form[[3]]<-formula("Surv(Time, Y) ~ DAPT")
dat[[3]]<-subset(surv_d$test, class=="no-benefit")


for (i in 1:3) {
  
  title.nm<-c("TIMI.MAJOR.BLEEDING", 
              "TIMI.MAJOR.BLEEDING for Benefit", 
              "TIMI.MAJOR.BLEEDING for No-benefit")
  
  y.max<-c(0.1, 0.1, 0.2)
  
  legend.lab<-c()
  legend.lab[[1]]<-c("Benefit(Short)", "No-benefit")
  legend.lab[[2]]<-c("DAPT 12 month", "DAPT 6 month")
  legend.lab[[3]]<-c("DAPT 12 month", "DAPT 6 month")
  
  
  palette<-c()
  palette[[1]]<-c("#7CAE00", "#FF6600")
  palette[[2]]<-"circulation"
  palette[[3]]<-"circulation"
  
  #plotting
  survplot[[i]]<-ggsurvplot(survfit(form[[i]], data =dat[[i]]),
                            conf.int = T,
                            pval = F,
                            censor=F,
                            break.time.by=90, 
                            break.y.by=0.01,
                            pval.method=TRUE,
                            risk.table.col = "strata", # Change risk table color by groups
                            ggtheme = theme_bw(), # Change ggplot2 theme
                            palette = palette[[i]],
                            fun = "cumhaz",
                            xlim = c(0, 360), ylim = c(0, y.max[i]),
                            risk.table = TRUE,
                            title = title.nm[i],#Cardiac death, MI, stent thrombosis, CVA
                            fontsize=4.0,     
                            risk.table.pos="out",
                            risk.table.y.text.col=FALSE,
                            #tables.theme = theme_cleantable(),
                            #ggtheme=theme_classic2(base_size=12, base_family="Arial"),
                            font.family="Arial",
                            legend.labs = legend.lab[[i]],
                            legend.title=""
                            
  ) 
  ggpar(survplot[[i]],
        font.main = c(13, "bold.italic", "darkblue"),
        font.x = c(11, "bold", "black"),
        font.y = c(11, "bold", "black"),
        font.tickslab = c(11, "plain", "black"), 
        font.legend = c(11, "bold") 
  )
}

survplot

# Log-Lank Test
km_diff<-c()
km_pval<-c()
tmp<-c()
cox_out<-c()
cox<-c()

dat[[1]]$class <- relevel(as.factor(dat[[1]]$class), ref = "no-benefit")
for (i in 1 : 3){
  
  km_diff[[i]] <-survdiff(form[[i]], data =dat[[i]])
  km_pval[i] <- 1 - pchisq(km_diff[[i]]$chisq, length(km_diff[[i]]$n) - 1)
  
  
  cox_out[[i]] <- coxph(form[[i]], data =dat[[i]]);
  cox[[i]]<-as.matrix(summary(eval(cox_out[[i]]))$conf.int)
  
  options(digits=3)
  
  print(cox_out[[i]])
  
  i=i+1
  
}

km_pval_all<-data.frame(rbind(km_pval[[1]], km_pval[[2]], km_pval[[3]]))
names(km_pval_all)<-c("bind.km_pval..1....km_pval..2....km_pval..3..."="log-lank")
cox_all<-data.frame(rbind(cox[[1]], cox[[2]], cox[[3]]))
Result<-cbind(km_pval_all, cox_all)
Result


# herral's c-index
source("C:/git/EVEREST-HTE/Final/Calibration.R")
c_indx(surv_d$test, "class")

# calibration plot(for ARR with PRECISE-DAPT model)
dapt.test.pred$Time<-ifelse(dapt.test.pred$DAYS.TIMI.MAJOR.BLEEDING<360, 
                            dapt.test.pred$DAYS.TIMI.MAJOR.BLEEDING, 360)
cb_x<-c()

source("C:/git/EVEREST-HTE/Final/Calibration.R")
cb_arr_prec <- Calib_ARR_dapt(dapt.test.pred, 5) 


model = glm(obs_rr1~pred_rr1, data=cb_arr_prec)
slope<-model$coefficients[[2]]#ideally 1
intercept<-model$coefficients[[1]]#ideally 0
summary(model)

Predicted.ARR<-cb_arr_prec[,1]
Observed.ARR<-cb_arr_prec[,3]

cb_arr_prec<-data.frame(Predicted.ARR, Observed.ARR)

plot(cb_arr_prec,pch=19, col="darkblue", 
     xlab="Predicted Risk-Probability", 
     ylab="Observed Risk-Probability",
     xlim=c(-0.1, 0.1), 
     ylim=c(-0.1, 0.1), 
     main="PRECISE DAPT SCORE")

abline(a=0,b=1,col="gray",lty=6)
abline(a=intercept, b=slope, col="darkblue", lty=5, lwd=1.8)

text(x=-0.04, y= 0.08, 
     labels = c(paste0("Slope: ", round(slope, digit=2), ",  Intercept: ", 
                                         round(intercept, digit=3))))

blp.summary <- lmtest::coeftest(model,
                                vcov = sandwich::vcovCL,
                                type = "HC3")

blp.summary


plt_arr_dapt<-ggplot(data=cb_arr_prec, aes(x=Predicted.ARR, y=Observed.ARR)) + 
  geom_point( size=2, colour="darkblue") + # shape 15: solid square
  ggtitle("Calibration for PRECISE-DAPT SCORE")+
  theme_test()+
  theme(axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11))+
  geom_abline(intercept= intercept, slope=slope, color='darkblue', size = 0.5, linetype="dashed")+
  geom_abline(intercept= 0, slope=1, color='darkgray', size = 0.5, linetype="solid")



plt_arr_dapt<-plt_arr_dapt+
  coord_cartesian(xlim=c(-0.1, 0.1), ylim=c(-0.1, 0.1))#+
#annotate("text", x=0.03, y=0.08, label="Versicolor", size=7)

