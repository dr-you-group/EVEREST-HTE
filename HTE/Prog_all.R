pacman::p_load("dplyr", "fBasics", "car", "remotes", "readr", "tidyr", "tibble", 
               "ggplot2", "aod", "evtree", "estimatr", "tidyverse", "gridExtra", 
               "gsubfn", "pdp", "ranger", "survival", "survminer", "moonBook", 
               "tableone", "Matching", "Boruta", "caret", "glmnet", "randomForest",
               "DALEX", "forcats", "rms", "survRM2")


dataFolder <- Sys.getenv("dataFolderHTE")

#data(raw data)
RCT <- read.csv(file.path(dataFolder,"EVEREST_ver3.2.csv")) 

#combine all names
all_variables_names<-c("STUDY.PATIENT.ID", "TIMI.MAJOR.BLEEDING", "DAYS.TIMI.MAJOR.BLEEDING", 
                       "CD.MI.ST.CVA", "DAYS.CD.MI.ST.CVA",
                       "AD.MI.ST.CVA.MAJOR.BLEEDING", "DAYS.AD.MI.ST.CVA.MAJOR.BLEEDING",
                       "DAPT.SHORT.LONG", "STUDY", "AGE", "SEX", "DIABETES", 
                       "DIABETES.INSULIN", "DYSLIPIDEMIA", "HTN", "PRIOR.MI", 
                       "PRIOR.PCI",	"PRIOR.STROKE", "CLINICAL.PRESENTATION", 
                       "ACUTE.CORONARY.SYNDROME","IMPLANTED.STENT.NO", "MVD.YN",
                       "TOTAL.STENT.LENGTH", "PCI.LM.YN", "BASELINE.HEMOGLOBIN", 
                       "WBC.COUNT")


df<-RCT %>% dplyr::select(one_of(all_variables_names))


# Rename variables
df <- df %>% rename(Y=TIMI.MAJOR.BLEEDING, W=DAPT.SHORT.LONG) #W=1: Short, W=2: Long

df$W2<-ifelse(df$W==1,1,0) #W2=0: long, W2=1: short
df$SEX2<-ifelse(df$SEX==1,0,1)
df$PRESENT_MI<-ifelse(df$CLINICAL.PRESENTATION == 3,1,0)

df<-df[,!names(df) %in% c("W")]


covariate_names<-c("AGE", "SEX2", "DIABETES", "DIABETES.INSULIN", 
                   "DYSLIPIDEMIA", "HTN", "PRIOR.MI", "PRIOR.PCI",
                   "PRIOR.STROKE", "ACUTE.CORONARY.SYNDROME",
                   "IMPLANTED.STENT.NO", "MVD.YN","PCI.LM.YN",
                   "TOTAL.STENT.LENGTH", "BASELINE.HEMOGLOBIN", 
                   "WBC.COUNT", "PRESENT_MI")



# Converting all columns to numerical
patient_id<-df$STUDY.PATIENT.ID
df <- data.frame(patient_id, lapply(df[c(2:28)], function(x) as.numeric(as.character(x))))

for(n_col in c(10:19, 21, 23, 26, 27, 28 )){
  df[, n_col] = as.factor(df[, n_col])
}

df_Missi<-subset(df,is.na(df$WBC.COUNT)==FALSE & is.na(df$BASELINE.HEMOGLOBIN)==FALSE)


n <- dim(df_Missi)[1]


set.seed(1)
df_train <- df_Missi[df_Missi$STUDY!=3,]
df_test <- df_Missi[df_Missi$STUDY==3,]#TICO



# predicte ARR through the X-learner
source("C:/git/EVEREST-HTE/Final/X-learner.R")

set.seed(1111)
train_x_learner_fit <- train_x_learner(df_train[,covariate_names], df_train$W2, 
                                       df_train$Y, ipcw=0.5) ##ipcw=rep(0.5, nrow(df_train))
predict_x_learner_fit <- predict_x_learner(df_test[,covariate_names], as.numeric(df_test$W2), 
                                           estimate_propensities = FALSE, predict_oob = TRUE)

test<-cbind(df_test[,-1], predict_x_learner_fit)

test$DAPT<-ifelse(test$W2==1,"Short", "Long")
test$Time<-ifelse(test$DAYS.TIMI.MAJOR.BLEEDING<360, test$DAYS.TIMI.MAJOR.BLEEDING, 360)
test$class<-ifelse(test$predict_x_learner_fit<=0, "benefit", "no-benefit")


# Interpretation of X-learner with Random Forests
source("C:/git/EVEREST-HTE/Final/Interpretation.R")

xf0 = readRDS(file.path(dataFolder, "xf0.rds"))
xf1 = readRDS(file.path(dataFolder, "xf1.rds"))

#Relative importance plot
plot_importance(xf0, xf1, df_train, covariate_names, 16)


pdps <- c()

cols = toupper(c(colnames(df_train[,covariate_names])))

for (i in 1:(length(cols))) {
  pdps[[i]] = plot_pdp(xf0, xf1, df_train, cols[i])
}

grid.arrange(pdps[[1]], pdps[[2]], pdps[[3]], pdps[[4]], pdps[[5]], pdps[[6]],
             pdps[[7]], pdps[[8]], pdps[[9]], pdps[[10]], pdps[[11]], pdps[[12]], 
             pdps[[13]], pdps[[14]], pdps[[15]], pdps[[16]], pdps[[17]], ncol = 5)


# survival plot

surv_func = function(test){
  
  test$DAPT<-ifelse(test$W2==1,"Short", "Long")
  test$DAPT2<-as.factor(test$DAPT)
  test$Time<-ifelse(test$DAYS.TIMI.MAJOR.BLEEDING<360, test$DAYS.TIMI.MAJOR.BLEEDING, 360)
  test$class<-ifelse(test$predict_x_learner_fit<=0, "benefit", "no-benefit")
  test$class2<-ifelse(test$predict_x_learner_fit<=0, 1, 0)
  Time<-ifelse(test$DAYS.TIMI.MAJOR.BLEEDING<360, test$DAYS.TIMI.MAJOR.BLEEDING, 360)
  
  
  fit1<-survfit(Surv(Time, Y) ~ class, data = test)
  fit2<-survfit(Surv(Time, Y) ~ DAPT, data = subset(test, class=="benefit"))
  fit3<-survfit(Surv(Time, Y) ~ DAPT, data = subset(test, class=="no-benefit"))
  
  
  
  summary(fit1);summary(fit2);summary(fit3)
  
  return(list(test=test, fit1=fit1, fit2=fit2, fit3=fit3, test=test))
  
}      

surv_d<-surv_func(test)

survplot_bld = c()

form<-c()
dat<-c()

form[[1]]<-formula("Surv(Time, Y) ~ class")
dat[[1]]<-surv_d$test

form[[2]]<-formula("Surv(Time, Y) ~ DAPT")
dat[[2]]<-subset(surv_d$test, class=="benefit")

form[[3]]<-formula("Surv(Time, Y) ~ DAPT")
dat[[3]]<-subset(surv_d$test, class=="no-benefit")


#Log-Lank Test
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

for (i in 1:3) {
  
  title.nm<-c("TIMI.MAJOR.BLEEDING", 
              "TIMI.MAJOR.BLEEDING for Benefit", 
              "TIMI.MAJOR.BLEEDING for No-benefit")
  y.max<-c(0.1, 0.1, 0.1)
  
  legend.lab<-c()
  legend.lab[[1]]<-c("No-benefit", "Benefit(Short)")
  legend.lab[[2]]<-c("12mo", "<6mo")
  legend.lab[[3]]<-c("12mo", "<6mo")
  
  
  palette<-c()
  palette[[1]]<-c("#FF6600", "#7CAE00")
  palette[[2]]<-"circulation"
  palette[[3]]<-"circulation"
  
  #plotting
  survplot_bld[[i]]<-ggsurvplot(survfit(form[[i]], data =dat[[i]]),
                                conf.int = T,
                                pval = F,
                                censor=F,
                                break.time.by=90, 
                                break.y.by=0.01,
                                pval.method=TRUE,
                                risk.table.col = "strata",
                                ggtheme = theme_bw(), 
                                palette = palette[[i]],
                                fun = "cumhaz",
                                xlim = c(0, 360), ylim = c(0, y.max[i]),
                                risk.table = TRUE,
                                tables.theme = theme_cleantable(),
                                title = title.nm[i],
                                fontsize=3.0,     
                                risk.table.pos="out",
                                risk.table.y.text.col=FALSE,
                                legend.labs = legend.lab[[i]],
                                legend.title=""
                                
  ) 
  
  survplot_bld[[i]]$plot <- survplot_bld[[i]]$plot + 
    ggplot2::annotate(
      "text",
      x=45, y=0.075,
      vjust = 1, hjust = 0.5,
      label=ifelse(i==1, paste0("Log-rank test \n p-value = ", round(Result[[1]][i], digit=3),
                                " \n\n HR(benefit) \n ", round(Result[[2]][i], digit=3), 
                                " (", round(Result[[4]][i], digit=3), 
                                " - ", round(Result[[5]][i], digit=3), ")"),
                         paste0("Log-rank test \n p-value = ", round(Result[[1]][i], digit=3),
                          " \n\n HR(short) \n ", round(Result[[2]][i], digit=3), 
                          " (", round(Result[[4]][i], digit=3), 
                          " - ", round(Result[[5]][i], digit=3), ")")),
      size=3.5 
      
    ) + 
    
    theme(
      
      plot.margin = unit(c(1, 5.5, 5.5, -8), "points")
          
          )
  
  
  survplot_bld[[i]]$table <- ggpubr::ggpar(survplot_bld[[i]]$table,
                                           font.title = list(size = 11, color = "black")
  )
  
  survplot_bld[[i]]$table <- survplot_bld[[i]]$table + 
    theme(
      plot.title = element_text(hjust = 0, vjust = 0),
      plot.margin = unit(c(1, 5.5, 5.5, -8), "points"))
}



survplot_bld



# demography
tb1=mytable(class ~ AGE+SEX2+DYSLIPIDEMIA+DIABETES+DIABETES.INSULIN+HTN+PRIOR.MI+
                    PRIOR.PCI+PRIOR.STROKE+WBC.COUNT+BASELINE.HEMOGLOBIN+
                    IMPLANTED.STENT.NO+MVD.YN+TOTAL.STENT.LENGTH+PCI.LM.YN+
                    PRESENT_MI, data=surv_d$test)

for(n_col in c(1, 3, 5, 7, 9:18, 20, 22, 25:28)){
  surv_tmp<-surv_d$test
  surv_tmp[, n_col]<- as.factor(surv_tmp[, n_col])
}


vars <- c("AGE", "SEX2", "DIABETES", "DIABETES.INSULIN", 
          "DYSLIPIDEMIA", "HTN", "PRIOR.MI", "PRIOR.PCI",
          "PRIOR.STROKE", "ACUTE.CORONARY.SYNDROME",
          "IMPLANTED.STENT.NO", "MVD.YN","PCI.LM.YN", 
          "PRESENT_MI", "TOTAL.STENT.LENGTH", 
          "BASELINE.HEMOGLOBIN", "WBC.COUNT", "W2")

# Construct a table
tabUnmatched <- CreateTableOne(vars = vars, 
                               strata = "class", 
                               data = surv_tmp, 
                               test = FALSE)

## Show table with SMD
print(tabUnmatched, smd = TRUE)



# herral's c-index for x-learner model
source("C:/git/EVEREST-HTE/Final/Calibration.R")
c_indx(surv_d$test, "class")


# calibration plot
## Form a calibration plot of predicted risk reduction, return the
## slope, intercept, predicted and observed ARRs.

test$Time<-ifelse(test$DAYS.TIMI.MAJOR.BLEEDING<360, 
                  test$DAYS.TIMI.MAJOR.BLEEDING, 360)
cb_x<-c()

cb_x <- Calib_arr_x(test, 5) 
cb_x<-as.data.frame(cb_x)

model = glm(cb_x[,3]~cb_x[,1], data=cb_x)
slope<-model$coefficients[[2]]#ideally 1
intercept<-model$coefficients[[1]]#ideally 0
summary(model)

Predicted.ARR<-cb_x[,1]
Observed.ARR<-cb_x[,3]

cb_p<-data.frame(Predicted.ARR, Observed.ARR)
plot(cb_p,pch=19, col="darkgreen", xlab="Predicted ARR", ylab="Observed ARR",
     xlim=c(-0.08, 0.08), ylim=c(-0.08, 0.08), main="Machine Learning for X-learner")
abline(a=0,b=1,col="gray",lty=6)
abline(a=intercept, b=slope, col="darkgreen", lty=5, lwd=1.8)
text(x=-0.03, y= 0.07, 
     labels = c(paste0("Slope: ", round(slope, digit=2), ",  Intercept: ", 
                       round(intercept, digit=3))))


blp.summary <- lmtest::coeftest(model,
                                vcov = sandwich::vcovCL,
                                type = "HC3")
blp.summary


plt_arr_x<-ggplot(data=cb_p, aes(x=Predicted.ARR, y=Observed.ARR)) + 
     geom_point( size=2, colour="darkgreen") + 
     ggtitle("Calibration for X-learner model")+
     theme_test()+
     theme(axis.text.x = element_text(size=11),
           axis.text.y = element_text(size=11))+
     geom_abline(intercept= intercept, slope=slope, color='darkgreen', size = 0.5, linetype="dashed")+
     geom_abline(intercept= 0, slope=1, color='darkgray', size = 0.5, linetype="solid")
     #stat_smooth(method="glm", level=0.95, color="darkgreen", size = 0.5)
  

plt_arr_x<-plt_arr_x+coord_cartesian(xlim=c(-0.08, 0.08), ylim=c(-0.08, 0.08))


# calibration plot(for Risk-Probability with x-learner model)
test$pred_Y<-ifelse(test$predict_x_learner_fit<=0,1,0)
test <- mutate(test, quantile_rank = ntile(test$predict_x_learner_fit, 25))
test$predict_x_learner_fit2<-test$predict_x_learner_fit*1

cb_precise <- Calib_risk_x(test, 5)
cb_precise<-as.data.frame(cb_precise)

model = glm(obs_prr~pred_prr, data=cb_precise)
slope<-model$coefficients[[2]]#ideally 1
intercept<-model$coefficients[[1]]#ideally 0
summary(model)


plot(cb_precise$pred_prr,cb_precise$obs_prr, 
     pch=19, col="darkgreen", 
     xlab="Predicted Risk-Probability", 
     ylab="Observed Risk-Probability", 
     xlim=c(0, 0.1), ylim=c(0, 0.1), 
     main="Machine Learning for X-learner")

abline(a=0,b=1,col="gray",lty=6)
abline(a=intercept, b=slope, col="darkgreen", lty=5, lwd=1.8)
text(x=0.033, y= 0.095, 
     labels = c(paste0("Slope: ", round(slope, digit=2), ",  Intercept: ", 
                       round(intercept, digit=3))))



plt_risk_x<-ggplot(data=cb_precise, aes(x=pred_prr, y=obs_prr)) + 
  geom_point( size=2, colour="darkgreen") + # shape 15: solid square
  ggtitle("Calibration for PRECISE-DAPT SCORE")+
  theme_test()+
  theme(axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11))+
  geom_abline(intercept= intercept, slope=slope, color='darkgreen', size = 0.5, linetype="dashed")+
  geom_abline(intercept= 0, slope=1, color='darkgray', size = 0.5, linetype="solid")



plt_risk_x<-plt_risk_x+
  coord_cartesian(xlim=c(0, 0.1), ylim=c(0, 0.1))#+
#annotate("text", x=0.03, y=0.08, label="Versicolor", size=7)



# Apply Model for other outcome(NACE, MACE)

# NACE
source("C:/git/EVEREST-HTE/Final/Apply.R")

result_ex<-aply("AD.MI.ST.CVA.MAJOR.BLEEDING")

form<-c()
dat<-c()

survplot_NACE <- c()

form[[1]]<-formula("Surv(Time, Y) ~ class_ex")
dat[[1]]<-result_ex$test_ex

form[[2]]<-formula("Surv(Time, Y) ~ DAPT")
dat[[2]]<-subset(result_ex$test_ex, class_ex=="benefit")

form[[3]]<-formula("Surv(Time, Y) ~ DAPT")
dat[[3]]<-subset(result_ex$test_ex, class_ex=="no-benefit")

#Log-Lank Test
km_diff<-c()
km_pval<-c()
tmp<-c()
cox_out<-c()
cox<-c()

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


for (i in 1:3) {
  
  title.nm<-c("AD.MI.ST.CVA.MAJOR.BLEEDING", 
              "AD.MI.ST.CVA.MAJOR.BLEEDING for Benefit", 
              "AD.MI.ST.CVA.MAJOR.BLEEDING for No-benefit")
  
  y.max<-c(0.08, 0.08, 0.1)
  
  legend.lab<-c()
  legend.lab[[1]]<-c("No-benefit", "Benefit(Short)")
  legend.lab[[2]]<-c("12mo", "<6mo")
  legend.lab[[3]]<-c("12mo", "<6mo")
  
  
  palette<-c()
  palette[[1]]<-c("#FF6600", "#7CAE00")
  palette[[2]]<-"circulation"
  palette[[3]]<-"circulation"
  
  
  #plotting
  survplot_NACE[[i]]<-ggsurvplot(survfit(form[[i]], data =dat[[i]]),
                                 conf.int = T,
                                 pval = T,
                                 censor=F,
                                 break.time.by=90, 
                                 break.y.by=0.01,
                                 pval.method=TRUE,
                                 risk.table.col = "strata",
                                 ggtheme = theme_bw(), 
                                 palette = palette[[i]],
                                 fun = "cumhaz",
                                 ylab = "",
                                 xlim = c(0, 360), ylim = c(0, y.max[i]),
                                 risk.table = TRUE,
                                 title = title.nm[i],
                                 fontsize=3.0,     
                                 risk.table.pos="out",
                                 risk.table.y.text.col=FALSE,
                                 tables.theme = theme_cleantable(),
                                 legend.labs = legend.lab[[i]],
                                 legend.title=""
                                 
  ) 
  
  survplot_NACE[[i]]$plot <- survplot_NACE[[i]]$plot + 
    ggplot2::annotate(
      "text",
      x=45, y=0.075,
      vjust = 1, hjust = 0.5,
      label=ifelse(i==1, paste0("Log-rank test \n p-value = ", round(Result[[1]][i], digit=3),
                                " \n\n HR(benefit) \n ", round(Result[[2]][i], digit=3), 
                                " (", round(Result[[4]][i], digit=3), 
                                " - ", round(Result[[5]][i], digit=3), ")"),
                         paste0("Log-rank test \n p-value = ", round(Result[[1]][i], digit=3),
                          " \n\n HR(short) \n ", round(Result[[2]][i], digit=3), 
                          " (", round(Result[[4]][i], digit=3), 
                          " - ", round(Result[[5]][i], digit=3), ")")
                   ),
      size=3.5 
      
    ) + theme(
      
      plot.margin = unit(c(1, 5.5, 5.5, -8), "points")
      
      )
  
  
  survplot_NACE[[i]]$table <- ggpubr::ggpar(survplot_NACE[[i]]$table,
                                            font.title = list(size = 11, color = "black")
  )
  
  survplot_NACE[[i]]$table <- survplot_NACE[[i]]$table + 
    theme(
      plot.title = element_text(hjust = 0, vjust = 0),
      plot.margin = unit(c(1, 5.5, 5.5, -8), "points"))
  
}



survplot_NACE

# herral's c-index for x-learner model
c_indx(result_ex$test_ex, "class_ex")


# MACE
result_ex<-aply("CD.MI.ST.CVA")

form<-c()
dat<-c()

survplot_MACE<-c()

form[[1]]<-formula("Surv(Time, Y) ~ class_ex")
dat[[1]]<-result_ex$test_ex

form[[2]]<-formula("Surv(Time, Y) ~ DAPT")
dat[[2]]<-subset(result_ex$test_ex, class_ex=="benefit")

form[[3]]<-formula("Surv(Time, Y) ~ DAPT")
dat[[3]]<-subset(result_ex$test_ex, class_ex=="no-benefit")

#Log-Lank Test
km_diff<-c()
km_pval<-c()
tmp<-c()
cox_out<-c()
cox<-c()

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

for (i in 1:3) {
  
  title.nm<-c("CD.MI.ST.CVA", "CD.MI.ST.CVA for Benefit", "CD.MI.ST.CVA for No-benefit")
  y.max<-c(0.08, 0.08, 0.08)
  
  legend.lab<-c()
  legend.lab[[1]]<-c("No-benefit", "Benefit(Short)")
  legend.lab[[2]]<-c("12mo", "<6mo")
  legend.lab[[3]]<-c("12mo", "<6mo")
  
  
  palette<-c()
  palette[[1]]<-c("#FF6600", "#7CAE00")
  palette[[2]]<-"circulation"
  palette[[3]]<-"circulation"
  
  #plotting
  survplot_MACE[[i]]<-ggsurvplot(survfit(form[[i]], data =dat[[i]]),
                                 conf.int = T,
                                 pval = T,
                                 censor=F,
                                 break.time.by=90, 
                                 break.y.by=0.01,
                                 pval.method=TRUE,
                                 risk.table.col = "strata", 
                                 ggtheme = theme_bw(), 
                                 palette = palette[[i]],
                                 fun = "cumhaz",
                                 ylab="",
                                 xlim = c(0, 360), ylim = c(0, y.max[i]),
                                 risk.table = TRUE,
                                 title = title.nm[i],
                                 fontsize=3.0,     
                                 risk.table.pos="out",
                                 risk.table.y.text.col=FALSE,
                                 tables.theme = theme_cleantable(),
                                 legend.labs = legend.lab[[i]],
                                 legend.title=""
                                 
  ) 
  
  survplot_MACE[[i]]$plot <- survplot_MACE[[i]]$plot + 
    ggplot2::annotate(
      "text",
      x=45, y=0.075,
      vjust = 1, hjust = 0.5,
      label=ifelse(i==1, paste0("Log-rank test \n p-value = ", round(Result[[1]][i], digit=3),
                                " \n\n HR(benefit) \n ", round(Result[[2]][i], digit=3), 
                                " (", round(Result[[4]][i], digit=3), 
                                " - ", round(Result[[5]][i], digit=3), ")"),
                         paste0("Log-rank test \n p-value = ", round(Result[[1]][i], digit=3),
                          " \n\n HR(short) \n ", round(Result[[2]][i], digit=3), 
                          " (", round(Result[[4]][i], digit=3), 
                          " - ", round(Result[[5]][i], digit=3), ")")),
      size=3.5 
      
    ) + theme(
      
      plot.margin = unit(c(1, 5.5, 5.5, -8), "points")
      
      )
  
  
  survplot_MACE[[i]]$table <- ggpubr::ggpar(survplot_MACE[[i]]$table,
                                            font.title = list(size = 11, color = "black")
  )
  
  survplot_MACE[[i]]$table <- survplot_MACE[[i]]$table + 
    theme(
      plot.title = element_text(hjust = 0, vjust = 0),
      plot.margin = unit(c(1, 5.5, 5.5, -8), "points"))
}


survplot_MACE


# herral's c-index for x-learner model
c_indx(result_ex$test_ex, "class_ex")



