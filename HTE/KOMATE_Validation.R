
dataFolder <- Sys.getenv("dataFolderHTE")
RWD<-read.csv("RCT_Merged_KOMATE.csv")

RWD$STUDY<-ifelse(RWD$STUDY == "KOMATE", 1, 0)
RWD$SMOKING<-ifelse(RWD$SMOKING == "Current" | RWD$SMOKING == "Former", 1, 
                    ifelse(RWD$SMOKING =="Never", 0, NA))

tb1=mytable(DAPT.SHORT.LONG ~ AGE+SEX+DYSLIPIDEMIA+SMOKING+DIABETES+
                              DIABETES.INSULIN+HTN+CKD+CHF+PRIOR.MI+PRIOR.PCI+
                              PRIOR.CABG+PRIOR.STROKE+SERUM.CREATININE+
                              PLATELET.COUNT+WBC.COUNT+BASELINE.HEMOGLOBIN+
                              CLINICAL.PRESENTATION+ACUTE.CORONARY.SYNDROME+
                              DISEASE.VESSEL.NO+STENT.VESSEL.NO+STENT.LESION.NO+
                              IMPLANTED.STENT.NO+MVD.YN+TOTAL.STENT.LENGTH+
                              MINIMAL.STENT.DIAMETER+PCI.LM.YN,data=RWD)


all_variables_names<-c("TIMI.MAJOR.BLEEDING", "DAYS.TIMI.MAJOR.BLEEDING", 
                       "CD.MI.ST.CVA", "DAYS.CD.MI.ST.CVA",
                       "AD.MI.ST.CVA.MAJOR.BLEEDING", 
                       "DAYS.AD.MI.ST.CVA.MAJOR.BLEEDING","AGE", "SEX", 
                       "DAPT.SHORT.LONG", "STUDY",	
                       "DIABETES", "DIABETES.INSULIN", "DYSLIPIDEMIA", 
                       "HTN", "PRIOR.MI", "PRIOR.PCI", "PRIOR.STROKE", 
                       "CLINICAL.PRESENTATION", "ACUTE.CORONARY.SYNDROME",
                       "IMPLANTED.STENT.NO", "MVD.YN", "TOTAL.STENT.LENGTH", 
                       "PCI.LM.YN", "BASELINE.HEMOGLOBIN", "WBC.COUNT")


rwd_df<-RWD %>%dplyr::select(one_of(all_variables_names))


miss<-c()
Miss_NA<-c()
name<-c()

for ( i in 1:length(all_variables_names)) {
  
  nm<-all_variables_names[i]
  miss[[i]]<-table(is.na(rwd_df[,nm]))[1]
  Miss_NA1<-ifelse(rwd_df[,nm] == "NaN.NaN" | rwd_df[,nm] == "ND",1,0)
  Miss_NA[[i]]<-table(Miss_NA1)[1]
  
  name[[i]]<-nm
  
  i+1
}

miss_chk<-cbind(miss, Miss_NA, name)
miss_var<-subset(miss_chk, as.numeric(miss) != as.numeric(Miss_NA))


# Rename variables
rwd_df <- rwd_df %>% dplyr::rename(Y=TIMI.MAJOR.BLEEDING , W=DAPT.SHORT.LONG)
rwd_df$W2<-ifelse(rwd_df$W==1,0,1)
rwd_df$SEX2<-ifelse(rwd_df$SEX==1,0,1)
rwd_df$PRESENT_MI<-ifelse(rwd_df$CLINICAL.PRESENTATION == 3,1,0)

rwd_df<-rwd_df[,!names(rwd_df) %in% c("W")]


#"DYSLIPIDEMIA", "HTN", "PRIOR.MI", #"PRIOR.PCI",,"PCI.LM.YN",
covariate_names<-c("AGE", "SEX2", "DIABETES", "DIABETES.INSULIN", 
                   "PRIOR.STROKE", "ACUTE.CORONARY.SYNDROME",
                   "IMPLANTED.STENT.NO", "MVD.YN",
                   "TOTAL.STENT.LENGTH", "BASELINE.HEMOGLOBIN", 
                   "WBC.COUNT", "PRESENT_MI")


# Converting all columns to numerical
rwd_df <- data.frame(lapply(rwd_df, function(x) as.numeric(as.character(x))))
sapply(rwd_df, FUN = "class")

rwd_df_f<-na.omit(rwd_df)


# Prediction for Bleeding
set.seed(1234)
validatation_fit<-predict_x_learner(rwd_df_f[,covariate_names], 
                                    as.numeric(rwd_df_f$W2), TRUE, TRUE)

summary(validatation_fit)

rwd_test<-cbind(rwd_df_f, validatation_fit)

rwd_test$DAPT=ifelse(rwd_test$W2==1,"Short", "Long")
rwd_test$DAPT2=as.factor(rwd_test$DAPT)
rwd_test$Time=ifelse(rwd_test$DAYS.TIMI.MAJOR.BLEEDING<360, 
                     rwd_test$DAYS.TIMI.MAJOR.BLEEDING, 360)
rwd_test$class=ifelse(rwd_test$validatation_fit<=0, "benefit", "no-benefit")
rwd_test$class2=ifelse(rwd_test$validatation_fit<=0, 1, 0)
  
table(rwd_test$class)

tb2=mytable(class ~ AGE+SEX+DIABETES+DIABETES.INSULIN+PRIOR.STROKE+
                    WBC.COUNT+BASELINE.HEMOGLOBIN+ACUTE.CORONARY.SYNDROME+ 
                    IMPLANTED.STENT.NO+MVD.YN+TOTAL.STENT.LENGTH+PCI.LM.YN+ 
                    PRESENT_MI, data=rwd_test)

tmp<-rwd_test

for(n_col in c(1, 3, 5, 7, 9:18, 20, 22, 25:28)){
  tmp[, n_col] = as.factor(tmp[, n_col])
}

vars <- c("AGE", "SEX2", "DIABETES", "DIABETES.INSULIN", 
          "DYSLIPIDEMIA", "HTN", "PRIOR.MI", "PRIOR.PCI",
          "PRIOR.STROKE", "ACUTE.CORONARY.SYNDROME",
          "IMPLANTED.STENT.NO", "MVD.YN","PCI.LM.YN", 
          "PRESENT_MI", "TOTAL.STENT.LENGTH", 
          "BASELINE.HEMOGLOBIN", "WBC.COUNT", "W2")

## Construct a table
tabUnmatched <- CreateTableOne(vars = vars, strata = "class", 
                               data = tmp, test = FALSE)

## Show table with SMD
print(tabUnmatched, smd = TRUE)



#survival fitting
fit1<-survfit(Surv(Time, Y) ~ class, data = rwd_test)

survplot1<-ggsurvplot(fit1,
                      pval = T,
                      censor=F,
                      break.time.by=90, 
                      break.y.by=0.01,
                      pval.method=TRUE,
                      risk.table.col = "strata", 
                      ggtheme = theme_bw(), 
                      palette = c("#7CAE00", "#FF6600"),
                      fun = "cumhaz",
                      xlim = c(0, 360), ylim = c(0, 0.1),
                      conf.int = TRUE,
                      risk.table = TRUE,
                      title = "Classification for Benefit/No-benefit group",
                      fontsize=4.0,     
                      risk.table.pos="out",
                      risk.table.y.text.col=FALSE,
                      legend.labs = c("Benefit(Short)", "No-benefit"),
                      legend.title=""
                      
) 
ggpar(survplot1,
      font.main = c(13, "bold.italic", "darkblue"),
      font.x = c(11, "bold", "black"),
      font.y = c(11, "bold", "black"),
      font.tickslab = c(11, "plain", "black"), 
      font.legend = c(11, "bold")
)
  
km_diff1 <- survdiff(Surv(Time, Y) ~ class, data = rwd_test);summary(km_diff1)
cox.out1<-coxph(Surv(Time,Y) ~ class, data = rwd_test);summary(cox.out1)$conf.int
  
