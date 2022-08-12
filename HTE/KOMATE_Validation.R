
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
                   "WBC.COUNT", "PRESENT_MI", 
                   "DYSLIPIDEMIA", "HTN", "PRIOR.MI", "PRIOR.PCI","PCI.LM.YN")


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


  
#survival fitting
surv_func = function(data){
  
  data$DAPT<-ifelse(data$W2==1,"Short", "Long")
  data$DAPT2<-as.factor(data$DAPT)
  data$Time<-ifelse(data$DAYS.TIMI.MAJOR.BLEEDING<360, data$DAYS.TIMI.MAJOR.BLEEDING, 360)
  data$class<-ifelse(data$validatation_fit<=0, "benefit", "no-benefit")
  data$class2<-ifelse(data$validatation_fit<=0, 1, 0)
  Time<-ifelse(data$DAYS.TIMI.MAJOR.BLEEDING<360, data$DAYS.TIMI.MAJOR.BLEEDING, 360)
  
  
  fit1<-survfit(Surv(Time, Y) ~ class, data = data)
  summary(fit1)
  
  return(list(test=data, fit1=fit1))
  
}      

surv_d<-surv_func(rwd_test)

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
                                pval.method=F,
                                risk.table.col = "strata",
                                ggtheme=theme_classic2(base_size=12, base_family="Arial"), 
                                palette = palette[[i]],
                                fun = "cumhaz",
                                xlim = c(0, 360), ylim = c(0, y.max[i]),
                                risk.table = TRUE,
                                #tables.theme = theme_cleantable(),
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



KOMATE_outf<-survplot_bld[[1]]

KOMATE_out <- read_pptx()
KOMATE_out <- add_slide(KOMATE_out)
KOMATE_out <- ph_with(x=KOMATE_out, value = dml(code = print(KOMATE_outf)), 
                 location=ph_location_type(type="body"), bg = "transparent" )

print(KOMATE_out, target = "KOMATE_out.pptx")



rwd_tmp<-surv_d$test
table(rwd_tmp$class)

tb2=mytable(class ~ AGE+SEX+DIABETES+DIABETES.INSULIN+HTN+DYSLIPIDEMIA+
                    PRIOR.MI+PRIOR.STROKE+PRIOR.PCI+WBC.COUNT+
                    BASELINE.HEMOGLOBIN+ACUTE.CORONARY.SYNDROME+ 
                    IMPLANTED.STENT.NO+MVD.YN+TOTAL.STENT.LENGTH+PCI.LM.YN+ 
                    PRESENT_MI, data=rwd_tmp)

mycsv(tb2,file="komate.csv")

for(n_col in c(1, 3, 5, 8:18, 20, 22, 25:27, 29)){
  rwd_tmp[, n_col] = as.factor(rwd_tmp[, n_col])
}

vars <- c("AGE", "SEX", "DIABETES", "DIABETES.INSULIN", 
          "HTN", "DYSLIPIDEMIA", "PRIOR.MI", "PRIOR.STROKE",
          "PRIOR.PCI", "WBC.COUNT", "BASELINE.HEMOGLOBIN",
          "ACUTE.CORONARY.SYNDROME", "MVD.YN", "IMPLANTED.STENT.NO", 
          "TOTAL.STENT.LENGTH", "PCI.LM.YN", "PRESENT_MI")

## Construct a table
tabUnmatched <- CreateTableOne(vars = vars, strata = "class", 
                               data = rwd_tmp, test = FALSE)

## Show table with SMD
print(tabUnmatched, smd = TRUE)





