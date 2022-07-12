# outcome for NACE, MACE

aply=function(name) {
  
  var<-df[,c(name)]
  
  df<-RCT %>%dplyr::select(one_of(all_variables_names))
  
  df_ex<-c()
  df_Missi_ex<-c()
  df_train_ex<-c()
  df_test_ex<-c()
  
  
  df_ex <- df %>% rename(Y=name, W=DAPT.SHORT.LONG)
  
  df_ex$W2<-ifelse(df_ex$W==1,1,0)
  df_ex$SEX2<-ifelse(df_ex$SEX==1,0,1)
  df_ex$PRESENT_MI<-ifelse(df_ex$CLINICAL.PRESENTATION == 3,1,0)
  df_ex<-df_ex[,!names(df) %in% c("W")]
  
  # Converting all columns to numerical
  patient_id<-df$STUDY.PATIENT.ID
  df_ex <- data.frame(patient_id, lapply(df_ex[c(2:28)], function(x) as.numeric(as.character(x))))
  
  df_Missi_ex<-subset(df_ex,is.na(df_ex$WBC.COUNT)==FALSE & is.na(df_ex$BASELINE.HEMOGLOBIN)==FALSE)
  
  n <- dim(df_Missi_ex)[1]
  
  set.seed(1)
  
  df_train_ex <- df_Missi_ex[df_Missi_ex$STUDY!=3,]
  df_test_ex <- df_Missi_ex[df_Missi_ex$STUDY==3,]#TICO
  
  
  
  class_ex<-surv_d$test[,c("class")]
  test_ex<-cbind(df_test_ex, class_ex)
  
  test_ex$DAPT<-ifelse(test_ex$W2==1,"Short", "Long")
  test_ex$DAPT2<-as.factor(test_ex$DAPT)
  
  Time<- paste0("DAYS.", paste(name))
  
  test_ex$Time<-ifelse(test_ex[,c(Time)]<360, test_ex[,c(Time)], 360)
  Time<-ifelse(test_ex[,c(Time)]<360, test_ex[,c(Time)], 360)
  
  fit1<-survfit(Surv(Time, Y) ~ class_ex, data = test_ex)
  fit3<-survfit(Surv(Time, Y) ~ DAPT, data = subset(test_ex, class_ex=="benefit"))
  fit4<-survfit(Surv(Time, Y) ~ DAPT, data = subset(test_ex, class_ex=="no-benefit"))
  
  return(list(fit1=fit1, fit3=fit3, fit4=fit4, test_ex=test_ex, name=name))
}
