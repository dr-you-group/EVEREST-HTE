library(dplyr)       # Data manipulation (0.8.0.1)
library(fBasics)     # Summary statistics (3042.89)
library(grf)         # Generalized random forests (0.10.2)
library(car)         # linear hypothesis testing for causal tree (3.0-2)
library(remotes)     # Install packages from github (2.0.1)
library(readr)       # Reading csv files (1.3.1)
library(tidyr)       # Database operations (0.8.3)
library(tibble)      # Modern alternative to data frames (2.1.1)
library(knitr)       # RMarkdown (1.21)
library(kableExtra)  # Prettier RMarkdown (1.0.1)
library(ggplot2)     # general plotting tool (3.1.0)
library(haven)       # read stata files (2.0.0)
library(aod)         # hypothesis testing (1.3.1)
library(evtree)      # evolutionary learning of globally optimal trees (1.0-7)
library(estimatr)    # simple interface for OLS estimation w/ robust std errors ()


#############################################################################################################
############################################################################################################

all_variables_names<-c("AD.MI.ST.CVA.MAJOR.BLEEDING", "DAYS.AD.MI.ST.CVA.MAJOR.BLEEDING", "DAPT.SHORT.LONG", "STUDY",	
                       "AGE", "SEX", 
                       "DIABETES", "DIABETES.INSULIN", 
                       "DYSLIPIDEMIA", 	
                       "HTN", "PRIOR.MI", "PRIOR.PCI", #"PRIOR.CABG", 	
                       "PRIOR.STROKE", #"CLINICAL.PRESENTATION", 
                       "ACUTE.CORONARY.SYNDROME","IMPLANTED.STENT.NO", "MVD.YN",
                       "TOTAL.STENT.LENGTH", "PCI.LM.YN", "BASELINE.HEMOGLOBIN", "WBC.COUNT")


df<-RCT %>%select(one_of(all_variables_names))


# Rename variables
df <- df %>% rename(Y=AD.MI.ST.CVA.MAJOR.BLEEDING, W=DAPT.SHORT.LONG)

df$W2<-ifelse(df$W==1,1,0)
table(df$W);table(df$W2)

df$SEX2<-ifelse(df$SEX==1,0,1)
table(df$SEX);table(df$SEX2)


df<-df[,!names(df) %in% c("W")]

covariate_names<-c("AGE", "SEX2", 
                   "DIABETES", "DIABETES.INSULIN", 
                   "DYSLIPIDEMIA", 	
                   "HTN", "PRIOR.MI", "PRIOR.PCI",
                   "PRIOR.STROKE", "ACUTE.CORONARY.SYNDROME",
                   "IMPLANTED.STENT.NO", "MVD.YN",
                   "TOTAL.STENT.LENGTH", "BASELINE.HEMOGLOBIN", "WBC.COUNT")#, "PCI.LM.YN"




# Converting all columns to numerical
df <- data.frame(lapply(df, function(x) as.numeric(as.character(x))))


#delete missing value for hemoglobin
df_Missi<-subset(df,is.na(df$WBC.COUNT)==FALSE & is.na(df$BASELINE.HEMOGLOBIN)==FALSE)

table(df_Missi$STUDY)

table(df_Missi$STUDY, df_Missi$Y)

n <- dim(df_Missi)[1]
set.seed(1)


df_train <- df_Missi[df_Missi$STUDY!=3,]
df_test <- df_Missi[df_Missi$STUDY==3,]#TICO



# Make a data.frame containing summary statistics of interest
summ_stats <- fBasics::basicStats(df)
summ_stats <- as.data.frame(t(summ_stats))


# Rename some of the columns for convenience
summ_stats <- summ_stats[c("Mean", "Stdev", "Minimum", "1. Quartile", "Median",  "3. Quartile", "Maximum")] %>% 
  rename("Lower quartile" = '1. Quartile', "Upper quartile"= "3. Quartile")

