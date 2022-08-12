

Var_sel <- function(df, cov, method, seed_num){  
  
  
  data<-df[,c(cov,"Y", "W2")]
  
  for(n_col in c(2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 17, 18)){
    data[, n_col] = as.factor(data[, n_col])
  }
  
  
  if(method=="boruta"){ 
    
    ## Boruta
    ## Feature ranking and selection algorithm based on random forests algorithm
    
    library(Boruta)

    # Perform Boruta search    
    set.seed(seed_num)
    boruta_output <- Boruta(Y ~ ., data=data, doTrace=0)  
    names(boruta_output)
    boruta_signif <- getSelectedAttributes(boruta_output, withTentative = TRUE)
    print(boruta_signif) 
    
    
    roughFixMod <- TentativeRoughFix(boruta_output)
    boruta_signif <- getSelectedAttributes(roughFixMod)
    print(boruta_signif)
    
    
    # Variable Importance Scores
    imps <- attStats(roughFixMod)
    imps2 = imps[imps$decision != 'Rejected', c('meanImp', 'decision')]
    head(imps2[order(-imps2$meanImp), ])  # descending sort
    
    
    # Plot variable importance
    # plot(boruta_output, cex.axis=.7, las=2, xlab="", main="Variable Importance") 
    
    return(boruta_signif)
    
    
  } else if (method=="lasso"){
    

    x <- data %>% dplyr::select(cov) %>% data.matrix()
    y <- as.double(as.matrix(ifelse(data$Y==0, 0, 1)))
    
   
    library(glmnet)
    
    # Fit the LASSO model (Lasso: Alpha = 1)
    set.seed(100)
    cv.lasso <- cv.glmnet(x, y, family='binomial', alpha=1, parallel=TRUE, standardize=TRUE, type.measure='auc')

    # Results
    plot(cv.lasso)
    
    # plot(cv.lasso$glmnet.fit, xvar="lambda", label=TRUE)
    cat('Min Lambda: ', cv.lasso$lambda.min, '\n 1Sd Lambda: ', cv.lasso$lambda.1se)
    df_coef <- round(as.matrix(coef(cv.lasso, s=cv.lasso$lambda.min)), 2)
    
    # See all contributing variables
    var_Imp<-df_coef[df_coef[, 1] != 0, ]
    
    return(var_Imp)
    
    
  } else if (method=="rf"){
    
    
    library(caret)
    
    set.seed(seed_num)
    rPartMod <- train(as.factor(Y) ~ ., data=data, method="RRF")
    rpartImp <- varImp(rPartMod)
    
    return(rpartImp)
    
    
  } else if (method=="stepwise"){
    
    # Define base intercept only model
    base.mod <- glm(as.factor(Y) ~ 1 , family = "binomial", data=data[,c(cov, "Y")])  
    
    # Full model with all predictors
    all.mod <- glm(as.factor(Y) ~ . , family = "binomial", data=data[,c(cov, "Y")])  
    
    # Perform step-wise algorithm. direction='both' implies both forward and backward stepwise
    stepMod <- stats::step(base.mod, scope = list(lower = base.mod, upper = all.mod), direction = "both", trace = 0, steps = 1000)  
    
    # Get the shortlisted variable.
    shortlistedVars <- names(unlist(stepMod[[1]])) 
    shortlistedVars <- shortlistedVars[!shortlistedVars %in% "(Intercept)"] # remove intercept
    
    # Show
    return(shortlistedVars)
    
    
  } else if (method=="annealing"){
    
    # Define control function
    sa_ctrl <- safsControl(functions = rfSA,
                           method = "repeatedcv",
                           repeats = 3,
                           improve = 5) # n iterations without improvement before a reset
    
    
    # Genetic Algorithm feature selection
    set.seed(seed_num)
    
    sa_obj_w <- safs(x=data[,cov], 
                     y=data$Y,
                     safsControl = sa_ctrl)
    
    sa_obj_w_var<- sa_obj_w$optVariables
    
    
    # Show
    return(sa_obj_w_var)
    
    
  } else if (method=="dalex"){
    
    
    library(randomForest)
    library(DALEX)
    
    data_tmp<-data[,-19]
    
    
    # Train random forest model
    set.seed(seed_num)
    rf_mod <- randomForest(Y ~ ., data=data_tmp, ntree=10)
    rf_mod
    
    
    # Variable importance with DALEX
    explained_rf <- explain(rf_mod, data=data_tmp[,cov], y=as.numeric(data_tmp$Y))
    
    # Get the variable importances
    varimps = variable_importance(explained_rf, type='raw')
    
    return(varimps)
    # plot(varimps)
    
    
  } else {
    
    error<-c("error: invalid method")
    
    return(error)
    
  }
  
}


method<-c("boruta", "lasso", "rf", "stepwise", "annealing", "dalex")
cov_list<-c()

n<-6

for(i in (1:n)){
  
  method[i]
  
  assign(method[i], Var_sel(df_train, covariate_names, method[i], 100))
  
  cov_list[[i]]<-get(method[i])
  
}

names(cov_list)<-method[1:n]

print(cov_list)



