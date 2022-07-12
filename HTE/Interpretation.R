# ----------------------------------------------------------------------------
#   Interpretation of X-learner with Random Forests
# ----------------------------------------------------------------------------

#Relative importance plot
plot_importance = function(xf0, xf1, data, col, col_n) {
  
  df_train_cov<-data[,col]
  cols = tolower(c(colnames(df_train_cov)))
  features = sprintf("x%d", seq(0,col_n))
  colnames(df_train_cov) <- features
  
  importances = (ranger::importance(xf0) + ranger::importance(xf1)) / 2
  imp_plt<-tibble(col = col, imp = importances)
  
  imp_plt %>%
    mutate(col = fct_reorder(col, imp)) %>%
    ggplot(aes(x=col, y=imp)) +
    geom_bar(stat = "identity", width=0.75) +
    labs(y = "Relative importance",
         x = "Covariate") +
    coord_flip() +
    theme_light()
  
}


plot_pdp = function(xf0, xf1, data, col_idx) {
  
    set.seed(1111)
  subsample = dplyr::sample_n(data, 1000)
  pdp_0 = as_tibble(partial(xf0, train = subsample, pred.var = col_idx))
  pdp_1 = as_tibble(partial(xf1, train = subsample, pred.var = col_idx))
  pdp = tibble(val = pdp_0[,1][[1]], yhat = -(pdp_0$yhat + pdp_1$yhat) / 2) 
  ggplot(pdp, aes(x = val, y = yhat)) +
    geom_point() +
    theme_light() +
    ylim(-0.1, 0.1) +
    labs(x = col_idx, y = "Predicted ARR")
  
}


