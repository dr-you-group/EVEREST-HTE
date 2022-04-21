# ----------------------------------------------------------------------------
#   Interpretation of X-learner with Random Forests
# ----------------------------------------------------------------------------

library(tidyverse)
library(gridExtra)
library(gsubfn)
library(ranger)
library(pdp)



xf0 = readRDS("/xf0.rds")
xf1 = readRDS("/xf1.rds")


#Relative importance plot
plot_importance = function(xf0, xf1, cols) {
  
  importances = (importance(xf0) + importance(xf1)) / 2
  
  ggplot(data = tibble(col = cols, imp = importances), aes(x=col, y=imp)) +
    geom_bar(stat = "identity", width=0.75) +
    labs(y = "Relative importance",
         x = "Covariate") +
    coord_flip() +
    theme_light()
}


plot_pdp = function(xf0, xf1, col_idx) {
  
  subsample = sample_n(data, 1000)
  pdp_0 = as_tibble(partial(xf0, train = subsample, pred.var = col_idx))
  pdp_1 = as_tibble(partial(xf1, train = subsample, pred.var = col_idx))
  pdp = tibble(val = pdp_0[,1][[1]], yhat = -(pdp_0$yhat + pdp_1$yhat) / 2)#ARR: Average Risk Reduction
  ggplot(pdp, aes(x = val, y = yhat)) +
    geom_point() +
    theme_light() +
    ylim(-0.1, 0.1) +
    labs(x = col_idx, y = "Predicted ARR")
}


plot_importance(xf0, xf1, cols)



pdps = c()



for (i in 1:(length(cols))) {
  pdps[[i]] = plot_pdp(xf0, xf1, covariate_names[i])
}

grid.arrange(pdps[[1]], pdps[[2]], pdps[[3]], pdps[[4]], pdps[[5]], pdps[[6]],
             pdps[[7]], pdps[[8]], pdps[[9]], pdps[[10]], pdps[[11]],
             pdps[[12]], pdps[[13]], pdps[[14]], ncol = 4)
