# ====================================================================================================
# Implementation of the X-learners meta-algorithm for estimating CATE using
# logistic regression. Predictions are *not* evaluated out-of-bag, so there
# may be overfitting if training and validating on the same data.
#
# Kunzel et. al. 2018: http://arxiv.org/abs/1706.03461.
# ====================================================================================================


library(glmnet)
library(ranger)

train_linear_learner = function(X, W, Y, iptw) {
  
  tf0 = glm(Y ~ ., data = data.frame(X[W == 0,], Y = Y[W == 0]),
            weights = rep(iptw, nrow(X[W == 0,])), family = "binomial")
  yhat0 = predict(tf0, X[W == 1,], type = "response")
  lf1 = glm(Y ~ ., data = data.frame(Y = Y[W == 1] - yhat0, X[W == 1,]),
            weights = rep(iptw, nrow(X[W == 1,])))
  
  tf1 = glm(Y ~ ., data = data.frame(X[W == 1,], Y = Y[W == 1]),
            weights = rep(iptw, nrow(X[W == 1,])), family = "binomial")
  yhat1 = predict(tf1, X[W == 0,], type = "response")
  lf0 = glm(Y ~ ., data = data.frame(Y = yhat1 - Y[W == 0], X[W == 0,]),
            weights = rep(iptw, nrow(X[W == 0,])))
  
  return(list(lf0 = lf0, lf1 = lf1))
}



predict_linear_learner = function(X, W, estimate_propensities) {
  
  preds_1 = predict(lf1, X)
  preds_0 = predict(lf0, X)
  
  if (estimate_propensities) {
    propf = ranger(W ~ ., data = data.frame(X, W = W),
                   min.node.size = 10)
    ehat = propf$predictions
    preds = (1 - ehat) * preds_1 + ehat * preds_0
  } else {
    preds = 0.5 * preds_1 + 0.5 * preds_0
  }
  
  return(preds)
}
