# ===========================================================================================
# Implementation of the X-learners meta-algorithm for estimating CATE using
# random forests. Predictions are evaluated using learners for which each
# observation is out-of-bag.
#
# Kunzel et. al. 2018: http://arxiv.org/abs/1706.03461.
#
# ===========================================================================================


library(ranger)


get_mapping_to_full_dataset = function(X, W, assignment) {
  
  mapping = vector("numeric", length = dim(X)[1])
  count = 1
  for (i in 1:dim(X)[1]) {
    if (W[i] == assignment) {
      mapping[i] = count
      count = count + 1
    }
  }
  return(mapping)
}

get_oob_predictions = function(X, forest, mapping) {
  
  raw_preds = predict(forest, X, predict.all = TRUE)$predictions
  final_preds = vector("numeric", length = dim(X)[1])
  inbag_counts = forest$inbag.counts
  
  for (i in 1:dim(X)[1]) {
    if (mapping[i] == 0 || i > length(mapping)) {
      final_preds[i] = mean(raw_preds[i,])
    } else {
      temp_preds = vector("list", length = forest$num.trees)
      for (j in 1:forest$num.trees) {
        if (inbag_counts[j][[1]][mapping[i]] == 0) {
          temp_preds[[j]] = raw_preds[i,j]
        }
      }
      final_preds[i] = mean(unlist(Filter(is.numeric, temp_preds)))
    }
  }
  return(final_preds)
}


train_x_learner <- function(X, W, Y, iptw, save = TRUE) {
  
  tf0 = ranger(Y ~ ., data = data.frame(X[W == 0,], Y = Y[W == 0]),
               num.trees = 2500, min.node.size = 50,
               case.weights = rep(iptw, nrow(X[W == 0,])))
  yhat0 = predict(tf0, X[W == 1,])$predictions
  xf1 = ranger(Y ~ ., data = data.frame(Y = Y[W == 1] - yhat0, X[W == 1,]),
               keep.inbag = TRUE, num.trees = 2500, min.node.size = 50,
               case.weights = rep(iptw, nrow(X[W == 1,])), importance = "impurity")
  mapping1 = get_mapping_to_full_dataset(X, W, 1)
  
  tf1 = ranger(Y ~ ., data = data.frame(X[W == 1,], Y = Y[W == 1]),
               num.trees = 2500, min.node.size = 50,
               case.weights = rep(iptw, nrow(X[W == 1,])))
  yhat1 = predict(tf1, X[W == 0,])$predictions
  xf0 = ranger(Y ~ ., data= data.frame(Y = yhat1 - Y[W == 0], X[W == 0,]),
               keep.inbag = TRUE, num.trees = 2500, min.node.size = 50,
               case.weights = rep(iptw, nrow(X[W == 0,])), importance = "impurity")
  mapping0 = get_mapping_to_full_dataset(X, W, 0)
  
  if (save) {
    saveRDS(xf0, "/xf0.rds")
    saveRDS(xf1, "/xf1.rds")
  }
  
  return(list(xf0 = xf0, xf1 = xf1, mapping0 = mapping0, mapping1 = mapping1))
}

predict_x_learner <- function(X, W, estimate_propensities, predict_oob) {
  
  if (predict_oob) {
    preds_1 = get_oob_predictions(X, train_x_learner_fit$xf1, train_x_learner_fit$mapping1)
    preds_0 = get_oob_predictions(X, train_x_learner_fit$xf0, train_x_learner_fit$mapping0)
  } else {
    preds_1 = predict(train_x_learner_fit$xf1, X)$predictions
    preds_0 = predict(train_x_learner_fit$xf0, X)$predictions
  }
  
  if (estimate_propensities) {
    propf = ranger(W ~ ., data = data.frame(X, W = W),
                   min.node.size = 1)
    ehat = propf$predictions
    preds = (1 - ehat) * preds_1 + ehat * preds_0
  } else {
    preds = 0.5 * preds_1 + 0.5 * preds_0
  }
  
  return(preds)
}



