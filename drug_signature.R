##Enet
x1 <- as.matrix(test_dd[, pre_var])
x2 <- as.matrix(test_dd$response)

set.seed(seed)
fit = cv.glmnet(x1, x2, alpha = alpha, nfolds = 10)
rs <- lapply(val_dd_list,function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit,type = 'link', newx = as.matrix(x[,-c(1,2)]), s = fit$lambda.min)))})


##GBM
set.seed(seed)
library(gbm)
?gbm()
fit <- gbm(formula = response~., data = est_dd, 
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = response~., data = est_dd, 
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 8)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})


##Lasso
x = model.matrix(response~.,trainDat[,-1])[,-1]
y = as.numeric(trainDat$response)-1

fit_lasso = glmnet(x,y,alpha = 1,
                   family = "binomial",
                   type.measure = "class")
cv_lasso = cv.glmnet(x,y,alpha = 1,
                     family = "binomial",
                     type.measure = "class",
                     nfolds = 100)

##RandomForest
fit_rf = randomForest(response~., data = trainDat[,-1],
                      importance = TRUE, ntree = 500)
set.seed(123)
ctrl <- trainControl(method = "cv",number = 100)
rf_model_cv = train(response~.,data = trainDat[,-1],
                    method = "rf",trControl = ctrl)

rf_pred = predict(rf_model_cv,newdata = testDat[,-1],method = "rf",
                  type = "raw")
rf_pred = as.character(rf_pred)
rf_pred = as.numeric(rf_pred)

rf_roc = roc(testDat$response,rf_pred)
auc_value = auc(rf_roc)

##Ridge
##Lasso
x = model.matrix(response~.,trainDat[,-1])[,-1]
y = as.numeric(trainDat$response)-1

fit_lasso = glmnet(x,y,alpha = 0,
                   family = "binomial",
                   type.measure = "class")
cv_lasso = cv.glmnet(x,y,alpha = 0,
                     family = "binomial",
                     type.measure = "class",
                     nfolds = 100)