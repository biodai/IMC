##signature for all
library(tidyverse)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBostatust)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(miscTools)
library(compareC)
library(ggplot2)
library(ggsci)
library(tidyr)
library(ggbreak)
library(data.table)

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

mm <- list(TCGA = tcga,
           GSE42743 = gse42743, GSE41613 = gse41613)

genelist1 = Reduce(intersect,list(genelist1,colnames(tcga),
                                  colnames(gse42743)))

com = c("sample","time","status")
test_list1 = list(TCGA = tcga[,c(com,genelist1)],
                  GSE42743 = gse42743[,c(com,genelist1)],
                  GSE41613 = gse41613[,c(com,genelist1)])

result <- data.frame()
# TCGA as training
est_data <- test_list1$TCGA
# GEO as validation
val_data_list <- test_list1
pre_var <- colnames(est_data)[-c(1:3)]
est_dd <- est_data[, c('time', 'status', pre_var)]
val_dd_list <- lapply(val_data_list, function(x){x[, c('time', 'status', pre_var)]})

rf_nodesize <- 5
seed <- 123

## 1-1.RSF
set.seed(seed)
fit <- rfsrc(Surv(time,status)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize,  
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = seed)
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS  = predict(fit, newdata = x)$predicted)})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- 'RSF'
result <- rbind(result, cc)

## 1-2.RSF + CoxBostatust
set.seed(seed)
fit <- rfsrc(Surv(time, status)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize,  
             splitrule = 'logrank', 
             importance = T, 
             proximity = T, 
             forest = T, 
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
est_dd2 <- est_data[, c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('time', 'status', rid)]})
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd2[, 'time'], est_dd2[, 'status'], as.matrix(est_dd2[, -c(1, 2)]), 
                            trace=TRUE, start.penalty = 500, parallel = T)
cv.CoxBoost()
cv.res <- cv.CoxBoost(est_dd2[, 'time'], est_dd2[, 'status'], as.matrix(est_dd2[, -c(1, 2)]), 
                      maxstepno = 500, K = 10, type = "verweij",  penalty = pen$penalty)
fit <- CoxBoost(est_dd2[, 'time'], est_dd2[, 'status'], as.matrix(est_dd2[, -c(1, 2)]), 
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1],  newstatus = x[, 2], type = "lp")))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + CoxBostatust')
result <- rbind(result, cc)

## 1-3.RSF + Enet
set.seed(seed)
fit <- rfsrc(Surv(time, status)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, 
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
est_dd2 <- est_data[, c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('time', 'status', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$time, est_dd2$status))
for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + ', 'Enet', '[α=', alpha, ']')
  result <- rbind(result, cc)
}

## 1-4.RSF + GBM
set.seed(seed)
fit <- rfsrc(Surv(time, status)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize,  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
est_dd2 <- est_data[, c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('time', 'status', rid)]})
set.seed(seed)
fit <- gbm(formula = Surv(time, status)~., data = est_dd2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)

# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(time,status)~., data = est_dd2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'GBM')
result <- rbind(result, cc)

## 1-5.RSF + Lasso
set.seed(seed)
fit <- rfsrc(Surv(time, status)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, 
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
est_dd2 <- est_data[, c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('time', 'status', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$time, est_dd2$status))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, 
                family = "cox", alpha = 1,
                type.measure = "class")
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'Lasso')
result <- rbind(result, cc)

## 1-6.RSF + plsRcox
set.seed(seed)
fit <- rfsrc(Surv(time, status)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize,
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
est_dd2 <- est_data[, c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('time', 'status', rid)]})
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$time, status = est_dd2$status), nt = 10, verbose = FALSE)
fit <- plsRcox(est_dd2[, rid], time = est_dd2$time, event = est_dd2$status, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[, -c(1, 2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'plsRcox')
result <- rbind(result, cc)

## 1-7.RSF + Ridge
set.seed(seed)
fit <- rfsrc(Surv(time, status)~., data = est_dd, 
             ntree = 1000, nodesize = rf_nodesize, 
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")

rid <- rid$topvars
est_dd2 <- est_data[, c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('time', 'status', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$time, est_dd2$status))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold=10, 
                family = "cox", alpha = 0,
                type.measure = "class")
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'Ridge')
result <- rbind(result, cc)

## 1-8.RSF + StepCox
set.seed(seed)
fit <- rfsrc(Surv(time,status)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, 
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")

rid <- rid$topvars
est_dd2 <- est_data[, c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('time', 'status', rid)]})
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(time, status)~., est_dd2), direction = direction)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS=predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('RSF + ', 'StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

## 1-9.RSF + SuperPC
  set.seed(seed)
fit <- rfsrc(Surv(time,status)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, 
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
est_dd2 <- est_data[, c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('time', 'status', rid)]})
data <- list(x = t(est_dd2[, -c(1, 2)]), y = est_dd2$time,
             censoring.status = est_dd2$status, 
             featurenames = colnames(est_dd2)[-c(1, 2)])
set.seed(seed)
fit <- superpc.train(data = data, type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                     n.fold = 10, 
                     n.components = 3, 
                     min.features = 5, 
                     max.features = nrow(data$x), 
                     compute.fullcv = TRUE, 
                     compute.preval = TRUE)
rs <- lapply(val_dd_list2, function(w){
  test <- list(x = t(w[, -c(1, 2)]), y = w$time, censoring.status=w$status, featurenames = colnames(w)[-c(1, 2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1, ])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[, 1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'SuperPC')
result <- rbind(result, cc)

## 1-10.RSF + survival-SVM
set.seed(seed)
fit <- rfsrc(Surv(time, status)~., data = est_dd,
             ntree = 1000, nodesize = rf_nodesize, 
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rid <- var.select(object = fit, conservative = "high")
rid <- rid$topvars
est_dd2 <- est_data[, c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('time', 'status', rid)]})
fit = survivalsvm(Surv(time, status)~., data= est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS=as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('RSF + ', 'survival-SVM')
result <- rbind(result,cc)

##2. Enet
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$time, est_dd$status))

for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  rs <- lapply(val_dd_list,function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit,type = 'link', newx = as.matrix(x[,-c(1,2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(time, status) ~ rs, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('Enet', '[α=', alpha, ']')
  result <- rbind(result, cc)
}

##3. StepCox
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(time,status)~., est_dd), direction = direction)
  rs <- lapply(val_dd_list,function(x){cbind(x[, 1:2], RS = predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}


## 4-1.CoxBoost
set.seed(seed)
pen <- optimCoxBoostPenalty(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]), 
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBostatust(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]), maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBostatust(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rs <- lapply(val_dd_list, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1,2)], newtime = x[,1], newstatus = x[,2], type = "lp")))})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost')
result <- rbind(result, cc)

## 4-2.CoxBoost + Enet
set.seed(seed)
pen <- optimCoxBostatustPenalty(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBostatust(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBostatust(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)`!=0), "id"]
est_dd2 <- est_data[, c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('time', 'status', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$time, est_dd2$status))
for (alpha in seq(0.1, 0.9, 0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost', ' + Enet', '[α=', alpha, ']')
  result <- rbind(result, cc)
}

## 4-3.CoxBoost + GBM
set.seed(seed)
pen <- optimCoxBostatustPenalty(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBostatust(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K= 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBostatust(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)`!=0), "id"]
est_dd2 <- est_data[, c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('time', 'status', rid)]})
set.seed(seed)
fit <- gbm(formula = Surv(time,status)~., data = est_dd2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(time,status)~., data = est_dd2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list2,function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'GBM')
result <- rbind(result, cc)

## 4-4.CoxBoost + Lasso
set.seed(seed)
pen <- optimCoxBostatustPenalty(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBostatust(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBostatust(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                stepno=cv.res$optimal.step, penalty=pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('time', 'status', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$time, est_dd2$status))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, 
                family = "cox", alpha = 1,
                type.measure = "class")
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'Lasso')
result <- rbind(result, cc)

## 4-5.CoxBoost + plsRcox
set.seed(seed)
pen <- optimCoxBostatustPenalty(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBostatust(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBostatust(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                stepno=cv.res$optimal.step, penalty=pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('time', 'status', rid)]})
cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[,rid], time = est_dd2$OS.time, status = est_dd2$OS), nt = 10, verbose = FALSE)
fit <- plsRcox(est_dd2[, rid], time = est_dd2$OS.time, event = est_dd2$OS, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type="lp", newdata = x[, -c(1,2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'plsRcox')
result <- rbind(result, cc)


## 4-6.CoxBoost + Ridge
set.seed(seed)
pen <- optimCoxBostatustPenalty(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBostatust(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K=10, type="verweij", penalty = pen$penalty)
fit <- CoxBostatust(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('time', 'status', rid)]})
x1 <- as.matrix(est_dd2[, rid])
x2 <- as.matrix(Surv(est_dd2$time, est_dd2$status))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold=10, 
                family = "cox", alpha = 0,
                type.measure = "class")
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBostatust + ', 'Ridge')
result <- rbind(result, cc)

## 4-7.CoxBoost + StepCox
set.seed(seed)
pen <- optimCoxBostatustPenalty(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBostatust(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBostatust(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('time', 'status', rid)]})
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(time,status)~., est_dd2), direction = direction)
  rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('CoxBoost + ', 'StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

## 4-8.CoxBoost + SuperPC
set.seed(seed)
pen <- optimCoxBostatustPenalty(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBostatust(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[,-c(1,2)]),
                      maxstepno = 500, K= 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBostatust(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[,-c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[,c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('time', 'status', rid)]})
data <- list(x = t(est_dd2[, -c(1,2)]), y = est_dd2$time, censoring.status = est_dd2$status,
             featurenames = colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data, type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit, data, n.threshold = 20, #default
                     n.fold = 10,
                     n.components = 3,
                     min.features = 5,
                     max.features = nrow(data$x),
                     compute.fullcv = TRUE,
                     compute.preval =TRUE)
rs <- lapply(val_dd_list2, function(w){
  test <- list(x=t(w[, -c(1,2)]), y = w$time, censoring.status = w$status, featurenames = colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBostatust + ', 'SuperPC')
result <- rbind(result, cc)

## 4-9.CoxBoost + survival-SVM
set.seed(seed)
pen <- optimCoxBostatustPenalty(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBostatust(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBostatust(est_dd[, 'time'], est_dd[, 'status'], as.matrix(est_dd[,-c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rid <- as.data.frame(coef(fit))
rid$id <- rownames(rid)
rid <- rid[which(rid$`coef(fit)` != 0), "id"]
est_dd2 <- est_data[, c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('time', 'status', rid)]})
fit = survivalsvm(Surv(time, status)~., data = est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('CoxBoost + ', 'survival-SVM')
result <- rbind(result, cc)

##5. plsCox
set.seed(seed)
cv.plsRcox.res = cv.plsRcox(list(x = est_dd[,pre_var], time = est_dd$time, status = est_dd$status), nt = 10, verbstatuse = FALSE)
fit <- plsRcox(est_dd[,pre_var], time = est_dd$time, event = est_dd$status, nt = as.numeric(cv.plsRcox.res[5]))
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit,type = "lp", newdata = x[, -c(1, 2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time,status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('plsRcox')
result <- rbind(result, cc)

##6. SuperPC
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit,type = "lp", newdata = x[, -c(1, 2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time,status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('plsRcox')
result <- rbind(result, cc)

##7. GBM
set.seed(seed)
fit <- gbm(formula = Surv(time,status)~., data = est_dd, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(time, status)~., data = est_dd, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 8)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('GBM')
result <- rbind(result, cc)

##8. survivalSVM
fit = survivalsvm(Surv(time,status)~., data = est_dd, gamma.mu = 1)
rs <- lapply(val_dd_list, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('survival - SVM')
result <- rbind(result, cc)

##9. Ridge
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$time, est_dd$status))
set.seed(seed)
fit = glmnet(x1, x2, family = "cox", alpha = 0, lambda = NULL)
cvfit = cv.glmnet(x1, x2,
                  nfold = 10, 
                  family = "cox",
                  type.measure = "class"
)

rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = cvfit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time,status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Ridge')
result <- rbind(result, cc)

##10. Lasso
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$time, est_dd$status))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, 
                family = "cox", alpha = 1,
                type.measure = "class")
rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'response', newx = as.matrix(x[, -c(1,2)]), s = fit$lambda.min)))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso')
result <- rbind(result, cc)

## 10.1.Lasso + CoxBoost
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$time, est_dd$status))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, 
                family = "cox", alpha = 1,
                type.measure = "class")
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid <- rid[-1]
est_dd2 <- est_data[, c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('time', 'status', rid)]})
set.seed(seed)
pen <- optimCoxBostatustPenalty(est_dd2[, 'time'], est_dd2[, 'status'], as.matrix(est_dd2[, -c(1,2)]),
                            trace = TRUE, start.penalty = 500, parallel = T)
cv.res <- cv.CoxBostatust(est_dd2[, 'time'], est_dd2[, 'status'], as.matrix(est_dd2[, -c(1,2)]),
                      maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty)
fit <- CoxBostatust(est_dd2[, 'time'], est_dd2[, 'status'], as.matrix(est_dd2[, -c(1,2)]),
                stepno = cv.res$optimal.step, penalty = pen$penalty)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, newdata = x[,-c(1,2)], newtime = x[,1], newstatus = x[,2], type = "lp")))})
cc <- data.frame(Cindex=sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + CoxBoost')
result <- rbind(result, cc)

## 10.2.Lasso + GBM
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$time, est_dd$status))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, 
                family = "cox", alpha = 1,
                type.measure = "class")
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid <- rid[-1]
est_dd2 <- est_data[,c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('time', 'status', rid)]})
set.seed(seed)
fit <- gbm(formula = Surv(time,status)~., data = est_dd2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(time,status)~., data = est_dd2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 8)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])}))%>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'GBM')
result <- rbind(result, cc)

## 10.3.Lasso + plsRcox
est_dd[1:4,1:4]
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$time, est_dd$status))
set.seed(seed)

fit = cv.glmnet(x1, x2,
                nfold = 10, 
                family = "cox",alpha =1,
                type.measure = "class")

fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid <- rid[-1]
est_dd2 <- est_data[,c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('time', 'status', rid)]})
set.seed(seed)

cv.plsRcox.res = cv.plsRcox(list(x = est_dd2[, rid], time = est_dd2$time, status = est_dd2$status), nt = 10, verbstatuse = FALSE)

fit <- plsRcox(est_dd2[, rid], time = est_dd2$time, event = est_dd2$status, nt = as.numeric(cv.plsRcox.res[5]))

rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = "lp", newdata = x[,-c(1,2)])))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'plsRcox')
result <- rbind(result, cc)

## 10.4.Lasso + RSF
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$time, est_dd$status))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10,
                family = "cox", alpha = 1,
                type.measure = "class")
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid<-rid[-1]
est_dd2 <- est_data[,c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('time', 'status', rid)]})
set.seed(seed)
fit <- rfsrc(Surv(time,status)~., data = est_dd2,
             ntree = 1000, nodesize = rf_nodesize, ##该值建议多调整
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso', ' + RSF')
result <- rbind(result, cc)

## 10.5.Lasso + stepcox
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$time, est_dd$status))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, 
                family = "cox", alpha = 1,
                type.measure = "class")
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid <- rid[-1]
est_dd2 <- est_data[, c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('time', 'status', rid)]})
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(time,status)~., est_dd2), direction = direction)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = predict(fit, type = 'risk', newdata = x))})
  cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
    rownames_to_column('ID')
  cc$Model <- paste0('Lasso + ', 'StepCox', '[', direction, ']')
  result <- rbind(result, cc)
}

## 10.6.Lasso + superPC
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$time, est_dd$status))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, 
                family = "cox", alpha = 1,
                type.measure = "class")
fit$lambda.min
myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid <- rid[-1]
est_dd2 <- est_data[,c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[, c('time', 'status', rid)]})
data <- list(x = t(est_dd2[,-c(1,2)]), y = est_dd2$time, censoring.status = est_dd2$status,
             featurenames = colnames(est_dd2)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival', s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20, #default
                     n.fold = 10,
                     n.components = 3,
                     min.features = 5,
                     max.features = nrow(data$x),
                     compute.fullcv = TRUE,
                     compute.preval = TRUE)
rs <- lapply(val_dd_list2, function(w){
  test <- list(x = t(w[,-c(1,2)]), y = w$time, censoring.status = w$status, featurenames = colnames(w)[-c(1,2)])
  ff <- superpc.predict(fit, data, test, threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])], n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[, 1:2], RS = rr)
  return(rr2)
})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'SuperPC')
result <- rbind(result, cc)

## 10.7.Lasso + survival-SVM
x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$time, est_dd$status))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, 
                family = "cox", alpha = 1,
                type.measure = "class")
fit$lambda.min

myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
rid <- rid[-1]
est_dd2 <- est_data[,c('time', 'status', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('time', 'status', rid)]})
fit = survivalsvm(Surv(time,status)~., data = est_dd2, gamma.mu = 1)
rs <- lapply(val_dd_list2, function(x){cbind(x[,1:2], RS = as.numeric(predict(fit, x)$predicted))})
cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
  rownames_to_column('ID')
cc$Model <- paste0('Lasso + ', 'survival-SVM')
result <- rbind(result, cc)







