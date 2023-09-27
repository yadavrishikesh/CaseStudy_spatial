rm(list = ls())
dataset.no<- 1
cross_val_no<- 1

setwd(this.path::here())
library(spNNGP)
source("conjugate-NNGP_function.R")
train <- readr::read_csv(paste0("KaustSimulatedData/2a_", i, "_cross_val_3/cross_val_",j, "/train.csv"))
test <- readr::read_csv(paste0("KaustSimulatedData/2a_", i, "_cross_val_3/cross_val_",j, "/test.csv"))
data_train<- data.frame(train)
data_train$z<- data_train$data
data_train<- data_train[,-3]
head(data_train)
test<- data.frame(test)
test_data<- test
test_data$z<- test$data
test_data<- test_data[,-3]
head(test_data)
library(spNNGP)
fit_conj_model<- conjNNGP(data_train=data_train, 
                          S=as.matrix(data_train[,1:2]), 
                          S.pred=as.matrix(test_data[,1:2]), 
                          X.0 = as.matrix(cbind(1,test_data[,-3])),
                          n.neighbors=15,
                          k.fold=5,
                          init.par=c(1,1,1,1),
                          n.omp.threads=1,
                          score.rule="rmspe",
                          n.samples=200,
                          discrete.level = 3)
save.image(file = paste0("CV_results_","dataset_n0=",dataset.no, "_cross_val_no=",cross_val_no,".RData"))


