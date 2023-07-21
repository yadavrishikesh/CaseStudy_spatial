rm(list=ls())
library(INLA)
setwd(this.path::here())
 for(i in 1:5) {
   for (j in 1:3) {
   #i=1;j=1
train <- readr::read_csv(paste0("C:/Users/11322929/Dropbox (KAUST)/KAUST data challenge 2023/Content-paper-JABES/KaustSimulatedData/2a_", i, "_cross_val_3/cross_val_",j, "/train.csv" ))
test<- readr::read_csv(paste0("C:/Users/11322929/Dropbox (KAUST)/KAUST data challenge 2023/Content-paper-JABES/KaustSimulatedData/2a_", i, "_cross_val_3/cross_val_",j, "/test.csv" ))

data_train<- data.frame(train)
data_train$z<- data_train$data
data_train<- data_train[,-3]
test<- data.frame(test)
test_data<- test
test_data$z<- test$data
test_data<- test_data[,-3]


loc_all<- rbind(data_train[,1:2], test_data[,1:2])
test_z<- rep(NA, nrow(test_data))
whole_z<- c(data_train$z, test_z)

Data_all<- data.frame(cbind(x=loc_all[,1], y=loc_all[,2], z= whole_z))


LonCentre <- mean(range(Data_all$x))
LatCentre <- mean(range(Data_all$y))

Data_all$LonC <- Data_all$x - LonCentre
Data_all$LatC <- Data_all$y - LatCentre
#Define mesh
pts <- as.matrix(Data_all[,1:2])
mesh <- inla.mesh.2d(loc.domain = pts, max.edge = c(0.05, 0.1),
                     offset = c(0.1, 0.2) )


# par(mfrow=c(1,1))
# plot(mesh, asp = 1, main = "")
# points(pts, col = 3)


#Create SPDE
Arsenic.spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
A.meuse <- inla.spde.make.A(mesh=mesh, loc = pts)
s.index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = Arsenic.spde$n.spde)

#Create data structure
Arsenic.stack <- inla.stack(data  = list(z = Data_all$z),
                            A = list(A.meuse, 1),
                            effects = list(c(s.index, list(Intercept = 1)),
                                           list(x = Data_all$x, y = Data_all$y)),
                            tag = "Aresenic.pred.data")




#Create data structure for prediction
A.pred <- inla.spde.make.A(mesh = mesh, loc = pts)
Arsenic.stack.pred <- inla.stack(data = list(z = NA),
                                 A = list(A.pred, 1),
                                 effects = list(c(s.index, list(Intercept = 1)),
                                                list(x = Data_all$x, y = Data_all$y)),
                                 tag = "Aresenic.pred")

#Join stack
join.stack <- inla.stack(Arsenic.stack, Arsenic.stack.pred)


#Fit model
run_time <- system.time({
form <- z ~ -1 + Intercept + x + y + f(spatial.field, model = spde)
m1 <- inla(form, data = inla.stack.data(join.stack, spde = Arsenic.spde),
           family = "gaussian", num.threads=1,
           control.predictor = list(A = inla.stack.A(join.stack), compute = TRUE),
           control.compute = list(cpo = TRUE, dic = TRUE), verbose = TRUE)
})
#Summary of results
#summary(m1)
#Get predicted data on grid
index.pred <- inla.stack.index(join.stack, "Aresenic.pred")$data
Data_all$post_mean <- m1$summary.linear.predictor[index.pred, "mean"]
Data_all$post_sd <- sqrt(m1$summary.linear.predictor[index.pred, "sd"]^2 + 1/m1$summary.hyperpar[1,"0.5quant"]) 

# Data_all$lCI<-  m1$summary.fitted.values[index.pred, "0.025quant"] 
# Data_all$uCI<-  m1$summary.fitted.values[index.pred, "0.975quant"] 
pred_mean_sd<- Data_all

save(pred_mean_sd, run_time, file=paste0("CV_results_","dataset_n0=",i, "_cross_val_no=",j,".RData"))

   }
   print(i);print(j)
}

