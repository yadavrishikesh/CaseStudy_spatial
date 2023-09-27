# Clear the current environment
rm(list = ls())

# Load the 'INLA' package for spatial modeling
library(INLA)
library(rSPDE)
library(inlabru)
# Set the working directory using the 'here' package (assuming it's already installed)
setwd(this.path::here())

# Loop through different datasets (i from 1 to 5)
for (i in 1:5) {
  # Loop through different cross-validation folds within each dataset (j from 1 to 3)
  for (j in 1:3) {
    #i=1;j=1
    # Read the train and test datasets for the current fold
    train <- readr::read_csv(paste0("C:/Users/11322929/Dropbox (KAUST)/KAUST data challenge 2023/Content-paper-JABES/KaustSimulatedData/2a_", i, "_cross_val_3/cross_val_",j, "/train.csv"))
    test <- readr::read_csv(paste0("C:/Users/11322929/Dropbox (KAUST)/KAUST data challenge 2023/Content-paper-JABES/KaustSimulatedData/2a_", i, "_cross_val_3/cross_val_",j, "/test.csv"))
    
    # Convert the train and test data to data frames
    data_train <- data.frame(train)
    
    # Define the mesh using 'inla.mesh.2d' based on the spatial locations
    pts <- as.matrix(data_train[,1:2])
    mesh <- inla.mesh.2d(loc.domain = pts, max.edge = c(0.05, 0.1),
                         offset = c(0.1, 0.2))
    
    rspde_model <- rspde.matern(
      mesh = mesh,
      nu.upper.bound = 2,
      parameterization = "spde"
    )
    
    
    df_rspde <- data.frame(coord1 = pts[,1],
                           coord2 = pts[,2],
                           z = as.vector(data_train$data))
    
    coordinates(df_rspde) <- c("coord1", "coord2")
    
    cmp <- z ~ -1 + field(coordinates, model = rspde_model)
    
    run_time <- system.time({
    rspde_bru_fit <-
      bru(cmp,
          data=df_rspde,
          options=list(
            family = "gaussian")
      )
    })
    
    # summary(rspde_bru_fit)
    # ##### prediction by kriging 
    # result_fit <- rspde.result(rspde_bru_fit, "field", rspde_model)
    # summary(result_fit)

    ### prediction at unknown spatil locations 
    test<- data.frame(test)
    pred_coords <- data.frame(x1 = test$x, x2 = test$y)
    coordinates(pred_coords) <- c("x1", "x2")
    field_pred <- predict(rspde_bru_fit, pred_coords, ~field)
    
  ### collecting all data and extracting the posterior predictive mean and standard errors
    Data_all <- data.frame(cbind(x=rbind(data_train[,1:2], test_data[,1:2])[,1], 
                                 y=rbind(data_train[,1:2], test_data[,1:2])[,2], 
                                 z=c(data_train$z, rep(NA, nrow(test)))))
    #Get predicted data on grid
    Data_all$post_mean <-  ifelse(is.na(Data_all$z), field_pred$mean, Data_all$z) #field_pred rspde_bru_fit$summary.linear.predictor[, "mean"]
    Data_all$post_sd <- ifelse(is.na(Data_all$z), sqrt(field_pred$sd^2 + 1/rspde_bru_fit$summary.hyperpar[1,"0.5quant"]), 0) 

    pred_mean_sd<- Data_all
    
    save(pred_mean_sd, run_time, file=paste0("CV_results_","dataset_n0=",i, "_cross_val_no=",j,".RData"))
    
  }
  #print(i);print(j)
}
