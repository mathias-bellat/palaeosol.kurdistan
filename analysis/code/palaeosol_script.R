####################################################################
# This script is the digital palaeosol mapping for Kurdistan area  #                       
# from 2,500 to 300 BCE                                            #                                                   
# Author: Mathias Bellat                                           #
# Affiliation : Tubingen University                                #
# Creation date : 18/08/2026                                       #
# E-mail: mathias.archaeology@gmail.com                            #
####################################################################


# 0 Environment setup ##########################################################

# 0.1 Prepare environment ======================================================

# Folder check
getwd()

# Clean up workspace
rm(list = setdiff(ls(), "save_dir"))

# 0.2 Install packages =========================================================

install.packages("pacman")        
#Install and load the "pacman" package (allow easier download of packages)
library(pacman)
pacman::p_load(tidyr, sf, mapview, terra, clhs, mapview, corrplot, viridis, usdm, ggplot2, caret, cli, rSDI, 
               MLmetrics, yardstick, pROC, gridExtra, ranger, randomForest, MASS, Matrix, blockCV,
               reshape2, dplyr, patchwork)

# 0.3 Show session infos =======================================================

sessionInfo()

# 01 cLHS ######################################################################

# 01.1 Import layers and AOI ===================================================
DEM <- rast("./analysis/data/raw_data/DEM_smooth.tif")
Soil <- rast("./analysis/data/raw_data/Soil.tif")
Geol <- rast("./analysis/data/raw_data/Geology.tif")
Geomorpho <- rast("./analysis/data/raw_data/Geomorphology.tif")
Precipitation <- rast("./analysis/data/raw_data/WorldClim_BIO12.tif")
Temperature <- rast("./analysis/data/raw_data/WorldClim_BIO01.tif")

# 1.2 Clean and prepare the covariates =========================================

# Transform the Soil
Soil <- as.factor(Soil)
levels(Soil) <- data.frame(ID = c(2:7,9:13),
                           classe = c("LP", "RG", "CCcm", "FL", "GY", "KZ", "VR", "CM", "CCrg", "CCvrfl", "CMvr"))

# Create a list of the different variables
stack <- c(DEM, Soil, as.factor(Geol), as.factor(Geomorpho), Precipitation, Temperature)
names(stack) <- c("DEM", "Soil", "Geology", "Geomorphology", "Precipitation", "Temperature")

# 1.3 Save and export the predictors ===========================================
writeRaster(stack ,filename = "./analysis/data/raw_data/predictors_clhs.tif", overwrite = TRUE)
plot(stack)

rm(list = setdiff(ls(), "save_dir"))

# 2 Create the Conditioned Latin Hypercube Sampling  ###########################
# 2.1 Import and prepare the files =============================================
predictors <- rast("./analysis/data/raw_data/predictors_clhs.tif")
crs(predictors) <- "EPSG:32638"

# Convert the file into Data Frame
PredForMap  <- as.data.frame(predictors, xy=TRUE) 
stack_frame <- PredForMap[complete.cases(PredForMap),] # Remove NA values
stack_frame <- na.omit(stack_frame) 

# 2.2 Conditioned Latin Hypercube sampling parameters ==========================

# Select the predictors for the Conditioned Latin Hypercube
preds <- c(names(stack_frame[,-c(1:2)])) 

# Set the size of the sampling 
c = 5000

# Set a seed
set.seed(1070) 

# Run the sampling
res <- clhs(stack_frame[, preds], size = c, iter = 20000, simple = FALSE, progress = TRUE)

pdf("./analysis/data/derived_data/model/CLHS_iter.pdf",    # File name
    width = 5, height = 6,  # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk")   # Color model 


plot(res)
dev.off()
# 2.3 Export the results =======================================================

CLHS_sampled_res <- stack_frame[res$index_samples, ] #fit the results with the line in the table
write.table(CLHS_sampled_res, "./analysis/data/derived_data/CLHS.csv", col.names = TRUE, sep = ";", row.names = FALSE, 
            fileEncoding = "UTF-8")

CLHS_sampled_sf <- st_as_sf(CLHS_sampled_res, coords = c("x", "y"), crs = 32638)
mapview(CLHS_sampled_sf)
rm(list = ls())

# 03 Import data sets for model ################################################

# 03.1 Import soils infos ======================================================
CLHS <- read.csv("./analysis/data/derived_data/CLHS.csv", sep=";")
DEM <- rast("./analysis/data/raw_data/DEM_smooth.tif")
Soil <- rast("./analysis/data/raw_data/Soil.tif")
covariates <- rast("./analysis/data/raw_data/CHELSA_present.tif")
covariates <- c(covariates, DEM)
names(covariates) <- c(names(covariates[[1:(nlyr(covariates)-1)]]), "DEM")

Soil <- as.factor(Soil)
levels(Soil) <- data.frame(ID = c(2:7,9:13),
                           classe = c("LP", "RG", "CCcm", "FL", "GY", "KZ", "VR", "CM", "CCrg", "CCvrfl", "CMvr"))

CLHS_soil <- terra::extract(Soil, CLHS[1:2], method = "simple")
CLHS_df <- cbind(CLHS[,1:2], CLHS_soil[,2])
CLHS_df[,3] <- as.factor(CLHS_df[,3]) 
plot(Soil)
plot(CLHS_df[,3])
CLHS_df[,3] <- droplevels(CLHS_df[,3])

# Extract the values of each band for the sampling location
df_cov <- raster::extract(covariates, CLHS_df[,1:2] , method = "bilinear")
df_cov <- as.data.frame(df_cov)
soil_infos <- cbind(CLHS_df, df_cov[,-1])
colnames(soil_infos) <- c("x", "y", "soil", colnames(soil_infos[4:22]), "DEM")
df_cov <- soil_infos[,-3]

# 03.2 Calculate euclidean distances ===========================================

# Extract raster coordinate
coords <- crds(covariates) 
center <- colMeans(coords, na.rm = TRUE)
corners <- matrix(c(
  ext(covariates)[1], ext(covariates)[3],
  ext(covariates)[1], ext(covariates)[3],
  ext(covariates)[2], ext(covariates)[4],
  ext(covariates)[2], ext(covariates)[4]), 
  ncol = 2, byrow = TRUE)

# Produce points Euclidian center distance
pts_coords <- data.frame(x = soil_infos$x, y = soil_infos$y)
dist_center <- euclidean(pts_coords[,1], pts_coords[,2], center[1], center[2])
dist_corner_NW <- euclidean(pts_coords[,1], pts_coords[,2], corners[1,1], corners[1,2])
dist_corner_SW <- euclidean(pts_coords[,1], pts_coords[,2], corners[2,1], corners[2,2])
dist_corner_SE <- euclidean(pts_coords[,1], pts_coords[,2], corners[3,1], corners[3,2])
dist_corner_NE <- euclidean(pts_coords[,1], pts_coords[,2], corners[4,1], corners[4,2])

distances <- data.frame(
  x = pts_coords[,1],
  y = pts_coords[,2],
  EDC = dist_center,
  EDC1 = dist_corner_NW,
  EDC2 = dist_corner_SW,
  EDC3 = dist_corner_SE,
  EDC4 = dist_corner_NE)

df_cov <- cbind(distances, df_cov[3:ncol(df_cov)])

# For the covariates maps

ED_cov <- as.data.frame(coords)
dist_center <- euclidean(ED_cov[,1], ED_cov[,2], center[1], center[2])
dist_corner_NW <- euclidean(ED_cov[,1], ED_cov[,2], corners[1,1], corners[1,2])
dist_corner_SW <- euclidean(ED_cov[,1], ED_cov[,2], corners[2,1], corners[2,2])
dist_corner_SE <- euclidean(ED_cov[,1], ED_cov[,2], corners[3,1], corners[3,2])
dist_corner_NE <- euclidean(ED_cov[,1], ED_cov[,2], corners[4,1], corners[4,2])

ED_cov <- data.frame(
  x = ED_cov[,1],
  y = ED_cov[,2],
  EDC = dist_center,
  EDC1 = dist_corner_NW,
  EDC2 = dist_corner_SW,
  EDC3 = dist_corner_SE,
  EDC4 = dist_corner_NE)

ED_rast <- rast(ED_cov, type = "xyz")
ED_rast$x <- ED_cov$x
ED_rast$y <- ED_cov$y

ED_rast <- resample(ED_rast, covariates, method ="bilinear")
covariates <- c(ED_rast, covariates)

covariates <- covariates[[colnames(df_cov)]]

# 03.3 Export and save data ====================================================
writeRaster(covariates, "./analysis/data/raw_data/predictors.tif", overwrite = TRUE)
write.csv(df_cov, "./analysis/data/raw_data/df_cov.csv")
save(df_cov, soil_infos, file = "./analysis/data/derived_data/save/Pre_process_model.RData")
rm(list = ls())

# 04 Check the data ############################################################

# 04.1 Import the data =========================================================
load(file = "./analysis/data/derived_data/save/Pre_process_model.RData")
covariates <- rast("./analysis/data/raw_data/predictors.tif")

# 04.2 Plot and export the correlation matrix ==================================

pdf("./analysis/data/derived_data/model/Correlation_matrix.pdf",    # File name
    width = 10, height = 10,  # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk")   # Color model 


# Correlation of the data
corrplot(cor(df_cov),  method = "color", col = viridis(200), 
         type = "upper", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, # Text label color and rotation
         number.cex = 0.7, # Size of the text labels
         cl.cex = 0.7) # Size of the color legend text) # Color legend limits

dev.off()

# 04.3 Select with VIF correlation =============================================

vif <- vifcor(df_cov, th=0.80)
vif_df <- as.data.frame(vif@results)

vif_plot <- ggplot(vif_df, aes(x = reorder(Variables, VIF), y = VIF)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "VIF Values for soil", x = "Variables", y = "VIF") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

vif_plot

ggsave("./analysis/data/derived_data/model/VIF_soil.png", vif_plot, width = 12, height = 8)
ggsave("./analysis/data/derived_data/model/VIF_soil.pdf", vif_plot, width = 12, height = 8)

write.table(vif_df, "./analysis/data/derived_data/model/vif_results.txt")

# 04.4 Select with RFE correlation =============================================

y <- soil_infos$soil
set.seed(1070)
ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   number = 10,
                   repeats = 3,
                   verbose = TRUE)


start_time <- Sys.time()

rfe_rf <- rfe(x = df_cov, y = y,
              sizes = c(3, 6, 9, 12, 15, 18, 20),
              rfeControl = ctrl)

end_time <- Sys.time()
cat("Time spend for RFE :", round(difftime(end_time, start_time, units="mins"), 2), "minutes\n")


print(rfe_rf)
predictors(rfe_rf)
plot(rfe_rf)

imp <- varImp(rfe_rf)
imp_df <- data.frame(Variable = rownames(imp),
                     Importance = imp$Overall)
imp_selected <- imp_df[imp_df$Variable %in% predictors(rfe_rf), ]

# Plot ggplot
rfe_plot <- ggplot(imp_selected, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Importance for RFE",
       x = "Variables",
       y = "Importance (Overall)") 

rfe_plot

ggsave("./analysis/data/derived_data/model/RFE_soil.png", rfe_plot, width = 12, height = 8)
ggsave("./analysis/data/derived_data/model/RFE_soil.pdf", rfe_plot, width = 12, height = 8)

write.table(imp_selected, "./analysis/data/derived_data/model/rfe_results.txt")
imp_selected$Variable <- sort(imp_selected$Variable)
covariates <- covariates[[imp_selected$Variable]]
writeRaster(covariates, "./analysis/data/derived_data/predictors_selected.tif", overwrite = TRUE)
save(rfe_rf, vif, file = "./analysis/data/derived_data/save/Features_selection.RData")

# 05 Pre models ################################################################

# 05.1 Pre-models preparation ==================================================

# Select the less correlated variables if you want to use RFE or VIF
SoilCovMLCon <- df_cov[, colnames(df_cov) %in% imp_selected$Variable] 

# Without VIF or RFE
#SoilCovMLCon <- soil_infos[, -c(1:3)] 

soil_infos[,c(3)] <- droplevels(soil_infos[,c(3)])
SoilCovMLCon <- cbind(SoilCovMLCon, soil_infos[,c(3)])
var_num <- length(SoilCovMLCon) - 1
colnames(SoilCovMLCon) <- c(colnames(SoilCovMLCon[,1:var_num]),  "soil")

preProcValues <- preProcess(SoilCovMLCon[,1:var_num], method = c("range"))
SoilCovMLConTrans <- predict(preProcValues, SoilCovMLCon[,1:var_num])
SoilCovMLConTrans <- cbind(SoilCovMLConTrans, soil_infos[,c(3)])
colnames(SoilCovMLConTrans) <- c(colnames(SoilCovMLCon[,1:var_num]),  "soil")

FormulaMLCon = as.formula(paste(names(SoilCovMLCon)[var_num+1]," ~ ",paste(names(SoilCovMLCon)[1:var_num],collapse="+")))

# 05.2 Develop pre-models ======================================================

# Define traincontrol
TrainControl <- trainControl(method="repeatedcv", 10, 3, allowParallel = TRUE, savePredictions=TRUE, verboseIter = TRUE)
seed=1070

cli_progress_bar(
  format = "  Running pre-models {cli::pb_bar} {cli::pb_percent} [{cli::pb_current}/{cli::pb_total}] | \ ETA: {cli::pb_eta} - Time elapsed: {cli::pb_elapsed_clock}",
  total = 5, clear = FALSE)

# Train different ML algorithms
#rpart (CART)
FitRpartCon <- train(FormulaMLCon, data=SoilCovMLCon, 
                     method="rpart", metric="Accuracy", trControl=TrainControl)
cli_progress_update()
print("CART done")

#Knn
FitKnnCon <- train(FormulaMLCon, data=SoilCovMLConTrans, 
                   method="knn", metric="Accuracy", trControl=TrainControl)
cli_progress_update()
print("Knn done")

# SVMr
FitSvrCon  <- train(FormulaMLCon, data=SoilCovMLConTrans, 
                    method="svmRadial", metric="Accuracy", trControl=TrainControl)
cli_progress_update()
print("SVM done")

# RF
FitRaFCon  <- train(FormulaMLCon, data=SoilCovMLCon, 
                    method="ranger", metric="Accuracy", trControl=TrainControl, importance = "impurity")
cli_progress_update()
print("RF done")

# ANN
FitANNCon  <- train(FormulaMLCon, data=SoilCovMLConTrans, 
                    method="nnet", metric="Accuracy", trControl=TrainControl, linout = FALSE)
cli_progress_update()
print("ANN done")
cli_progress_done()


# 05.3 Combine models statistics ===============================================
# Look at the primary results of ML

ModelConList <- list(CART=FitRpartCon, Knn=FitKnnCon,SVM=FitSvrCon, RF=FitRaFCon, ANN=FitANNCon)
ResultsModelCon<- resamples(ModelConList)
SummaryModelCon <- summary(ResultsModelCon)
ScalesMolel <- list(x=list(relation="free"), y=list(relation="free"))
BwplotModelCon <- bwplot(ResultsModelCon, scales=ScalesMolel, main = "Comparative models")
plot(BwplotModelCon)

pdf("./analysis/data/derived_data/model/Models_results.pdf",    # File name
    width = 10, height = 5,  # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk")   # Color model 
BwplotModelCon

dev.off()

# Calculate Error indices
Error1Con <- NaN*seq(length(FormulaMLCon))
for(j in 1:(2 * length(ModelConList))) { 
  Error1Con[j] <- mean(SummaryModelCon$values[[j]])
}
Error <- as.data.frame(Error1Con[c(1,3,5,7,9)])
Error <- cbind(Error, as.data.frame(Error1Con[c(2,4,6,8,10)]))
colnames(Error) <- c("Accuracy", "Kappa")
rownames(Error) <- c("CART", "Knn", "SVMr", "RF", "ANN")

write.table(Error, "./analysis/data/derived_data/model/Models_errors.txt")
Error

save(SoilCovMLCon, var_num, imp_selected, FormulaMLCon, file = "./analysis/data/derived_data/save/First_models.RData")
rm(list = ls())

# 06 Gaussian noice ############################################################

# 06.1 Noice values ============================================================

v.ed <- 15.1 # The variability for Euclidean distance is set at the the distance to a cells from the center (originally 30m cells)
v.prec <- 120 # variability prec CHELSA = 0.71 RMSE 120
v.tmax <- 3  # variability tasmax CHELSA = 0.37 RMSE 3
v.tmin <- 4.2  # variability tasmin CHELSA = 0.37 RMSE 4.2
v.tmean <- mean(c(v.tmax ,v.tmin))
v.prec.month <- v.prec/12
v.prec.quarter <- v.prec/4
v.dem <- 4 # Vertical error of the GLO Copernicus measurement

gaus.noice <- as.data.frame(t(c(v.ed, v.ed, v.ed, v.ed, v.ed, v.ed, v.ed, v.tmean, v.tmean, v.tmean, v.tmean, v.tmax, v.tmin, v.tmean, v.tmean, v.tmean, v.tmean, v.tmean, 
                                v.prec, v.prec.month, v.prec.month, v.prec.month, v.prec.quarter, v.prec.quarter, 
                                v.prec.quarter, v.prec.quarter, v.dem)))


bio_names <- sprintf("bio%02d", 1:19)
names(gaus.noice) <- c("x", "y", "EDC", "EDC1", "EDC2", "EDC3", "EDC4", bio_names, "DEM")

n_iter = 1000
set.seed(1070)
seed.list <- sample(1000:9999, n_iter, replace=FALSE)

save(gaus.noice, n_iter, seed.list, file = "./analysis/data/derived_data/save/Gaussian_noice.RData")
rm(list = ls())

# 07 Random forest final model #################################################

# 07.1 Prepare the BlockCV =====================================================
load("./analysis/data/derived_data/save/First_models.RData")
covariates <- rast("./analysis/data/raw_data/predictors.tif")

start_time <- Sys.time()

sac1 <- cv_spatial_autocor(r = covariates, num_sample = 5000)

ggsave("./analysis/data/derived_data/model/blockCV/BlockCV_autocorrelation_range.png", sac1$plots$barchart, width = 12, height = 8)
ggsave("./analysis/data/derived_data/model/blockCV/BlockCV_autocorrelation_range.pdf", sac1$plots$barchart, width = 12, height = 8)

ggsave("./analysis/data/derived_data/model/blockCV/BlockCV_autocorrelation_spatial_block.png", sac1$plots$map_plot, width = 12, height = 8)
ggsave("./analysis/data/derived_data/model/blockCV/BlockCV_autocorrelation_spatial_block.pdf", sac1$plots$map_plot, width = 12, height = 8)


points <- sf::st_as_sf(SoilCovMLCon, coords = c("x", "y"), crs = crs(covariates))

points$soil <- as.numeric(points$soil)
points <- points[,c(14,15)]


sac2 <- cv_spatial_autocor(x = points, 
                           column =  "soil")

png("./analysis/data/derived_data/model/blockCV/BlockCV_autocorrelation_spatial_block.png", width = 2000, height = 1500, res = 300)
pdf("./analysis/data/derived_data/model/blockCV/BlockCV_autocorrelation_spatial_block.pdf", width = 8, height = 6)

plot(sac2$variograms[[1]])

dev.off()

# For a visualisation of the blockCV size (11,650 m here)
cv_block_size(x = points, 
              column =  "soil",
              r = covariates,
              min_size = 5000,
              max_size = 50000)


# 07.2 Set the BlockCV =========================================================

#Split the data
set.seed(1070)
var_num <- length(SoilCovMLCon) - 1
split <- createDataPartition(SoilCovMLCon[,var_num+1], p = 0.8, list = FALSE, times = 1)
Train_data <- SoilCovMLCon[ split,]
Test_data  <- SoilCovMLCon[-split,]

points <- sf::st_as_sf(Train_data, coords = c("x", "y"), crs = crs(covariates))

points$soil <- as.numeric(points$soil)
points <- points[,c(14,15)]

# Spatial blocking
n_repeats <- 3
seed.blockCV <- c(1070, 2019, 8215)
all_folds <- list()

for (i in 1:n_repeats) {
  
  pdf(paste0("./analysis/data/derived_data/model/blockCV/Block_CV_folds_raw_iter_",i,".pdf"), width = 8, height = 6)
  
  set.seed(seed.blockCV[i])
  sb <- cv_spatial(
    x = points, 
    column = "soil",
    r = covariates,
    size = 11000,
    k = 10,
    selection = "random",
    iteration = 100
  )
  dev.off()
  
  all_folds[[i]] <- sb$folds_list
  
  pdf(paste0("./analysis/data/derived_data/model/blockCV/Block_CV_folds_visualisation_iter_",i,".pdf"), width = 8, height = 6)
  cv_plot(cv = sb, 
          x = points)
  dev.off()
}

end_time <- Sys.time()
cat("Time spend for BlocCV :", round(difftime(end_time, start_time, units="mins"), 2), "minutes\n")


# 07.3 Prepare the model tuning ================================================

index <- list()
indexOut <- list()

counter <- 1

for (rep in 1:n_repeats) {
  for (fold in 1:length(all_folds[[rep]])) {
    
    index[[counter]]    <- all_folds[[rep]][[fold]][[1]]
    indexOut[[counter]] <- all_folds[[rep]][[fold]][[2]]
    
    counter <- counter + 1
  }
}

# Training control and hyperparameters
TrainControl <- trainControl(
  method = "cv",
  index = index,
  indexOut = indexOut,
  search = "random",
  allowParallel = TRUE,
  savePredictions = TRUE,
  verboseIter = TRUE,
  classProbs = TRUE,
  summaryFunction = multiClassSummary
)

num_trees_values <- c(500, 750, 1000)

grid <- expand.grid(
  mtry = seq(1, ncol(Train_data) - 1),
  splitrule = "extratrees",
  min.node.size = 1
)

# 07.4 Run the models ==========================================================

results <- list()
FinalRFCon <- list()
cli_progress_bar(
  format = "  Running final models {cli::pb_bar} {cli::pb_percent} [{cli::pb_current}/{cli::pb_total}] | \ ETA: {cli::pb_eta} - Time elapsed: {cli::pb_elapsed_clock}",
  total = length(num_trees_values), clear = FALSE)

# Test for every num.trees
for (num_trees in num_trees_values) {
  cat("Testing num.trees =", num_trees, "\n")
  set.seed(1070)
  model  <- train(FormulaMLCon, data = Train_data, 
                  method="ranger", metric="Accuracy", trControl = TrainControl, 
                  importance = "impurity", tuneGrid = grid, num.trees = num_trees)
  
  FinalRFCon[[as.character(num_trees)]] <- model
  model$results$num.trees <- num_trees
  results[[as.character(num_trees)]] <- model$results
  cli_progress_update()
}
cli_progress_done()

final_results <- do.call(rbind, results)


# 07.5 Plot model tuning parameters ============================================

head(final_results[order(-final_results$Accuracy), ], 5)

gg1 <- ggplot(final_results, aes(x = mtry, y = Accuracy, color = factor(num.trees))) +
  geom_line() +
  geom_point() +
  labs(
    title = "Model performance regarding mtry and num.trees",
    x = "Number of variables used (mtry)",
    y = "Accuracy",
    color = "num.trees"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),  
    legend.position = "bottom"          
  )

ggsave("./analysis/data/derived_data/model/Tuning_parameters.png", gg1, width = 12, height = 8)
ggsave("./analysis/data/derived_data/model/Tuning_parameters.pdf", gg1, width = 12, height = 8)
gg1

RF_final_model <- FinalRFCon[[1]]

save(all_folds, FinalRFCon, RF_final_model, file = "./analysis/data/derived_data/save/Models.RData")
rm(list=ls())

# 07.6 Noice generation ========================================================
load("./analysis/data/derived_data/save/First_models.RData")
load("./analysis/data/derived_data/save/Models.RData")
covariates <- rast("./analysis/data/raw_data/predictors.tif")
load("./analysis/data/derived_data/save/Gaussian_noice.RData")

# Import data and create Sigma correlation matrix
set.seed(1070)
var_num <- length(SoilCovMLCon) - 1
split <- createDataPartition(SoilCovMLCon[,var_num+1], p = 0.8, list = FALSE, times = 1)
Train_data <- SoilCovMLCon[ split,]
Test_data  <- SoilCovMLCon[-split,]

R <- cor(Train_data[,1:var_num], use = "complete.obs")
sd_vec <- as.numeric(gaus.noice[1, colnames(R)])
names(sd_vec) <- colnames(R)

n_obs <- nrow(Test_data)
num_cols <- names(Test_data)[sapply(Test_data, is.numeric)]
classes <- levels(Test_data$soil)

y_true <- as.character(Test_data$soil)

Sigma_base <- diag(sd_vec) %*% R %*% diag(sd_vec)
Sigma_base <- as.matrix(nearPD(Sigma_base)$mat)

# Create the functions for evaluation and Monte-Carlo

simulate_mc <- function(alpha) {
  
  Sigma <- alpha * Sigma_base
  
  results <- array(NA, dim = c(n_obs, length(classes), n_iter))
  
  cli_progress_bar(
    format = "Monte Carlo {alpha}  {cli::pb_bar} {cli::pb_percent} [{cli::pb_current}/{cli::pb_total}] | \ ETA: {cli::pb_eta} - Time elapsed: {cli::pb_elapsed_clock}",
    total = n_iter, clear = FALSE)
  
  for (i in 1:n_iter) {
    
    set.seed(seed.list[i])
    noise <- mvrnorm(n_obs, mu = rep(0, length(num_cols)), Sigma = Sigma)
    colnames(noise) <- num_cols
    
    df_noisy <- Test_data
    df_noisy[, num_cols] <- df_noisy[, num_cols] + noise
    
    pred <- predict(
      FinalRFCon[[1]]$finalModel,
      data = df_noisy,
      type = "response"
    )$predictions
    
    pred <- pred[, classes, drop = FALSE]
    
    results[,,i] <- pred
    cli_progress_update()
  }
  cli_progress_done()
  colnames(results) <- classes
  return(results)
}

evaluate_mc <- function(results) {
  
  mean_pred <- apply(results, c(1,2), mean)
  pred_class <- colnames(mean_pred)[max.col(mean_pred)]
  
  acc <- mean(pred_class == y_true)
  
  pred_matrix <- apply(results, 3, function(mat) {
    max.col(mat)
  })
  
  stability <- mean(apply(pred_matrix, 1, function(x) {
    length(unique(x)) == 1
  }))
  
  return(list(
    accuracy = acc,
    stability = stability
  ))
}

# Run Monte Carlo with 5 different alpha values

alpha_values <- c(0.00005, 0.0001, 0.0002, 0.0003, 1)

simulation_results <- list()

results_alpha <- data.frame(
  alpha = alpha_values,
  accuracy = NA,
  stability = NA
)

for (j in seq_along(alpha_values)) {
  
  cat("Running alpha =", alpha_values[j], "\n")
  
  res <- simulate_mc(alpha_values[j])
  simulation_results[[j]] <- res
  eval <- evaluate_mc(res)
  
  results_alpha$accuracy[j] <- eval$accuracy
  results_alpha$stability[j] <- eval$stability
}

save(simulation_results, Sigma_base, results_alpha, file = "./analysis/data/derived_data/save/Monte_Carlo.RData")
# 07.7 Monte Carlo evaluation ==================================================

# The 1 value for delta is to low so we removed it from the plots
results_alpha

gg1 <- ggplot(results_alpha[-5,], aes(x = alpha, y = accuracy)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Accuracy vs Noise Level (alpha)",
       x = "Alpha",
       y = "Accuracy")

gg1
ggsave("./analysis/data/derived_data/model/Monte_Carlo_accuracy.pdf", gg1, width = 12, height = 8)
ggsave("./analysis/data/derived_data/model/Monte_Carlo_accuracy.png", gg1, width = 12, height = 8)

gg1 <-ggplot(results_alpha[-5,], aes(x = alpha, y = stability)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Prediction Stability vs Noise",
       x = "Alpha",
       y = "Stability")

gg1
ggsave("./analysis/data/derived_data/model/Monte_Carlo_stability.pdf", gg1, width = 12, height = 8)
ggsave("./analysis/data/derived_data/model/Monte_Carlo_stability.png", gg1, width = 12, height = 8)


# The selected values of delta is 1e-04
results_selected <- simulation_results[[2]]

# Create the "raw" prediction set
raw_pred <- predict(FinalRFCon[[1]]$finalModel, data = Test_data, type = "response")$predictions
raw_pred <- raw_pred[, classes, drop = FALSE]
pred_class <- colnames(raw_pred)[max.col(raw_pred)]

# Create the delta for accuracy
acc_ref <- mean(pred_class == y_true)
delta_acc_mean <- results_alpha[2,2] - acc_ref

# Histrogram of accuracy derivation

pred_matrix <- apply(results_selected, 3, function(mat) {
  colnames(mat)[max.col(mat)]
})

conf_list <- lapply(1:ncol(pred_matrix), function(i) {
  table(pred = pred_matrix[, i],true = y_true)})

acc_iter <- sapply(1:n_iter, function(i) {
  mean(pred_matrix[, i] == y_true)
})


df <- data.frame(acc = acc_iter)
mean(acc_iter)
sd(acc_iter)
quantile(acc_iter, c(0.025, 0.5, 0.975))

gg1 <- ggplot(df, aes(x = acc)) +
  geom_histogram(bins = 16, fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  geom_vline(xintercept = acc_ref, color = "red", size = 1) +
  ggtitle("Monte Carlo distribution accuracy") +
  xlab("Accuracy") +
  annotate("text",
           x = min(df$acc),
           y = Inf,
           label = paste0("\U0394 Acc = ", round(delta_acc_mean,3),
                          "\nMean = ", round(mean(acc_iter),3),
                          "\nSD = ", round(sd(acc_iter),3)),
           hjust = 0,
           vjust = 1.5,
           size = 3.5) +
  theme_minimal()

gg1
ggsave("./analysis/data/derived_data/model/Monte_Carlo_distribution.pdf", gg1, width = 8, height = 6)
ggsave("./analysis/data/derived_data/model/Monte_Carlo_distribution.png", gg1, width = 8, height = 6)

# 07.8 Model visualisation =====================================================

# Confert the list into an array
conf_array <- simplify2array(conf_list)
conf_mean <- apply(conf_array, c(1,2), mean)
conf_sd   <- apply(conf_array, c(1,2), sd)

labels <- matrix(
  paste0(
    round(conf_mean, 1),
    "\n±",
    round(conf_sd, 1)
  ),
  nrow = nrow(conf_mean)
)

# transformer en format long
df_mean <- melt(conf_mean)
df_sd   <- melt(conf_sd)

df <- cbind(df_mean, df_sd$value)
colnames(df) <- c("Pred", "True", "Mean", "SD")

# Create a label integrating sd
df$label <- paste0(round(df$Mean,1), "\n±", round(df$SD,1))
df$label[df$label == "0\n±0"] <- "0"


# Confusion matrix plot
p1 <- ggplot(df, aes(x = True, y = Pred, fill = Mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 2.5) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(title = "Confusion Matrix",
       x = "Actual Class",
       y = "Predicted Class") +
  scale_y_discrete(limits = rev) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1
ggsave("./analysis/data/derived_data/model/Confusion_matrice.pdf", p1, width = 12, height = 8)
ggsave("./analysis/data/derived_data/model/Confusion_matrice.png", p1, width = 12, height = 8)

predictions <- as.data.frame(apply(results_selected, c(1,2), mean))
predictions$predicted <- colnames(predictions)[max.col(predictions, ties.method = "random")]
predictions$real <- y_true

# Function to plot ROC curves for all labels
plot_roc_curves <- function(data) {
  # Get probability columns
  prob_cols <- c("LP", "RG", "CCcm", "FL", "KZ", "VR", "CM", "CCrg")
  
  roc_data_list <- list()
  auc_values <- numeric()
  
  for(label in prob_cols) {
    true_binary <- data$real == label
    
    roc_obj <- roc(true_binary, data[[label]])
    auc_values[label] <- auc(roc_obj)
    
    roc_data_list[[label]] <- data.frame(
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      Label = label
    )
  }
  
  
  all_roc_data <- do.call(rbind, roc_data_list)
  
  # Create AUC annotation text
  auc_text <- paste(names(auc_values), 
                    sprintf("%.3f", auc_values), 
                    sep = ": ", 
                    collapse = "\n")
  
  # Plot ROC curves
  ggplot(all_roc_data, aes(x = FPR, y = TPR, color = Label)) +
    geom_line() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    theme_minimal() +
    labs(title = "ROC Curves for Each Class",
         x = "False Positive Rate",
         y = "True Positive Rate") +
    annotate("text", x = 0.75, y = 0.25, 
             label = paste("AUC values:\n", auc_text),
             hjust = 0, vjust = 0) +
    theme(legend.position = "right")
}

# Function to plot prediction confidence
plot_prediction_confidence <- function(data) {
  confidence_data <- data %>%
    dplyr::select(LP, RG, CCcm, FL, KZ, VR, CM, CCrg) %>%
    as.matrix() %>%
    apply(1, max) %>%
    data.frame(Confidence = ., 
               Correct = data$predicted == data$real)
  
  ggplot(confidence_data, aes(x = Confidence, fill = Correct)) +
    geom_histogram(position = "dodge", bins = 30) +
    scale_fill_manual(values = c("red", "green")) +
    theme_minimal() +
    labs(title = "Prediction Confidence Distribution",
         x = "Confidence Score",
         y = "Count",
         fill = "Correct Prediction")
}

sd_iter <- as.data.frame(apply(results_selected, c(1,2), sd))
pred_class_idx <- max.col(predictions[,c(1:11)])
predictions$sd <- sd_iter[cbind(1:nrow(sd_iter), pred_class_idx)]
prob_cols <- colnames(predictions[1:11])

# Function to plot class distribution
plot_class_distribution <- function(data) {
  actual_counts <- table(factor(data$real, levels = prob_cols))
  predicted_counts <- table(factor(data$predicted, levels = prob_cols))
  sd_counts <- predictions %>%
    group_by(predicted) %>%
    summarise(SD_mean = mean(sd, na.rm = TRUE),.groups = "drop")
  
  distribution_data <- data.frame(
    Class = prob_cols,
    Actual = as.vector(actual_counts),
    Predicted = as.vector(predicted_counts),
    SD = sd_counts$SD_mean
  ) %>%
    pivot_longer(cols = c("Actual", "Predicted"),
                 names_to  = "Type", 
                 values_to = "Count")
  
  pos <- position_dodge(width = 0.9)
  ggplot(distribution_data, aes(x = Class, y = Count, fill = Type)) +
    geom_bar(stat = "identity", position = pos) +
    
    geom_errorbar(
      data = subset(distribution_data, Type == "Predicted"),
      aes(
        ymin = Count - (SD*100),
        ymax = Count + (SD*100)
      ),
      position = pos,
      width = 0.2,
      color = "grey30"
    ) +
    
    geom_point(
      data = subset(distribution_data, Type == "Predicted"),
      aes(y = Count),
      position = pos,
      size = 2,
      color = "grey20"
    ) +
    scale_fill_manual(values = c("palevioletred", "steelblue")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Class Distribution: Actual vs Predicted",
         x = "Class",
         y = "Count")
}

p2 <- plot_roc_curves(predictions)  
p3 <- plot_prediction_confidence(predictions)
p4 <- plot_class_distribution(predictions)

# Arrange plots in a grid
finalplot <- grid.arrange(p1, p2, p3, p4, ncol = 2)

ggsave("./analysis/data/derived_data/model/Plots_results_RF.pdf", finalplot, width = 12, height = 8)
ggsave("./analysis/data/derived_data/model/Plots_results_RF.png", finalplot, width = 12, height = 8)

# 07.9 Model metrics =====================================================

# Calculate and print accuracy
accuracy <- sum(predictions$predicted == predictions$real) / nrow(predictions)
F1 <- F1_Score(predictions$predicted, predictions$real)

F1_iter <- sapply(1:n_iter, function(i) {
  F1_Score(pred_matrix[, i], y_true)
})

cat("Overall Accuracy:", round(accuracy * 100, 2), " \U2213 ", round(sd(acc_iter)* 100, 3), 
    " and F1 score", round(F1 * 100, 2), " \U2213 ", round(sd(F1_iter)* 100, 3), "%\n")

# Calculate per-class metrics
class_metrics <- data.frame(
  Class = unique(predictions$real),
  Precision = NA,
  Recall = NA,
  F1_Score = NA
)

class_metrics_iter <- data.frame(
  Class = unique(predictions$real),
  Precision = NA,
  Recall = NA,
  F1_Score = NA
)

for(cls in class_metrics$Class) {
  true_pos <- sum(predictions$predicted == cls & predictions$real == cls)
  true_pos_iter <- sapply(1:n_iter, function(i) {
    sum(pred_matrix[, i] == cls & predictions$real == cls)
  })
  
  false_pos <- sum(predictions$predicted == cls & predictions$real != cls)
  false_pos_iter <- sapply(1:n_iter, function(i) {
    sum(pred_matrix[, i] == cls & predictions$real != cls)
  })
  
  false_neg <- sum(predictions$predicted != cls & predictions$real == cls)
  false_neg_iter <- sapply(1:n_iter, function(i) {
    sum(pred_matrix[, i] != cls & predictions$real == cls)
  })
  
  true_neg <- sum(predictions$predicted != cls & predictions$real != cls)
  true_neg_iter <- sapply(1:n_iter, function(i) {
    sum(pred_matrix[, i] != cls & predictions$real != cls)
  })
  
  precision <- true_pos / (true_pos + false_pos)
  precision_iter <- true_pos_iter / (true_pos_iter + false_pos_iter)
  
  recall <- true_pos / (true_pos + false_neg)
  recall_iter <- true_pos_iter / (true_pos_iter + false_neg_iter)
  
  f1 <- 2 * (precision * recall) / (precision + recall)
  f1_iter <- 2 * (precision_iter * recall_iter) / (precision_iter + recall_iter)
  
  accuracy <- sum(true_pos + true_neg) / (true_pos + false_pos + true_neg + false_neg)
  accuracy_iter <- (true_pos_iter + true_neg_iter) /(true_pos_iter + false_pos_iter + true_neg_iter + false_neg_iter)
  
  class_metrics[class_metrics$Class == cls, c("Precision", "Recall", "F1_Score", "Accuracy")] <- 
    c(precision, recall, f1, accuracy)
  
  class_metrics_iter[class_metrics_iter$Class == cls, c("Precision", "Recall", "F1_Score", "Accuracy")] <- 
    c(sd(precision_iter), sd(recall_iter), sd(f1_iter), sd(accuracy_iter))
}

row.names(class_metrics) <- class_metrics[,1]
class_metrics <- class_metrics[,-1]
write.table(class_metrics, "./analysis/data/derived_data/model/Metrics_RF.txt")
print(round(class_metrics, 3))

row.names(class_metrics_iter) <- class_metrics_iter[,1]
class_metrics_iter <- class_metrics_iter[,-1]
write.table(class_metrics_iter, "./analysis/data/derived_data/model/Metrics_SD_RF.txt")
print(round(class_metrics_iter, 3))

class_metrics <- as.matrix(class_metrics)
class_metrics_iter <- as.matrix(class_metrics_iter)

formatted_table <- matrix(
  paste0(
    sprintf("%.3f", class_metrics),
    " ± ",
    sprintf("%.3f", class_metrics_iter)
  ),
  nrow = nrow(class_metrics),
  dimnames = dimnames(class_metrics)
)
write.table(formatted_table, "./analysis/data/derived_data/model/Metrics_summed.txt")

# 07.10 Look at models covariates influences ====================================

importance_df <- as.data.frame(FinalRFCon[[1]]$finalModel$variable.importance)


# Convert row names (variables) to a column for ggplot
importance_df$scale <- (importance_df[,1]/sum(importance_df[,1])*100)
importance_df$Variable <- rownames(importance_df)

# Plot using ggplot2
gg1 <- ggplot(importance_df, aes(x = reorder(Variable, scale), y = scale)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  coord_flip() +
  xlab("Covariates") +
  ylab("Importance scaled (%)") +
  ggtitle("Variable Importance from RF Model") +
  theme_minimal()

gg1
ggsave("./analysis/data/derived_data/model/Variables_importance.png", gg1, width = 12, height = 8)
ggsave("./analysis/data/derived_data/model/Variables_importance.pdf", gg1, width = 12, height = 8)

rm(list=ls())

# 08 Prediction ################################################################
# 08.1 Import the data =========================================================
period <- c("EB", "MB", "LB", "IA")
covariates <- rast("./analysis/data/raw_data/predictors.tif")
load("./analysis/data/derived_data/save/Models.RData")
load("./analysis/data/derived_data/save/Gaussian_noice.RData")
load("./analysis/data/derived_data/save/Monte_Carlo.RData")

# If you made a VIF or RFE selection
selection <- colnames(FinalRFCon[[1]]$trainingData[,2:ncol(FinalRFCon[[1]]$trainingData)])


# 08.2 Set maps ================================================================

start_time <- Sys.time()
for (era in period) {
  
  r <- rast(paste0("./analysis/data/raw_data/CHELSA_PAST/cov_", era,".tif"))
  raster_stack <- c(covariates[[1:7]], r, covariates$DEM)
  bio_names <- sprintf("bio%02d", 1:19)
  names(raster_stack) <- c(names(covariates[[1:7]]), bio_names, "DEM")
  raster_stack <- raster_stack[[selection]]
  raster_stack <- aggregate(raster_stack, fact = 10)
  terra::writeRaster(raster_stack, paste0("./analysis/data/raw_data/Cov_selected_",era,".tif"), overwrite = TRUE)
  gaus.noice <- gaus.noice[selection]
  rm(raster_stack, r)
  
  # 08.3 Run predictions ========================================================= 
  Cov_rast <- rast(paste0("./analysis/data/raw_data/Cov_selected_",era,".tif"))
  
  x_df <- as.data.frame(Cov_rast, xy = TRUE)
  predicted_block <- predict(FinalRFCon[[1]]$finalModel, data = x_df[,-c(1:2)])
  
  class_names <- colnames(predicted_block$predictions)
  sum_probs <- matrix(0, nrow = nrow(x_df), ncol = length(class_names))
  colnames(sum_probs) <- class_names
  
  sum_probs <- predicted_block$predictions
  
  cli_progress_bar(
    format = paste0("Running predictions {era} {cli::pb_bar} {cli::pb_percent} [{cli::pb_current}/{cli::pb_total}] | \ ETA: {cli::pb_eta} - Time elapsed: {cli::pb_elapsed_clock}"),
    total = n_iter, clear = FALSE)
  
  for (i in 1:n_iter) {
    features <- as.data.frame(Cov_rast)
    Sigma <- Sigma_base * 0.0001
    
    n_obs <- nrow(features)
    num_cols <- names(features)[sapply(features, is.numeric)]
    
    set.seed(seed.list[i])
    noise <- mvrnorm(n_obs, mu = rep(0, length(selection)), Sigma = Sigma)
    colnames(noise) <- selection
    
    df_noisy <- features
    df_noisy[, num_cols] <- features[, num_cols] + noise
    
    pred <- predict(FinalRFCon[[1]]$finalModel, data = df_noisy)
    probs <- pred$predictions
    
    sum_probs <- sum_probs + probs
    cli_progress_update()  
  }
  
  mean_probs <- sum_probs / n_iter
  
  final_class <- colnames(mean_probs)[max.col(mean_probs, ties.method = "random")]
  uncertainty <- 1 - apply(mean_probs, 1, max)
  x_pred <- cbind(x_df[,1:2], final_class, uncertainty)
  
  mapping <- c("LP" = 2, "RG" = 3, "CCcm" = 4, "FL" = 5, "GY" = 6, "KZ" = 7, "VR" = 9, "CM" = 10, "CCrg" = 11, "CCvrfl" = 12, "CMvr" = 13 )
  x_pred$final_class <- as.numeric(mapping[x_pred$final_class])
  x_raster <- rast(x_pred, type = "xyz")
  crs(x_raster) <- "EPSG:32638"
  plot(x_raster$final_class)
  plot(x_raster$uncertainty)
  terra::writeRaster(x_raster, paste0("./analysis/data/derived_data/maps/Palaeosol_",era,".tif"), overwrite = TRUE)
  cli_progress_done()
}

end_time <- Sys.time()
cat("Time spend for mapping :", round(difftime(end_time, start_time, units="mins"), 2), "minutes\n")

# 08.4 Render the images ======================================================= 

x_df <- data.frame()

for (era in period) {
  
  r <- rast(paste0("./analysis/data/derived_data/maps/Palaeosol_",era,".tif")) 
  names(r) <- c("prediction", "uncertainty")
  r_df <- as.data.frame(r, xy = TRUE)
  r_df$era <- era
  x_df <- rbind(x_df, r_df)
  
}


soil_levels <- c(2,3,4,5,7,9,10,13)

labels <- c(
  "2" = "Leptosol",
  "3" = "Regosol",
  "4" = "Calcisol\nminor Cambisol",
  "5" = "Fluvisol",
  "7" = "Kastanozem",
  "9" = "Vertisol",
  "10" = "Cambisol",
  "13" = "Cambisol\nminor Vertisol")

soil_colors <- c(
  "2" = "#440154",
  "3" = "#482475",
  "4" = "#22a884",
  "5" = "#7ad151",
  "7" = "#bddf26",
  "9" = "#fde725",
  "10" = "#44bf70",
  "13" = "#414487")

fig_list_pred <- list()
fig_list_uncer <- list()

for (era in period) {
  x_era <- x_df[x_df$era == era,] # Select the period
  h_dist <- ifelse(era == "IA", -1.25, -0.95) # Small problem with the h_dist for IA wihtout reason
  gg1 <- ggplot(x_era, aes(x = x, y = y, fill = factor(prediction))) +
    geom_raster() +
    coord_equal() +
    xlab("") +
    ylab("") +
    scale_fill_manual(
      name   = "Soil classes",
      values = soil_colors,
      breaks = soil_levels,
      labels = labels,
      drop   = FALSE
    ) +
    annotate("text",
             x = -Inf,
             y = Inf,
             label = era,
             hjust = h_dist,
             vjust = 2.5,
             size = 6,
             color = "white",
             fontface = "bold") +
    theme_void()
  
  gg2 <- ggplot(x_era, aes(x = x, y = y, fill = uncertainty)) +
    geom_raster() +
    coord_equal() +
    xlab("") +
    ylab("") +
    scale_fill_viridis_c(name   = "Uncertainty", begin = 0, end = 0.86, option = "E", na.value = "grey50") +
    annotate("text",
             x = -Inf,
             y = Inf,
             label = era,
             hjust = h_dist,
             vjust = 2.5,
             size = 6,
             color = "white",
             fontface = "bold") +
    theme_void()
  
  fig_list_pred[[era]] <- gg1
  fig_list_uncer[[era]] <- gg2
}


pred_panel <- wrap_plots(fig_list_pred, nrow = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

pred_panel
ggsave("./analysis/data/derived_data/maps/Palaeosol_pred_final_fig.pdf", pred_panel, width = 12, height = 8)
ggsave("./analysis/data/derived_data/maps/Palaeosol_pred_final_fig.png", pred_panel, width = 12, height = 8)


uncer_panel <- wrap_plots(fig_list_uncer, nrow = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave("./analysis/data/derived_data/maps/Palaeosol_uncer_final_fig.pdf", uncer_panel, width = 12, height = 8)
ggsave("./analysis/data/derived_data/maps/Palaeosol_uncer_final_fig.png", uncer_panel, width = 12, height = 8)


rm(list = ls())

# 09 Prediction for GIF images #################################################
# 09.1 Import the data =========================================================
period <- seq(-4400, -2300, by = 100)
covariates <- rast("./analysis/data/raw_data/predictors.tif")
load("./analysis/data/derived_data/save/Models.RData")
load("./analysis/data/derived_data/save/Gaussian_noice.RData")
load("./analysis/data/derived_data/save/Monte_Carlo.RData")

# If you made a VIF or RFE selection
selection <- colnames(FinalRFCon[[1]]$trainingData[,2:ncol(FinalRFCon[[1]]$trainingData)])


# 09.2 Set maps ================================================================

start_time <- Sys.time()
for (era in period) {
  r <- rast(paste0("./analysis/data/raw_data/CHELSA_PAST/GIF_seq/", era,".tif"))
  raster_stack <- c(covariates[[1:7]], r, covariates$DEM)
  bio_names <- sprintf("bio%02d", 1:19)
  names(raster_stack) <- c(names(covariates[[1:7]]), bio_names, "DEM")
  raster_stack <- raster_stack[[selection]]
  raster_stack <- aggregate(raster_stack, fact = 10)
  
  # 09.3 Run predictions =========================================================
  
  x_df <- as.data.frame(raster_stack, xy = TRUE)
  predicted_block <- predict(FinalRFCon[[1]]$finalModel, data = x_df[,-c(1:2)])
  
  class_names <- colnames(predicted_block$predictions)
  sum_probs <- matrix(0, nrow = nrow(x_df), ncol = length(class_names))
  colnames(sum_probs) <- class_names
  
  sum_probs <- predicted_block$predictions
  
  cli_progress_bar(
    format = paste0("Running predictions {era} {cli::pb_bar} {cli::pb_percent} [{cli::pb_current}/{cli::pb_total}] | \ ETA: {cli::pb_eta} - Time elapsed: {cli::pb_elapsed_clock}"),
    total = n_iter, clear = FALSE)
  
  for (i in 1:n_iter) {
    features <- as.data.frame(raster_stack)
    Sigma <- Sigma_base * 0.0001
    
    n_obs <- nrow(features)
    num_cols <- names(features)[sapply(features, is.numeric)]
    
    set.seed(seed.list[i])
    noise <- mvrnorm(n_obs, mu = rep(0, length(selection)), Sigma = Sigma)
    colnames(noise) <- selection
    
    df_noisy <- features
    df_noisy[, num_cols] <- features[, num_cols] + noise
    
    pred <- predict(FinalRFCon[[1]]$finalModel, data = df_noisy)
    probs <- pred$predictions
    
    sum_probs <- sum_probs + probs
    cli_progress_update()  
  }
  
  mean_probs <- sum_probs / n_iter
  
  final_class <- colnames(mean_probs)[max.col(mean_probs, ties.method = "random")]
  uncertainty <- 1 - apply(mean_probs, 1, max)
  x_pred <- cbind(x_df[,1:2], final_class, uncertainty)
  
  mapping <- c("LP" = 2, "RG" = 3, "CCcm" = 4, "FL" = 5, "GY" = 6, "KZ" = 7, "VR" = 9, "CM" = 10, "CCrg" = 11, "CCvrfl" = 12, "CMvr" = 13 )
  x_pred$final_class <- as.numeric(mapping[x_pred$final_class])
  x_raster <- rast(x_pred, type = "xyz")
  crs(x_raster) <- "EPSG:32638"
  plot(x_raster$final_class)
  plot(x_raster$uncertainty)
  terra::writeRaster(x_raster, file = paste0("./analysis/data/derived_data/maps/GIF/Palaeosol_",era,".tif"), overwrite = TRUE)
  cli_progress_done()
}

end_time <- Sys.time()
cat("Time spend for GIF :", round(difftime(end_time, start_time, units="mins"), 2), "minutes\n")

rm(list = ls())

# 09.4 Render the images =======================================================
period <- seq(-4400, -2300, by = 100)

files <- list.files("./analysis/data/derived_data/maps/GIF", full.names = TRUE)
maps <- rast(rev(files)) 

x_df <- data.frame()

for (era in period) {
  
  r <- rast(paste0("./analysis/data/derived_data/maps/GIF/Palaeosol_",era,".tif")) 
  names(r) <- c("prediction", "uncertainty")
  r_df <- as.data.frame(r, xy = TRUE)
  r_df$era <- era
  x_df <- rbind(x_df, r_df)
  
}

soil_levels <- c(2,3,4,5,6,7,9,10,11,12,13)

labels <- c(
  "2" = "Leptosol",
  "3" = "Regosol",
  "4" = "Calcisol\nminor Cambisol",
  "5" = "Fluvisol",
  "6" = "Gleyosol",
  "7" = "Kastanozem",
  "9" = "Vertisol",
  "10" = "Cambisol",
  "11" = "Calcisol\nminor Cambisol",
  "12" = "Calcisol minor \nVertisol and Fluvisol",
  "13" = "Cambisol\nminor Vertisol")

soil_colors <- c(
  "2" = "#440154",
  "3" = "#482475",
  "4" = "#22a884",
  "5" = "#7ad151",
  "6" = "#355f8d",
  "7" = "#bddf26",
  "9" = "#fde725",
  "10" = "#44bf70",
  "11" = "#2a788e",
  "12" = "#21918c",
  "13" = "#414487")

for (era in period) {
  x_era <- x_df[x_df$era == era,]
  gg1 <- ggplot(x_era, aes(x = x, y = y, fill = factor(prediction))) +
    geom_raster() +
    coord_equal() +
    xlab("") +
    ylab("") +
    scale_fill_manual(
      name   = "Soil classes",
      values = soil_colors,
      breaks = soil_levels,
      labels = labels,
      drop   = FALSE
    ) +
    annotate("text",
             x = -Inf,
             y = Inf,
             label = era,
             hjust = -0.5,
             vjust = 2.5,
             size = 6,
             color = "black",
             fontface = "bold") +
    theme_void()
  
  gg2 <- ggplot(x_era, aes(x = x, y = y, fill = uncertainty)) +
    geom_raster() +
    coord_equal() +
    xlab("") +
    ylab("") +
    scale_fill_viridis_c(name   = "Uncertainty", begin = 0, end = 0.86, option = "E", na.value = "grey50") +
    annotate("text",
             x = -Inf,
             y = Inf,
             label = era,
             hjust = -0.5,
             vjust = 2.5,
             size = 6,
             color = "black",
             fontface = "bold") +
    theme_void()
  
  gg <- gg1 + gg2
  ggsave(paste0("./analysis/data/derived_data/maps/GIF_fig/",era, ".png"), gg, width = 12, height = 8)
}
