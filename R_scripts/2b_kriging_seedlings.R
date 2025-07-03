### kriging of seedling density ####


library("xlsx")
library("tidyverse")
library("graph4lg")
library("sf")
library("sp")
library("raster")
library("ggplot2")
library("tidyverse")
library("terra")
library("GeoModels")
require("fields")
set.seed(89)

path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"

### Allenwiller ####
all_data=fread(paste0(path,"Data_clean/Allenwiller_circleplot_info_withcoord.csv"))
hist(all_data$`CountedC1(new)`)
hist(log(all_data$`CountedC1(new)`)) # data transformation
## add a small constant to 0 to avoid -Inf
epsilon <- 0.0001  # You can also try smaller values like 0.1 or 0.01
log_data <- log(all_data$`CountedC1(new)` + epsilon)

seedlings=matrix(log_data) # vector of data to krige
coords=data.matrix(all_data[,c("x", "y") ]) # coordinates of the Circle plots (1 row / circle plot)
plot(coords)

# correlation model and correlation parameters (relationship between the points, the relationship between 2 points is exponential)
corrmodel="Exp"
CorrParam(corrmodel) # parameters to set in the model
scale=5 # spatial window to take the data for the estimation of the parameters
NuisParam("Gaussian") # parameters to set (?)
mean=mean(log_data) 
nugget=0.2
sill=150 
start=list(mean=mean, sill=sill, scale=scale) # parameters estimated with the regression
fixed=list(nugget=nugget) # fixed parameters

# running the function separately (fitting the model to the points)
fitML <- GeoModels::GeoFit(data=seedlings, # data to krige 
                           coordx=coords, # coordinates
                           corrmodel=corrmodel, # exponential
                           model="Gaussian", # random field involving Gaussian probability density function of the variables (numbers are drawn from a Gaussian probability distribution)
                           start=start, # starting parameters (mean, scale, sill)
                           maxdist = 200, # may dist used for the model
                           fixed=fixed) # fixed parameters (nugget)

# info about the max likelihood approach
fitML # if close to 0 is good (probability? that the observation is matching/close to the fit)
# checking model assumptions - residuals
res=GeoResiduals(fitML) 
# ST variogram model
vario = GeoVariogram(data=res$data,coordx=coords,maxdist = 200) # max dist?
# marginal distribution (residuals)
GeoQQ(res)
GeoCovariogram(res,vario=vario,show.vario=TRUE,pch=20)
GeoQQ(res,type="D",ylim=c(0,0.5),breaks=20)

# Prediction
grid_al <- vect(paste0(path,"Data_clean/Allenwiller_grid_r4m_9844patches.shp"), crs = "EPSG:3035")
centr <- centroids(grid_al)
plot(centr)
centr_m <- crds(centr)

# estimated parameters
param_est=as.list(c(fitML$param,fixed))

### from here
# optimal linear prediction
pr = GeoKrig(data=seedlings,
             coordx=coords,  
             corrmodel=corrmodel,
             sparse=FALSE,
             model="Gaussian",
             mse=TRUE,
             loc=centr_m,
             param=param_est)

# plot
par(mfrow=c(1,2))
quilt.plot(coords[,1],coords[,2],exp(seedlings),main ="Seedling density: Observed - Allenwiller")
quilt.plot(centr_m[,1],centr_m[,2],exp(pr$pred),main ="Seedling density: Prediction - Allenwiller") # data transformed back to original scale (exp(log))
dev.off()

# saving kriging results
kriging_res=data.frame(x_centroid=centr_m[,1],y_centroid=centr_m[,2],predicted_log=pr$pred,predicted=exp(pr$pred))
write.csv(kriging_res,paste0(path,"Data_clean/Allenwiller_kriging_results.csv"))


### Waldi ###
all_data=fread(paste0(path,"Data_clean/Waldi_circleplot_info_withcoord.csv"))
hist(all_data$NClass1Class3)
hist(log(all_data$NClass1Class3)) # data transformation
## add a small constant to 0 to avoid -Inf
epsilon <- 0.0001  # You can also try smaller values like 0.1 or 0.01
log_data <- log(all_data$NClass1Class3 + epsilon)

seedlings=matrix(log_data) # vector of data to krige
coords=data.matrix(all_data[,c("x", "y") ]) # coordinates of the Circle plots (1 row / circle plot)
plot(coords)

## check for missing coordinates
anyNA(coords)
all_data[!complete.cases(coords), ]
# calculate coordinates from mother trees coords
#convert orientation degrees to radians
theta_rad <- as.numeric(all_data$OrientationDegree) * pi / 180

# get x and y offsets (note: North = 0°, East = 90°, etc.)
all_data[, x := ifelse(is.na(x),x_mothertree + as.numeric(DistanceM) * sin(theta_rad),x)]
all_data[, y := ifelse(is.na(y),y_mothertree + as.numeric(DistanceM) * cos(theta_rad),y)]
anyNA(all_data[, c("x", "y")]) 

# rebuild coords
coords=data.matrix(all_data[,c("x", "y") ]) # coordinates of the Circle plots (1 row / circle plot)
plot(coords)

# correlation model and correlation parameters (relationship between the points, the relationship between 2 points is exponential)
corrmodel="Exp"
CorrParam(corrmodel) # parameters to set in the model
scale=5 # spatial window to take the data for the estimation of the parameters
NuisParam("Gaussian") # parameters to set 
mean=mean(log_data) 
nugget=0.2
sill=150 
start=list(mean=mean, sill=sill, scale=scale) # parameters estimated with the regression
fixed=list(nugget=nugget) # fixed parameters

# running the function separately (fitting the model to the points)
fitML <- GeoModels::GeoFit(data=seedlings, # data to krige 
                           coordx=coords, # coordinates
                           corrmodel=corrmodel, # exponential
                           model="Gaussian", # random field involving Gaussian probability density function of the variables (numbers are drawn from a Gaussian probability distribution)
                           start=start, # starting parameters (mean, scale, sill)
                           maxdist = 200, # max dist used for the model
                           fixed=fixed) # fixed parameters (nugget)

# info about the max likelihood approach
fitML # if close to 0 is good (probability? that the observation is matching/close to the fit)
# checking model assumptions - residuals
res=GeoResiduals(fitML) 
# ST variogram model
vario = GeoVariogram(data=res$data,coordx=coords,maxdist = 200) # max dist?
# marginal distribution (residuals)
GeoQQ(res)
GeoCovariogram(res,vario=vario,show.vario=TRUE,pch=20)
GeoQQ(res,type="D",ylim=c(0,0.5),breaks=20)

# Prediction
grid_wal <- vect(paste0(path,"Data_clean/Waldi_grid_4m_5396patches.shp"), crs = "EPSG:3035")
centr <- centroids(grid_wal)
plot(centr)
centr_m <- crds(centr)

# estimated parameters
param_est=as.list(c(fitML$param,fixed))

### from here
# optimal linear prediction
pr = GeoKrig(data=seedlings,
             coordx=coords,  
             corrmodel=corrmodel,
             sparse=FALSE,
             model="Gaussian",
             mse=TRUE,
             loc=centr_m,
             param=param_est)

# plot
par(mfrow=c(1,2))
quilt.plot(coords[,1],coords[,2],exp(seedlings),main ="Seedling density: Observed - Waldi")
quilt.plot(centr_m[,1],centr_m[,2],exp(pr$pred),main ="Seedling density: Prediction - Waldi") # data transformed back to original scale (exp(log))
dev.off()

# saving kriging results
kriging_res=data.frame(x_centroid=centr_m[,1],y_centroid=centr_m[,2],predicted_log=pr$pred,predicted=exp(pr$pred))
write.csv(kriging_res,paste0(path,"Data_clean/Waldi_kriging_results.csv"))


