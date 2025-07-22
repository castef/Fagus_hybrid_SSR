# Convert Coordinates from original to WGS84 and metric system #

library(xlsx)
library(terra)

rm(list=ls())

setwd(dir="~/Data_clean")

### Waldi 
# circle plots
waldi_cp_o <- read.xlsx("20221206_CoordinatesCirclePlots_Wäldi.xlsx", sheetIndex = 1)
waldi_cp_o = waldi_cp_o[, c(2,5,6)]
which(duplicated(waldi_cp_o$CirclePlotID))
# define original coordinate system
waldi_cp_v <- vect(waldi_cp_o, geom = c("X2..CirclePlot.", "Y2..CirclePlot."), crs = "+init=epsg:21781")
# conversion to metric system
waldi_m = project(x=waldi_cp_v,y= "+init=epsg:3035")
# conversion to WGS84
waldi_wgs84 <- project(x=waldi_cp_v,y= "+proj=longlat +datum=WGS84")
# saving 
write.csv(as.data.frame(waldi_m,geom="XY"), file = "Waldi_circlecplots_Coordinates_metric.csv", quote=F)

# adults
waldi_a_o <- read.xlsx("20221202_CoordinatesAdults_Wäldi.xlsx", sheetIndex = 1)
waldi_a_o = waldi_a_o[, c(1:3)]
waldi_a_o=subset(waldi_a_o,ProbeID!="na")
# define original coordinate system
waldi_a_v <- vect(waldi_a_o, geom = c("Easting","Northing"), crs = "+init=epsg:21781")
# conversion to metric system
waldi_m = project(x=waldi_a_v,y= "+init=epsg:3035")
# conversion to WGS84
waldi_wgs84 <- project(x=waldi_a_v,y= "+proj=longlat +datum=WGS84")
# saving 
write.csv(as.data.frame(waldi_m,geom="XY"), file = "Waldi_adults_Coordinates_metric.csv", quote=F)

### Allenwiller 
# circle plots
allen_cp_o <- read.xlsx("20221206_CoordinatesCirclePlots_Allenwiller.xlsx", sheetIndex = 1)
allen_cp_o = allen_cp_o[, c(4:6)]
which(duplicated(allen_cp_o$CirclePlotID))

# define original coordinate system
allen_cp_v <- vect(allen_cp_o, geom = c("X2..CirclePlot.", "Y2..CirclePlot."), crs = "+init=epsg:32632")
# conversion to metric system
allen_m = project(x=allen_cp_v,y= "+init=epsg:3035")
# conversion to WGS84
allen_wgs84 <- project(allen_m, "+proj=longlat +datum=WGS84")
# saving 
write.csv(as.data.frame(allen_m,geom="XY"), file = "Allenwiller_circlecplots_Coordinates_metric.csv", quote=F)

# adults
allen_a_o <- read.xlsx("20221202_CoordinatesAdults_Allenwiller.xlsx", sheetIndex = 1)
allen_a_o = allen_a_o[, c(1:3)]
allen_a_o=subset(allen_a_o,ProbeID!="na")
# remove circle plots and other mistakes
allen_a_o = allen_a_o[grepl("\\d+$", allen_a_o$ProbeID),]

# define original coordinate system
allen_a_v <- vect(allen_a_o, geom = c("Easting","Northing"), crs = "+init=epsg:32632")
# conversion to metric system
allen_m = project(x=allen_a_v,y= "+init=epsg:3035")
# conversion to WGS84
allen_wgs84 <- project(x=allen_a_v,y= "+proj=longlat +datum=WGS84")
# saving 
write.csv(as.data.frame(allen_m,geom="XY"), file = "Allenwiller_adults_Coordinates_metric.csv", quote=F)


