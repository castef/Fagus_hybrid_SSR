### script to create input files for NEMO simulation #####

library(dplyr)
library(grid)
library(raster)
library(terra)
library(ggplot2)
library(readr)
library(xlsx)
library(adegenet)
library(tidyr)
library(data.table)

rm(list=ls())

# change paths as needed
setwd("~/Data_clean/")
#setwd("C:/Users/Camilla/Dropbox (Old (1))/Dropbox/WSL_PhD/Projects/Hybridization/Data_clean/")

# path to get previous infos
path="~/Data_clean/"

# path where to create and store the input files
input_path = "~/NEMO/Input_files/"

# path to run the simulation (where to move the final input files)
exe_path = "~/nemoage0.32.6/exe/"

####  sampling scheme #####

# WALDI #

wal_sam <- vect("~/waldi_sampling_points.shp")
wal_sam=wal_sam[!is.na(geom(wal_sam)[,"x"]),] #remove a point without coordinates
# file of adults and circle plots coordinates
trees= fread("~/Waldi_ACp_Coordinates_21781.csv")
# convert to SpatVector 
trees_v <- vect(trees, geom = c("X", "Y"), crs = crs(wal_sam))

#check before coordiante transformation
plot(wal_sam)
points(trees_v,col="red")

# conversion to metric system
wal_sam = project(x=wal_sam,y= "epsg:3035")
trees_v = project(x=trees_v,y= "epsg:3035")
#check after coordiante transformation
plot(wal_sam)
points(trees_v,col="red")
# create buffer with the radius of the circle plots
wal_sam_buffer=buffer(wal_sam,wal_sam$sam_radius)
plot(wal_sam_buffer, main = "waldi sampling scheme")
points(trees_v,col="red", cex = 0.4)

# save shapefile of sampling scheme
writeVector(wal_sam_buffer, "Waldi_sampling_scheme.shp", overwrite = T)

# ALLENWILLER #

# read shapefile with circleplots
al_sam <- vect("~/Allenwiller_sampling_points.shp")
al_sam=al_sam[!is.na(geom(al_sam)[,"x"]),] #remove a point without coordinates
# file of adults and circle plots coordinates
trees1= read.xlsx("~/Allenwiller_ACp_Coordinates_UTM32N_32632.xlsx", sheetIndex = 1)
# convert to SpatVector 
trees1_v <- vect(trees1, geom = c("Easting", "Northing"), crs = crs("EPSG:32632"))
plot(trees1_v)

#check before coordiante transformation
plot(al_sam)
points(trees1_v,col="red")

# conversion to metric system
al_sam = project(x=al_sam,y= "epsg:3035")
trees1_v = project(x=trees1_v,y= "epsg:3035")
#check after coordiante transformation
plot(al_sam)
points(trees1_v,col="red")

# add radius of the circle plots
al_sam$sam_radius = 2
# buffer the shapefile with the circle plot radius
al_sam_buffer=buffer(al_sam,al_sam$sam_radius)
plot(al_sam_buffer, main = "Allenwiller sampling scheme")
points(trees1_v,col="red")

##NB shpfiles have a limit of colnames of 10 characters

# write shapefile in Data clean
writeVector(al_sam_buffer, "Allenwiller_sampling_scheme.shp", overwrite = T)

# plots of the sampling scheme
library(grid)

# function to keep the same style 
plot_sampling_scheme <- function(buffer_obj, title, add_scale = TRUE) {
  plot(buffer_obj,main = title,col = "lightblue",border = "darkblue",lwd = 0.5)
  # scale bar
  if (add_scale) {
    ext <- ext(buffer_obj)  # length of the scale bar
    width_m <- max(ext[2] - ext[1])
    scale_length <- if (width_m > 1000) 100 else 50
    sbar(d = scale_length,  # length in meters
         xy = c("bottomleft"),
         divs = 2,  # number of divisions
         labels = paste(scale_length,"m"),  # labels for start and end
         style = "bar",
         lwd = 2,
         cex = 0.8)
  }
}

par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))

# Plot Waldi sampling scheme
plot_sampling_scheme(wal_sam_buffer, "Waldi Sampling Scheme")
points(trees_v,col="red", cex = 0.4)

# Plot Allenwiller sampling scheme
plot_sampling_scheme(al_sam_buffer, "Allenwiller Sampling Scheme")
points(trees1_v,col="red", cex = 0.4)

#### grid for simulation ####

# ALLENWILLER #

# load data and sampling scheme
obs <- read.csv("Newhybrids/Allenwiller_indiv_info_class.csv")
obs_v <- vect(obs,geom = c("x", "y"), crs ="EPSG:3035" )
plot(obs_v)
sampling_v <- vect("Allenwiller_sampling_scheme.shp")
plot(sampling_v, add = T)

# create grid over the extent of the data + add some meters to create a buffer area to absorbe dispersal
grid_extent <- ext(obs_v)
grid_extent <- grid_extent + 40

# set the resolution 
resolution = 4

## add a margin of few meters to include all the trees
grid <- rast(ext = grid_extent, resolution = resolution, crs = "EPSG:3035")
ncell(grid)
dim(grid)
# assign sequential IDs to the raster cells (starting from 1)
values(grid) <- 1:ncell(grid)
names(grid) <-  "patch.ID"
# convert the raster to polygons, retaining the cell IDs
grid_v <- as.polygons(grid, dissolve = FALSE)
print(head(grid_v$patch.ID))

# visual check 
plot(grid_v)
plot(obs_v, add = T)
plot(sampling_v, add = T)
sbar(d = 50)

writeVector(grid_v, paste0("Allenwiller_grid_r", resolution, "m_", ncell(grid), "patches.shp"), overwrite = TRUE)

# check if patch.ID is stored
final_grid <- vect( paste0("Allenwiller_grid_r", resolution, "m_", ncell(grid), "patches.shp"))
print(head(final_grid$patch.ID))
plot(final_grid)
text(final_grid, "patch.ID", cex=0.4)

# WALDI #

obs <- read.csv("Newhybrids/Waldi_indiv_info_class.csv")
obs_v <- vect(obs,geom = c("x", "y"), crs ="EPSG:3035" )
plot(obs_v)
sampling_v <- vect("Waldi_sampling_scheme.shp")
plot(sampling_v, add = T)
# create grid over the extent of the data -add some meters to cover all the points
grid_extent <- ext(obs_v)
grid_extent <- grid_extent + 40

# set the resolution 
resolution = 4

## add a margin of few meters to include all the trees
grid <- rast(ext = grid_extent, resolution = resolution, crs = "EPSG:3035")
ncell(grid)
dim(grid)
# assign sequential IDs to the raster cells (starting from 1)
values(grid) <- 1:ncell(grid)
names(grid) <-  "patch.ID"
# convert the raster to polygons, retaining the cell IDs
grid_v <- as.polygons(grid, dissolve = FALSE)
print(head(grid_v$patch.ID))

# visual check 
plot(grid_v)
plot(obs_v, add = T)
plot(sampling_v, add = T)
sbar(d = 50)

writeVector(grid_v, paste0("Waldi_grid_", resolution, "m_", ncell(grid), "patches.shp"), overwrite = TRUE)

# check if cellID is stored
final_grid <- vect( paste0("Waldi_grid_",resolution, "m_", ncell(grid), "patches.shp"))
print(head(final_grid$patch.ID))
plot(final_grid)
text(final_grid, "patch.ID", cex=0.4)

#### simulate individuals  #### 

# for each stand, create a locus list with allele frequencies
# NB: in theory, we need just one locus list (it should be the same for the 2 stands), however, there are few orientalis in Waldi and we want to reproduce the same allele frequencies as the reality

# ALLENWILLER #

# load individuals infos
obs <- read.csv("Newhybrids/Allenwiller_indiv_info_class.csv")

# read genepop file (but txt), read header and loci
lines <- readLines("Newhybrids/Allenwiller_genepop_all.txt")
loci_end <- which(lines == "Pop")[1]
header <- lines[1]
loci <- lines[2:(loci_end-1)]
# get genotypes (after ("pop"))
genotypes <- lines[(loci_end+1):length(lines)]
# extract individuals
ind_ids <- sapply(strsplit(genotypes, ","), function(x) trim(x[1]))

# in the genotype file, take only the adults (from the observed dataframe) and divide per subspecies
adults_ori <- obs %>% filter(generation == "Adult",NH_class_t80 == "Pure orientalis") %>%pull(SampleID)
adults_syl <- obs %>% filter(generation == "Adult",NH_class_t80 == "Pure sylvatica") %>%pull(SampleID)
ori_data <- genotypes[ind_ids %in% adults_ori]
syl_data <- genotypes[ind_ids %in% adults_syl]
# combine the adults into 2 populations, and write the genind file
genind <- c(header,loci,"Pop",syl_data,"Pop", ori_data)
#writeLines(genind, "Newhybrids/Allenwiller_genind_adults.gen")

genind <- read.genepop("Newhybrids/Allenwiller_genind_adults.gen", ncode = 3L, quiet = TRUE)
popNames(genind)

# convert to genepop and compute allele frequecies
genepop <- genind2genpop(genind)
allele_freq <- as.data.frame(makefreq(genepop))
rownames(allele_freq)[1]= "Sylvatica"
rownames(allele_freq)[2]= "Orientalis"

# dataframe
allele_freq_df = allele_freq %>% tibble::rownames_to_column("Species") %>%
  pivot_longer(cols =-Species, names_to = c("Locus", "Alleles"), names_sep = "\\.", values_to = "Freq")

# visual check 
mycols <- c("Orientalis" = "#482173FF", "Sylvatica" = "#ffcc00")
colSet <- scale_colour_manual(values = mycols)
fillSet <- scale_fill_manual(values = mycols)

ggplot(allele_freq_df, aes(x = Alleles, y = Freq, fill = Species)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  facet_wrap(~ Locus, scales = "free_x") +  
  fillSet+
  labs(
    x = "Alleles",
    y = "Frequency"
  ) +
  theme_minimal()

# creating a list of dataframes where each dataframe is a locus (Locus - Allele - Allelic frequency)
locus_list_a <-  split(allele_freq_df, allele_freq_df$Locus)
head(locus_list_a)[1]

# WALDI #

# load the individuals info
obs <- read.csv("Newhybrids/Waldi_indiv_info_class.csv")

# read genepop file (but txt), read header and loci
lines <- readLines("Newhybrids/Waldi_genepop_all.txt")
loci_end <- which(lines == "Pop")[1]
header <- lines[1]
loci <- lines[2:(loci_end-1)]
# get genotypes (after ("pop"))
genotypes <- lines[(loci_end+1):length(lines)]
# extract individuals
ind_ids <- sapply(strsplit(genotypes, ","), function(x) trim(x[1]))

# in the genotype file, take only the adults (from the observed dataframe) and divide per subspecies
adults_ori <- obs %>% filter(generation == "Adult",NH_class_t80_corrected == "Pure orientalis") %>%pull(SampleID)
adults_syl <- obs %>% filter(generation == "Adult",NH_class_t80_corrected == "Pure sylvatica") %>%pull(SampleID)
ori_data <- genotypes[ind_ids %in% adults_ori]
syl_data <- genotypes[ind_ids %in% adults_syl]
# combine the adults into 2 populations, and write the genind file
genind <- c(header,loci,"Pop",syl_data,"Pop", ori_data)
#writeLines(genind, "Newhybrids/Waldi_genind_adults.gen")

genind <- read.genepop("Newhybrids/Waldi_genind_adults.gen", ncode = 3L, quiet = TRUE)
popNames(genind)

# convert to genepop and compute allele frequecies
genepop <- genind2genpop(genind)
allele_freq <- as.data.frame(makefreq(genepop))
rownames(allele_freq)[1]= "Sylvatica"
rownames(allele_freq)[2]= "Orientalis"

# dataframe
allele_freq_df = allele_freq %>% tibble::rownames_to_column("Species") %>%
  pivot_longer(cols =-Species, names_to = c("Locus", "Alleles"), names_sep = "\\.", values_to = "Freq")

# visual check 
mycols <- c("Orientalis" = "#482173FF", "Sylvatica" = "#ffcc00")
colSet <- scale_colour_manual(values = mycols)
fillSet <- scale_fill_manual(values = mycols)

ggplot(allele_freq_df, aes(x = Alleles, y = Freq, fill = Species)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  facet_wrap(~ Locus, scales = "free_x") +  
  fillSet+
  labs(
    x = "Alleles",
    y = "Frequency"
  ) +
  theme_minimal()

# creating a list of dataframes where each dataframe is a locus (Locus - Allele - Allelic frequency)
locus_list_w <-  split(allele_freq_df, allele_freq_df$Locus)
head(locus_list_w)[1]



###### function to simulate an individual by sampling 2 alleles for each locus weighted by allele freuqnencis, for each spacies
simul_indiv <- function(locus_list, species) {
  # create an individual df where each row has the locus name
  indiv <- data.frame(matrix(ncol = 2, nrow = length(locus_list)))      
  colnames(indiv) <- c("Allele1", "Allele2")
  row.names(indiv) <- names(locus_list)          
  # for each locus sample two random alleles with probability = allele frequency
  for (locus_name in names(locus_list)) {
    locus_data <- locus_list[[locus_name]]  # take data for specific locus
    # take only data for the chosen species
    species_data <- subset(locus_data, Species == species)
    # take two alleles with probability == frequnency
    indiv[locus_name, 1:2] <- sample(
      species_data$Alleles, 
      size = 2, 
      replace = TRUE, 
      prob = species_data$Freq
    )
  }
  
  # merge alleles into genotype column
  indiv$Genotype <- paste0(indiv$Allele1, indiv$Allele2)
  return(indiv)
}


# function to create a population of a defined size for a given species
make_pop <- function(pop.size, locus_list, species) {
  # use the simul_indiv function as many times as the population size
  individuals <- lapply(seq_len(pop.size), function(x) simul_indiv(locus_list, species))
  # extract the genotype column for all individuals
  genotypes <- lapply(individuals, function(indiv) indiv$Genotype)
  # combnie genotypes into a single dataframe
  pop <- do.call(cbind, genotypes) # create array
  pop <- data.frame(pop)    # convert to df where each row is a locus       

  row.names(pop) <- names(locus_list)
  # name each individuals with the first 3 letters of its species
  colnames(pop) <- paste0(substr(tolower(species), 1, 3), "_sim_", seq_len(ncol(pop)))
  
  return(pop)
}

 
# test creating individulas, specify locus list and species
simul_indiv(locus_list_a, species = "Sylvatica")
simul_indiv(locus_list_w, species = "Sylvatica")
simul_indiv(locus_list_a, species = "Orientalis")
simul_indiv(locus_list_w, species = "Orientalis")
# test creating populations
ori_sim_10 <- make_pop(10, locus_list_a, species = "Orientalis")
ori_sim_10

# ----- calculate the proportions of each species taking the real situation

# ALLENWILLER #

# calculate proportions of orientalis / sylvatica adults in observed data
obs <- read.csv("Newhybrids/Allenwiller_indiv_info_class.csv")

prop <- obs %>%
  group_by(generation) %>%
  mutate(tot = n()) %>%
  group_by(generation, NH_class_t80) %>%
  summarise(class_n = n(), class_prop = n()/first(tot), .groups = 'drop')
prop

# load grid and plot observed individuals
grid <- vect("Allenwiller_grid_r4m_9844patches.shp")
obs_v <- vect(obs, geom = c("x", "y"), crs = "EPSG:3035")

plot(grid)
points(obs_v, col = as.factor(obs$NH_class_t80), cex =  obs$Age*0.01 )

#create polygons for simulation
polygon_ori <- rbind(
  c(4126715, 2840045),
  c(4126800, 2840120),
  c(4126705, 2840152),
  c(4126715, 2840045)  # Closing the polygon
)
polygon_ori<- vect(polygon_ori, type = "polygons", crs = "EPSG:3035")
plot(polygon_ori, add = T)

polygon_syl <- rbind(
  c(4126800, 2840120),
  c(4126870, 2840200),
  c(4126705, 2840260),
  c(4126705, 2840152), 
  c(4126800, 2840120) # Closing the polygon
)
polygon_syl<- vect(polygon_syl, type = "polygons", crs = "EPSG:3035")
plot(polygon_syl, add = T)

# count how many adults individuals are in each polygon, for each class and check the circumference
ori_points <- relate(obs_v, polygon_ori, "intersects")
ori_in_polygon <- obs[ori_points, ]
table(ori_in_polygon$generation, ori_in_polygon$NH_class_t80)
ori_in_polygon <- subset(ori_in_polygon, generation == "Adult" & NH_class_t80 == "Pure orientalis")

# sylvatica 
syl_points <- relate(obs_v, polygon_syl, "intersects")
syl_in_polygon <- obs[syl_points, ]
table(syl_in_polygon$generation, syl_in_polygon$NH_class_t80)
syl_in_polygon <- subset(syl_in_polygon, generation == "Adult" & NH_class_t80 == "Pure sylvatica")

# how many adult trees today in the polygons
n_ori <- nrow(ori_in_polygon)
n_ori
n_syl <- nrow(syl_in_polygon)
n_syl

# proportion of orientalis / sylvatica
n_ori/n_syl

# check how many adult trees today in the polygons with more than 60 years
nrow(subset(syl_in_polygon, Age >  60 ))
nrow(subset(ori_in_polygon, Age >  60 ))

# -------- simulating individuals 
syl_points=centroids(grid)
precision=13.23
# keep only the points that fall into the polydon
syl_in_poly <- relate(syl_points, polygon_syl, "intersects")
syl_simulated <- syl_points[syl_in_poly]
obj=as.data.frame(syl_simulated,geom="XY")
obj$X2=plyr::round_any(obj$x,precision)
obj$Y2=plyr::round_any(obj$y,precision)
obj2=obj %>% group_by(X2,Y2) %>% summarise(x=x[which.min(abs(x-median(x)))],y=y[which.min(abs(y-median(y)))])
syl_simulated=vect(obj2,geom=c("x","y"))
plot(grid)
plot(polygon_syl, add = TRUE)
plot(syl_simulated, add = TRUE, col = "gold")
nrow(syl_simulated)

ori_points=centroids(grid)
precision=10.96
# keep only the points that fall into the polydon
ori_in_poly <- relate(ori_points, polygon_ori, "intersects")
ori_simulated <- ori_points[ori_in_poly]
obj=as.data.frame(ori_simulated,geom="XY")
obj$X2=plyr::round_any(obj$x,precision)
obj$Y2=plyr::round_any(obj$y,precision)
obj2=obj %>% group_by(X2,Y2) %>% summarise(x=x[which.min(abs(x-median(x)))],y=y[which.min(abs(y-median(y)))])
ori_simulated=vect(obj2,geom=c("x","y"))
plot(grid)
plot(polygon_ori, add = TRUE)
plot(ori_simulated, add = TRUE, col = "purple")
nrow(ori_simulated)

# compare the observed with the simulated numbers
n_ori
dim(ori_simulated)[1]
n_syl
dim(syl_simulated)[1]

# if it is good, save the coordinates and the patch ID
ori_coords <- crds(ori_simulated)
ori_patches <- terra::extract(grid, ori_simulated)

syl_coords <- crds(syl_simulated)
syl_patches <- terra::extract(grid, syl_simulated)

# simulate individuals with functions make.pop
# load locus list for Allenwiller
ori_sim <- make_pop(52, locus_list_a, species = "Orientalis")
syl_sim <-  make_pop(97, locus_list_a, species = "Sylvatica")

# assign coordinates to the simulated individuals
ori_sim2 <- cbind(t(ori_sim), ori_coords, ori_patches$patch.ID)
syl_sim2 <- cbind(t(syl_sim), syl_coords, syl_patches$patch.ID)

# save the simulated population
sim_pop <- as.data.frame(rbind(ori_sim2, syl_sim2))
sim_pop <- tibble::rownames_to_column(sim_pop, "SampleID")
colnames(sim_pop)[20] = "patch_id"
#write.table(sim_pop, "NEMO/Allenwiller_simulated_spatial.txt", row.names = F,col.names = T, quote = FALSE)

# WALDI #

# calculate proportions of orientalis / sylvatica adults in observed data
obs <- read.csv("Newhybrids/Waldi_indiv_info_class.csv")

prop <- obs %>%
  group_by(generation) %>%
  mutate(tot = n()) %>%
  group_by(generation, NH_class_t80_corrected) %>%
  summarise(class_n = n(), class_prop = n()/first(tot), .groups = 'drop')
prop

# load grid and plot observed individuals
grid <- vect("Waldi_grid_4m_5396patches.shp")
obs_v <- vect(obs, geom = c("x", "y"), crs = "EPSG:3035")

plot(grid)
points(obs_v, col = as.factor(obs$NH_class_t80), cex =  obs$Age*0.04 )

#create polygon for orientalis 
polygon_ori <- rbind(
  c(4253611.693750093,	2724337.9876462617	),
  c(4253676.985703846,	2724374.137993538),
  c(4253700.778534452,	2724321.5724375495),
  c(4253622.575742385,	2724286.7131741047) , 
  c(4253611.140428451,	2724337.9876462617	), 
  c(4253611.693750093,	2724337.9876462617)# Closing the polygon
)

polygon_ori<- vect(polygon_ori, type = "polygons", crs = "EPSG:3035")
plot(polygon_ori, add = T)

# count how many adults individuals are in each polygon, for each class and check the circumference
ori_points <- relate(obs_v, polygon_ori, "intersects")
ori_in_polygon <- obs[ori_points, ]
table(ori_in_polygon$generation, ori_in_polygon$NH_class_t80)
ori_in_polygon <- subset(ori_in_polygon, generation == "Adult" & NH_class_t80 == "Pure orientalis")

polygon_syl <- rbind(
  c(4253611.048208177,	2724342.229778851	),
  c(4253666.564812922,	2724397.026103004	),
  c(4253766.716030122,	2724365.6537283612	),
  c(4253835	,2724229.1589825405	),
  c(4253679.475651235,2724156.791	),
  c(4253623.2212843,	2724285.606530821	),
  c(4253701.424076367	,2724322.3101997394	),
  c(4253677.4468052145,	2724374.1379935383),	
  c( 4253611.601529819,	2724342.96754104),	
  c(4253611.048208177,	2724342.229778851)
)
polygon_syl<- vect(polygon_syl, type = "polygons", crs = "EPSG:3035")
plot(polygon_syl, add = T)

# sylvatica 
syl_points <- relate(obs_v, polygon_syl, "intersects")
syl_in_polygon <- obs[syl_points, ]
table(syl_in_polygon$generation, syl_in_polygon$NH_class_t80)
syl_in_polygon <- subset(syl_in_polygon, generation == "Adult" & NH_class_t80 == "Pure sylvatica")

# how many adult trees today in the polygons
n_ori <- nrow(ori_in_polygon)
n_ori
n_syl <- nrow(syl_in_polygon)
n_syl

# proportion of orientalis / sylvatica
n_ori/n_syl

# check how many adult trees today in the polygons with more than 60 years
nrow(subset(syl_in_polygon, Age >  60 ))
nrow(subset(ori_in_polygon, Age >  60 ))

# --------- simulate spatial scenario

syl_points=centroids(grid)
precision=13.238
# keep only the points that fall into the polydon
syl_in_poly <- relate(syl_points, polygon_syl, "intersects")
syl_simulated <- syl_points[syl_in_poly]
obj=as.data.frame(syl_simulated,geom="XY")
obj$X2=plyr::round_any(obj$x,precision)
obj$Y2=plyr::round_any(obj$y,precision)
obj2=obj %>% group_by(X2,Y2) %>% summarise(x=x[which.min(abs(x-median(x)))],y=y[which.min(abs(y-median(y)))])
syl_simulated=vect(obj2,geom=c("x","y"))
plot(grid)
plot(polygon_syl, add = TRUE)
plot(syl_simulated, add = TRUE, col = "gold")
nrow(syl_simulated)

ori_points=centroids(grid)
precision=35
# keep only the points that fall into the polydon
ori_in_poly <- relate(ori_points, polygon_ori, "intersects")
ori_simulated <- ori_points[ori_in_poly]
obj=as.data.frame(ori_simulated,geom="XY")
obj$X2=plyr::round_any(obj$x,precision)
obj$Y2=plyr::round_any(obj$y,precision)
obj2=obj %>% group_by(X2,Y2) %>% summarise(x=x[which.min(abs(x-median(x)))],y=y[which.min(abs(y-median(y)))])
ori_simulated=vect(obj2,geom=c("x","y"))
plot(grid)
plot(polygon_ori, add = TRUE)
plot(ori_simulated, add = TRUE, col = "purple")
nrow(ori_simulated)

# compare the observed with the simulated numbers
n_ori
dim(ori_simulated)[1]
n_syl
dim(syl_simulated)[1]

# if it is good, save the coordinates and the patch ID
ori_coords <- crds(ori_simulated)
ori_patches <- terra::extract(grid, ori_simulated)

syl_coords <- crds(syl_simulated)
syl_patches <- terra::extract(grid, syl_simulated)

# simulate individuals with functions make.pop

# load locus list for Waldi
ori_sim <- make_pop(8, locus_list_w, species = "Orientalis")
syl_sim <-  make_pop(193, locus_list_w, species = "Sylvatica")

# assign coordinates to the simulated individuals
ori_sim2 <- cbind(t(ori_sim), ori_coords, ori_patches$patch.ID)
syl_sim2 <- cbind(t(syl_sim), syl_coords, syl_patches$patch.ID)

# save the simulated population
sim_pop <- as.data.frame(rbind(ori_sim2, syl_sim2))
sim_pop <- tibble::rownames_to_column(sim_pop, "SampleID")
colnames(sim_pop)[20] = "patch_id"
#write.table(sim_pop, "NEMO/Waldi_simulated_spatial.txt", row.names = F,col.names = T, quote = FALSE)


#### FSTAT files ####

# individuals_file must contain SampleID, Genotype (loci coded like 123123), x,y, coordinates and patch_id in the simulation grid

build_fstat_file <- function(path, individuals_file, grid, output_file, age_range = age_range, stage = stage, sex = 1, pedigree = 0) {
  
  # read simulated/real population
  fstat_df <- read_table(individuals_file, col_names = T)
  # create FSTAT extended version
  # remove the coordintes
  fstat_ext <- subset(fstat_df, select = -c(x, y))
  # set other columns
  fstat_ext$Age <- sample(age_range[1]:age_range[2], nrow(fstat_ext), replace = TRUE)
  fstat_ext$Stage <- stage
  fstat_ext$Sex <- sex
  fstat_ext$Pedigree <- pedigree
  # reorder columns
  fstat_ext <- fstat_ext %>% 
    relocate(Stage, .after = sfc_1143) %>% 
    relocate(Age, .after = Stage) %>% 
    relocate(Sex, .after = Age) %>% 
    relocate(Pedigree, .after = Sex) %>% 
    relocate(patch_id, .after = Pedigree)
  # remove first column (individual names)
  fstat_ext_final <- fstat_ext[, -1]
  
  # replace '000000' with '001001' (ONLY FOR TRUE INDIVIDUALS; SIMULATED DONT HAVE MISSING INFORMATION)
  fstat_ext_final <- apply(fstat_ext_final, 2, function(x) gsub("000000", "001001", x))
  fstat_ext_final <- as.data.frame(fstat_ext_final)
  
  # set first column as patch_id
  fstat_ext_final = cbind(patch_id = fstat_ext_final[, ncol(fstat_ext_final)], fstat_ext_final)
  
  # add header (n patches, col number - 1, max number of alleles per locus, number of digits for each genotype)
  header <- c(
    paste(ncell(grid),ncol(fstat_ext_final)-1, 238, 3),
    "casolfagus_29",
    "concat14",
    "csolfagus_05",
    "csolfagus_06",
    "csolfagus_19",
    "csolfagus_31",
    "DE576",
    "DUKCT",
    "DZ447",
    "EEU75",
    "EJV8T",
    "EMILY",
    "ERHBI",
    "FS1_15",
    "sfc_0036",
    "sfc_1143",
    "stage",
    "age",
    "sex",
    "ped",
    "origin"
  )

  # write header and after the data
  writeLines(header, output_file)
  write.table(fstat_ext_final, output_file, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  return(fstat_ext_final)
}

# creating FSTAT files for both spatial and non spatial scenario
build_fstat_file(
  age_range = c(50:51), #all the trees are adults
  stage = 3,
  grid = vect("Allenwiller_grid_r4m_9844patches.shp"),
  individuals_file =  "NEMO/Allenwiller_simulated_spatial.txt",
  output_file = "NEMO/Allenwiller_FSTAT_simulated_spatial.txt"
)

build_fstat_file(
  age_range = c(50:51), #all the trees are adults
  stage = 3,
  grid = vect("Waldi_grid_4m_5396patches.shp"),
  individuals_file =  "NEMO/Waldi_simulated_spatial.txt",
  output_file = "NEMO/Waldi_FSTAT_simulated_spatial.txt"
)



# files with adult from age 30
build_fstat_file(
  age_range = c(30:31), #all the trees are adults
  stage = 3,
  grid = vect("Allenwiller_grid_r4m_9844patches.shp"),
  individuals_file =  "NEMO/Allenwiller_simulated_spatial.txt",
  output_file = "NEMO/Allenwiller_FSTAT_simulated_spatial_30.txt"
)

build_fstat_file(
  age_range = c(30:31), #all the trees are adults
  stage = 3,
  grid = vect("Waldi_grid_4m_5396patches.shp"),
  individuals_file =  "NEMO/Waldi_simulated_spatial.txt",
  output_file = "NEMO/Waldi_FSTAT_simulated_spatial_30.txt"
)

#### quanti_init_freq matrices ####

# matrix defining the species of the individuals taken from the FSTAT file
# 1 column, nrows = n patches of the simulation grid
# each patch is 0 if contains sylvatica, 1 if contains orientalis

# function to write matrix in NEMO format (change such that NA result in emtpy rows)
write.matrix.nemo <- function(mat, outfile) {
  rows <- dim(mat)[1]
  cols <- dim(mat)[2]
  cat("{", file = outfile)
  
  for(i in 1:rows) {
    cat("{", file = outfile, append = TRUE)
    cat(mat[i,], sep = ",", file = outfile, append = TRUE)
    cat("}\n", file = outfile, append = TRUE)
  }
  
  cat("}\n", file = outfile, append = TRUE)
}

# ALLENWILLER #

# importing simulation grid 
grid <- vect("Allenwiller_grid_r4m_9844patches.shp")

# --------- spatial scenario 

# importing simulated population with species information and patch ID (check the path)
indiv <- fread("NEMO/Allenwiller_simulated_spatial.txt")
# create a species column (1 orientalis, 0 sylvatica)
indiv$species_code = ifelse(grepl("^ori_sim", indiv$SampleID), 1, 0)

# creating matrix with 0 as default value (empty patches will be ignored automatically, but they need to appear in the matrix)
quanti_init_freq = matrix(0, ncol = 1, nrow = ncell(grid))
# assign patch value based on patch id in individuals file
quanti_init_freq[ c(indiv$patch_id), ] = indiv$species_code

# write quanti_init matrix to file (or whether you want to store them) and copy to the nemo folder 
write.matrix.nemo(quanti_init_freq, paste0(input_path,"/NEMO_quanti_init_freq_a_spatial.txt"))

# ----------- non spatial scenario

# importing simulated population with species information and patch ID (check the path)
indiv <- fread(paste0(input_path, "/Allenwiller_simulated_nonspatial.txt"))
# create a species column (1 orientalis, 0 sylvatica)
indiv$species_code = ifelse(grepl("^ori_sim", indiv$SampleID), 1, 0)


# creating matrix with 0 as default value (empty patches will be ignored automatically, but they need to appear in the matrix)
quanti_init_freq = matrix(0, ncol = 1, nrow = ncell(grid))
# assign patch value based on patch id in individuals file
quanti_init_freq[ c(indiv$patch_id), ] = indiv$species_code


# write quanti_init matrix to file (or whether you want to store them) and copy to the nemo folder 
outfile_quanti_init <- paste0(input_path,"/NEMO_quanti_init_freq_a_nonspatial.txt")
write.matrix.nemo(quanti_init_freq, outfile_quanti_init)
file.copy(from=paste0(input_path,"/NEMO_quanti_init_freq_a_nonspatial.txt"),
          to=paste0(exe_path, "NEMO_quanti_init_freq_a_nonspatial.txt"), overwrite = T)


# WALDI #

# importing simulation grid 
grid <- vect("Waldi_grid_4m_5396patches.shp")

# ------------- spatial scenario 

# importing simulated population with species information and patch ID
indiv <- fread("NEMO/Waldi_simulated_spatial.txt")
# create a species column (1 orientalis, 0 sylvatica)
indiv$species_code = ifelse(grepl("^ori_sim", indiv$SampleID), 1, 0)

# creating matrix with 0 as default value (empty patches will be ignored automatically, but they need to appear in the matrix)
quanti_init_freq = matrix(0, ncol = 1, nrow = ncell(grid))
# assign patch value based on patch id in individuals file
quanti_init_freq[ c(indiv$patch_id), ] = indiv$species_code

# write quanti_init matrix to file (or whether you want to store them) and copy to the nemo folder 

write.matrix.nemo(quanti_init_freq, paste0(input_path,"/NEMO_quanti_init_freq_w_spatial.txt"))


# ------------ non spatial scenario 

# importing simulated population with species information and patch ID
indiv <- fread(paste0(input_path, "/Waldi_simulated_nonspatial.txt"))
# create a species column (1 orientalis, 0 sylvatica)
indiv$species_code = ifelse(grepl("^ori_sim", indiv$SampleID), 1, 0)

# creating matrix with 0 as default value (empty patches will be ignored automatically, but they need to appear in the matrix)
quanti_init_freq = matrix(0, ncol = 1, nrow = ncell(grid))
# assign patch value based on patch id in individuals file
quanti_init_freq[ c(indiv$patch_id), ] = indiv$species_code

# write quanti_init matrix to file (or whether you want to store them) and copy to the nemo folder 
outfile_quanti_init <- paste0(input_path,"/NEMO_quanti_init_freq_w_nonspatial.txt")
write.matrix.nemo(quanti_init_freq, outfile_quanti_init)
file.copy(from=paste0(input_path,"/NEMO_quanti_init_freq_w_nonspatial.txt"),
          to=paste0(exe_path, "NEMO_quanti_init_freq_w_nonspatial.txt"), overwrite = T)


#### carrying_capacity matrices (patch_nbfem) ####

# creating matrices where patches without beech have carrying capacity 0, and patches with beech can have a defined carrying capacity

# Function to write matrix in NEMO format
write.matrix.nemo <- function(mat, outfile) {
  rows <- dim(mat)[1]
  cols <- dim(mat)[2]
  cat("{", file = outfile)
  
  for(i in 1:rows) {
    cat("{", file = outfile, append = TRUE)
    cat(mat[i,], sep = ",", file = outfile, append = TRUE)
    cat("}\n", file = outfile, append = TRUE)
  }
  
  cat("}\n", file = outfile, append = TRUE)
}

# ALLENWILLER #

# loading species classification (from remote sensing, 0 = beech) and project to metric system
class <- rast("~/Allen_classification.tif")
class <- project(x=class,y= "epsg:3035")
plot(class)
# load simulation grid
grid <- vect("Allenwiller_grid_r4m_9844patches.shp",crs=crs("+init=epsg:3035"))

# for each cell of the grid, extract the value of classification in a dataframe (min to be conservative)
class_extract = terra::extract(class, grid, fun = "min", cells = T, xy = T, layer = 1)
nrow(class_extract)

# (OPTIONAL) check the highest number of seedlings counted into circle plots
seed_n <- read.csv("Allenwiller_circleplot_info.csv")
max(seed_n$CountedC1.new.)

# define vector of carrying capacity values to test
car_cap = c(200,500,800)

# set the carrying capacity to 0 on the patches at the borders, to avoid edge effects: 
grid_ext <- ext(grid)
coords <- terra::centroids(grid) |> crds()
border <- 20
border_patch <- coords[,1] <= (grid_ext[1] + border) | 
  coords[,1] >= (grid_ext[2] - border) |
  coords[,2] <= (grid_ext[3] + border) | 
  coords[,2] >= (grid_ext[4] - border)

plot(grid, col = border_patch)

# for each value to test, create a vector of max carrying capacity observed in patches with presence of beech (0), and 0 where there is not beech
for(i in seq_along(car_cap)){
  # create column for each value to test
  class_extract[paste0("car_cap_", car_cap[i])] <- ifelse(class_extract$value <= 0.5, car_cap[i],0)
  # set border patches to 0
  class_extract[border_patch, paste0("car_cap_", car_cap[i])] <- 0
  # assign to the grid the values for visual check
  grid$car_cap <- class_extract[paste0("car_cap_", car_cap[i])]
  # vector to matrix
  car_cap_matrix <- matrix(grid$car_cap, 
                           nrow = 1, 
                           ncol = ncell(grid$patch.ID))
  
  # Write carrying capacity matrix to file
  write.matrix.nemo(car_cap_matrix, paste0(input_path, "/NEMO_patch_nbfem_a_spatial_", car_cap[i], ".txt"))
  
  }

# visual check
plot(grid, col = grid$car_cap)

# WALDI #
# importing NIR data for Waldi (fro remote sensing data) to set carrying cap only in the stand
ndvi_0317 = rast("~//Data/NDVI_classification/Waldi_2017_2019/Waldi_2017/PSScene/20170316_093506_0e1f_3B_AnalyticMS_SR_clip.tif")
ndvi_0317 <- project(x=ndvi_0317,y= "epsg:3035")
plot(ndvi_0317$nir)

# import simulation grid
grid <- vect("Waldi_grid_4m_5396patches.shp",crs=crs("+init=epsg:3035"))

# for each cell of the grid, extract the value of classification in a dataframe (min to be conservative)
ndvi_extract = terra::extract(ndvi_0317, grid, fun = "mean", cells = T, xy = T, layer = "nir")
nrow(ndvi_extract)

# (OPTIONAL) check the highest number of seedlings counted into circle plots
seed_n <- read.csv("Waldi_circleplot_info.csv")
max(seed_n$NClass1Class3)

# define vector of carrying capacity values to test
car_cap = c(50,200,500)

# set the carrying capacity to 0 on the patches at the borders, to avoid edge effects: 
grid_ext <- ext(grid)
coords <- terra::centroids(grid) |> crds()
border <- 20
border_patch <- coords[,1] <= (grid_ext[1] + border) | 
  coords[,1] >= (grid_ext[2] - border) |
  coords[,2] <= (grid_ext[3] + border) | 
  coords[,2] >= (grid_ext[4] - border)

plot(grid, col = border_patch)

# for each value to test, create a vector of max carrying capacity observed in cells with presence of beech (0), and 0 where there is not beech
for(i in seq_along(car_cap)){
  # create column for each value to test
  ndvi_extract[paste0("car_cap_", car_cap[i])] <- ifelse(ndvi_extract$value <= 2000, car_cap[i],0)
  # set border patches to 0
  ndvi_extract[border_patch, paste0("car_cap_", car_cap[i])] <- 0
  # assign to the grid the values for visual check
  grid$car_cap <- ndvi_extract[paste0("car_cap_", car_cap[i])]
  # vector to matrix
  car_cap_matrix <- matrix(grid$car_cap, 
                           nrow = 1, 
                           ncol = ncell(grid$patch.ID))
  
  # Write carrying capacity matrix to file
  #write.matrix.nemo(car_cap_matrix, paste0(input_path, "/NEMO_patch_nbfem_w_spatial_", car_cap[i], ".txt"))

}

# visual check
plot(grid, col = grid$car_cap)

# check how many patches can be colonized
length(grid[grid$car_cap != 0])

#### seed dispersal and pollen dispersal matrices ####

# dispersal function (implement such that can be customized)

exp_power_kernel <- function(d, b, x) {
  #  scale parameter a
  a <- d * exp(lgamma(2 / b) - lgamma(3 / b))
  
  #  kernel function K_jk
  return(b * exp(-(x^b / a^b))  / (2* pi * a^2 * exp(lgamma(2/b))))
}

# function to build both reduced dispersal matrices and connectivity matrix
# num_patch = number of patch of the simulation grid
# distance_matrix = pairwise distance matrix between the patches
# d = mean dispersal distance for the kernel
# b = shape of the kernel
# d_thresh = minimum dispersal rate under which dispersal rate is set to 0

build.reduced.dispersal.matrices <- function(num_patch, distance_matrix,  d, b, d_thresh = threshold) {
  stopifnot(dim(distance_matrix)[1] == num_patch)
  
  dispersal_matrix <- matrix(0, nrow = num_patch, ncol = num_patch)
  rate_matrix <- matrix(0, nrow = num_patch, ncol = num_patch)
  connectivity_matrix <- matrix(NA, nrow = num_patch, ncol = num_patch)
  
  for (i in 1:num_patch) {
    for (j in 1:num_patch) {
      # calculate dispersal probability from patch i to patch j
      rate_matrix[i, j] <- exp_power_kernel(d, b, distance_matrix[i, j])
    }
  }
  # define full disp matrix 
  dispersal_matrix <- rate_matrix
  
  for (i in 1:num_patch) {
    # Get indices of the patches sorted by dispersal probability in descending order
    ord <- order(rate_matrix[i, ], decreasing = TRUE)
    
    # re-order the values in row i of reduced rate matrix, connectivity matrix
    rate_matrix[i, ] <- rate_matrix[i, ord]
    connectivity_matrix[i, ] <- ord

    # remove elements below the threshold by setting them equal to NA 
    to_remove <- which(rate_matrix[i, ] < d_thresh)
    rate_matrix[i, to_remove] <- NA
    connectivity_matrix[i, to_remove] <- NA
    
    # normalize dispersal rates to sum to 1 (NB: THIS IS CHANGING SLIGHTLY THE VALUES ACROSS ROWS
    sum_d <- sum(rate_matrix[i, ], na.rm = TRUE)
    rate_matrix[i, !is.na(rate_matrix[i, ])] <- rate_matrix[i, !is.na(rate_matrix[i, ])] / sum_d
  }
  
  return(list(dispersal_matrix = dispersal_matrix, connectivity_matrix = connectivity_matrix, rate_matrix = rate_matrix))
}

# function to write (general) matrix in NEMO format 
write.matrix.nemo <- function(mat, outfile) {
  rows <- dim(mat)[1]
  cols <- dim(mat)[2]
  cat("{", file = outfile)
  
  for (i in 1:rows) {
    cat("{", file = outfile, append = TRUE)
    cat(mat[i, ], sep = ",", file = outfile, append = TRUE)
    cat("}\n", file = outfile, append = TRUE)
  }
  
  cat("}\n", file = outfile, append = TRUE)
}

## original function for the dipersal matrix 
# write a reduced dispersal matrix, removing elements flagged with NA
# works for both the connectivity and dispersal rate matrices
# mat        : the input matrix (connectivity or rate)
# outfile    : name of the output text file
# discard.NA : if NA values are discarded or not, if not, calls 'write.matrix.nemo'

write.dispersal.matrix = function(mat, outfile, discard.NA=TRUE) {
  
  if(discard.NA) {
    
    cat("{",file=outfile)
    
    # write every row at once, removing NA elements
    for(i in 1:dim(mat)[1]) {
      cat("{",file=outfile, append=TRUE)
      cat(mat[i,which(!is.na(mat[i,]))],sep=",", file=outfile, append=TRUE)
      cat("}\n",file=outfile, append=TRUE)
    }
    cat("}\n",file=outfile, append=TRUE)
    
  } else {
    write.matrix.nemo(mat, file)
  }
  
}

# modified function with connectivity matrix starting always with focal patch, and rate matrix starting with 1 when the row is empty
# to ensure no empty rows
write.red.dispersal.matrix <- function(connectivity_matrix, rate_matrix, conn_file, rate_file) {
  # Write connectivity matrix
  cat("{", file = conn_file)
  for (i in 1:dim(connectivity_matrix)[1]) {
    cat("{", file = conn_file, append = TRUE)
    # Start with the focal patch number
    cat(i, file = conn_file, append = TRUE)
    if (!all(is.na(connectivity_matrix[i, ]))) {
      cat(",", file = conn_file, append = TRUE)
      # Write the non-NA elements of the row, excluding the focal patch number itself
      non_na_elements <- connectivity_matrix[i, which(!is.na(connectivity_matrix[i, ]))]
      non_na_elements <- non_na_elements[non_na_elements != i]
      cat(non_na_elements, sep = ",", file = conn_file, append = TRUE)
    }
    cat("}\n", file = conn_file, append = TRUE)
  }
  cat("}\n", file = conn_file, append = TRUE)
  
  # Write rate matrix
  cat("{", file = rate_file)
  for (i in 1:dim(rate_matrix)[1]) {
    cat("{", file = rate_file, append = TRUE)
    if (all(is.na(rate_matrix[i, ]))) {
      # If the entire row is NA, write 1
      cat(1, file = rate_file, append = TRUE)
    } else {
      # Write only the non-NA elements of the row
      cat(rate_matrix[i, which(!is.na(rate_matrix[i, ]))], sep = ",", file = rate_file, append = TRUE)
    }
    cat("}\n", file = rate_file, append = TRUE)
  }
  cat("}\n", file = rate_file, append = TRUE)
}


# ALLENWILLER #

# importing grid 
grid <- vect("Allenwiller_grid_r4m_9844patches.shp",crs=crs("+init=epsg:3035"))

# build distance matrix between each patch (using centroid of the patch)
centroids = centroids(grid)
plot(centroids)
patch_coords = crds(centroids)
colnames(patch_coords) = c("X", "Y")
patch_coords = as.data.frame(patch_coords)  
patch_coords$patch_ID = grid$patch.ID # attach patch ID
distance_matrix = as.matrix(distance(centroids))

# define number of patches
num_patch = ncell(grid)

# --------------   SEED dispersal parameters (here we use Waldi dispersal parameters)
b <- c(0.53606, 0.53606) # shape parameter from NMpi2 
d <- c(28.95, 10)  # mean dispersal distance from NMpi2 and reduced dispersal 

## check the curves
x <- 1:100 
results <- data.frame()
for (i in 1:length(d)) {
  d_i <- d[i]
  b_i <- b[i]
  y_values <- exp_power_kernel(d_i, b_i, x)
  temp <- data.frame(x = x, K = y_values, 
                     d = as.factor(d_i), 
                     b = as.factor(b_i))
  results <- rbind(results, temp)
}

ggplot(results, aes(x = x, y = K, group =d, color = d))+
  geom_line(size = 1) +
  labs(x = "Euclidean Distance",y = "Probability of dispersal")

# final parameters
# define threshold for rate going to 0 at tot meters from focal tree
d_t <- 50
d <- 28.95
b <- 0.53606
threshold = exp_power_kernel(d,b, d_t) ## otherwise, just keep d
threshold

# building matrices and plot
dispersal_matrices = build.reduced.dispersal.matrices(num_patch, distance_matrix, d, b, threshold)

dispersal_matrix <- dispersal_matrices$dispersal_matrix
connectivity_matrix <- dispersal_matrices$connectivity_matrix
rate_matrix <- dispersal_matrices$rate_matrix

dim(dispersal_matrix)
dim(connectivity_matrix)
dim(rate_matrix)

# check that rate for one (random) row sums to 1
sum(rate_matrix[1500, ], na.rm = T)

# check dispresal of one random focal patch
grid_r <- terra::rast(ext = ext(grid), resolution = 4, crs = "EPSG:3035")
# grid dimensions 
n_rows = dim(grid_r)[1]
n_cols = dim(grid_r)[2]

library(plot.matrix)
par(mar = c(4, 4, 4, 10))
plot(matrix(dispersal_matrix[,5000],ncol=n_cols,byrow= T))

# check distance against dispersal
distance_v <- distance_matrix[lower.tri(distance_matrix, diag = T)][order(row(distance_matrix)[lower.tri(row(distance_matrix), diag = T)])]
dispersal_v <- dispersal_matrix[lower.tri(dispersal_matrix, diag = T)][order(row(dispersal_matrix)[lower.tri(row(dispersal_matrix), diag = T)])]
plot(distance_v, dispersal_v)

## visualize the NEMO matrices
        par(mfrow = c(1,3))
        disp_raster = rast(dispersal_matrix)
        plot(disp_raster$lyr.1, main = "Dispersal matrix")
        conn_raster = rast(connectivity_matrix)
        plot(conn_raster$lyr.1, main = "Connectivity matrix")
        rate_raster = rast(rate_matrix)
        plot(rate_raster$lyr.1, main = "Reduced matrix")
        
# (OPTIONAL) Write full dispersal matrix to file (if needed)
# write.dispersal.matrix(dispersal_matrix, paste0(input_path, "NEMO_seed_dispersal_matrix_a.txt"), discard.NA = TRUE)

# write connectivity and rate matrix to files 
# save d = mean dispersal distance of the kernel
write.red.dispersal.matrix(connectivity_matrix, rate_matrix, 
                           paste0(input_path, "NEMO_seed_connectivity_matrix_a_d", d, ".txt"), 
                           paste0(input_path, "NEMO_seed_rate_matrix_a_d", d, ".txt"))





# ---------------- POLLEN dispersal parameters (here we use Waldi dispersal parameters)
d <- c(109.41, 50, 20) # mean pollen dispersal distance from NMpi2 and reduced version
b <- c(0.85, 0.85,0.85) # shape parameter from NMpi2 ML seed dispersal 

## define the d_threshold
x <- 1:200 
results <- data.frame()
for (i in 1:length(d)) {
  d_i <- d[i]
  b_i <- b[i]
  y_values <- exp_power_kernel(d_i, b_i, x)
  temp <- data.frame(x = x, K = y_values, 
                     d = as.factor(d_i), 
                     b = as.factor(b_i))
  results <- rbind(results, temp)
}

ggplot(results, aes(x = x, y = K, group =d, color = d))+
  geom_line(size = 1) +
  labs(x = "Euclidean Distance",y = "Probability of dispersal")

# final parameters
# define threshold for rate going to 0 at tot meters from focal tree
#d = 109.41
#d_t = 200
#b = 0.85

#d = 50
#b = 0.85
#d_t = 100

#d = 20
#b = 0.85
#d_t = 50
threshold = exp_power_kernel(d,b, d_t) ## otherwise, just keep d
threshold

# building matrices and plot
dispersal_matrices = build.reduced.dispersal.matrices(num_patch, distance_matrix, d, b, threshold)

dispersal_matrix <- dispersal_matrices$dispersal_matrix
connectivity_matrix <- dispersal_matrices$connectivity_matrix
rate_matrix <- dispersal_matrices$rate_matrix

# check dimansions
dim(dispersal_matrix)
dim(connectivity_matrix)
dim(rate_matrix)

# check that rate for one (random) row sums to 1
sum(rate_matrix[1500, ], na.rm = T)

# check dispresal of one random focal patch
library(plot.matrix)
plot(matrix(dispersal_matrix[,5000],ncol=n_cols,byrow= T))

# check distance against dispersal
distance_v <- distance_matrix[lower.tri(distance_matrix, diag = T)][order(row(distance_matrix)[lower.tri(row(distance_matrix), diag = T)])]
dispersal_v <- dispersal_matrix[lower.tri(dispersal_matrix, diag = T)][order(row(dispersal_matrix)[lower.tri(row(dispersal_matrix), diag = T)])]
plot(distance_v, dispersal_v)

        ## visualize the NEMO matrices
        par(mfrow = c(1,3))
        disp_raster = rast(dispersal_matrix)
        plot(disp_raster$lyr.1, main = "Dispersal matrix")
        conn_raster = rast(connectivity_matrix)
        plot(conn_raster$lyr.1, main = "Connectivity matrix")
        rate_raster = rast(rate_matrix)
        plot(rate_raster$lyr.1, main = "Reduced matrix")
        
# (OPTIONAL) Write full dispersal matrix to file (if needed)
# write.dispersal.matrix(dispersal_matrix, paste0(input_path, "NEMO_pollen_dispersal_matrix_a.txt"), discard.NA = TRUE)

# write connectivity and rate matrix to files 
write.red.dispersal.matrix(connectivity_matrix, rate_matrix, 
                           paste0(input_path, "NEMO_pollen_connectivity_matrix_a_d", d, ".txt"), 
                           paste0(input_path, "NEMO_pollen_rate_matrix_a_d", d, ".txt"))

# WALDI #

# importing grid 
grid <- vect("Waldi_grid_4m_5396patches.shp",crs=crs("+init=epsg:3035"))

# build distance matrix between each patch (centroid of the patch)
# get centroid of the patches and coordinates
centroids = centroids(grid)
plot(centroids)
patch_coords = crds(centroids)
colnames(patch_coords) = c("X", "Y")
patch_coords = as.data.frame(patch_coords)  # coordinates as a data frame
patch_coords$patch_ID = grid$patch.ID # attaching patch ID
distance_matrix = as.matrix(distance(centroids))

# define number of patches
num_patch = ncell(grid)

# ------------------ SEED dispersal parameters 
b <- c(0.53606, 0.53606) # shape parameter from NMpi2 
d <- c(28.95, 10)  # mean dispersal distance from NMpi2 and reduced dispersal 

## check the curves
x <- 1:100 
results <- data.frame()
for (i in 1:length(d)) {
  d_i <- d[i]
  b_i <- b[i]
  y_values <- exp_power_kernel(d_i, b_i, x)
  temp <- data.frame(x = x, K = y_values, 
                     d = as.factor(d_i), 
                     b = as.factor(b_i))
  results <- rbind(results, temp)
}

ggplot(results, aes(x = x, y = K, group =d, color = d))+
  geom_line(size = 1) +
  labs(x = "Euclidean Distance",y = "Probability of dispersal")

# final parameters
# define threshold for rate going to 0 at tot meters from focal tree
d_t <- 50
d <- 10
b <- 0.53606
threshold = exp_power_kernel(d,b, d_t) ## otherwise, just keep d
threshold

# building matrices and plot
dispersal_matrices = build.reduced.dispersal.matrices(num_patch, distance_matrix, d, b, threshold)

dispersal_matrix <- dispersal_matrices$dispersal_matrix
connectivity_matrix <- dispersal_matrices$connectivity_matrix
rate_matrix <- dispersal_matrices$rate_matrix

dim(dispersal_matrix)
dim(connectivity_matrix)
dim(rate_matrix)

# check that rate for one (random) row sums to 1
sum(rate_matrix[1500, ], na.rm = T)

# check dispresal of one random focal patch
grid_r <- terra::rast(ext = ext(grid), resolution = 4, crs = "EPSG:3035")
# grid dimensions 
n_rows = dim(grid_r)[1]
n_cols = dim(grid_r)[2]

library(plot.matrix)
plot(matrix(dispersal_matrix[,1500],ncol=n_cols,byrow= T))

# check distance against dispersal
distance_v <- distance_matrix[lower.tri(distance_matrix, diag = T)][order(row(distance_matrix)[lower.tri(row(distance_matrix), diag = T)])]
dispersal_v <- dispersal_matrix[lower.tri(dispersal_matrix, diag = T)][order(row(dispersal_matrix)[lower.tri(row(dispersal_matrix), diag = T)])]
plot(distance_v, dispersal_v)

## visualize the NEMO matrices
      par(mfrow = c(1,3))
      disp_raster = rast(dispersal_matrix)
      plot(disp_raster$lyr.1, main = "Dispersal matrix")
      conn_raster = rast(connectivity_matrix)
      plot(conn_raster$lyr.1, main = "Connectivity matrix")
      rate_raster = rast(rate_matrix)
      plot(rate_raster$lyr.1, main = "Reduced matrix")

# (OPTIONAL) Write full dispersal matrix to file (if needed)
# write.dispersal.matrix(dispersal_matrix, paste0(input_path, "NEMO_seed_dispersal_matrix_w.txt"), discard.NA = TRUE)

# write connectivity and rate matrix to files 
write.red.dispersal.matrix(connectivity_matrix, rate_matrix, 
                           paste0(input_path, "NEMO_seed_connectivity_matrix_w_d", d, ".txt"), 
                           paste0(input_path, "NEMO_seed_rate_matrix_w_d", d, ".txt"))





# ---------------POLLEN dispersal parameters (here we use Waldi dispersal parameters)
d <- c(109.41, 50, 20) # mean pollen dispersal distance from NMpi2 and reduced version
b <- c(0.85, 0.85, 0.85) # shape parameter from NMpi2 ML seed dispersal 

## define the d_threshold
x <- 1:200 
results <- data.frame()
for (i in 1:length(d)) {
  d_i <- d[i]
  b_i <- b[i]
  y_values <- exp_power_kernel(d_i, b_i, x)
  temp <- data.frame(x = x, K = y_values, 
                     d = as.factor(d_i), 
                     b = as.factor(b_i))
  results <- rbind(results, temp)
}

ggplot(results, aes(x = x, y = K, group =d, color = d))+
  geom_line(size = 1) +
  labs(x = "Euclidean Distance",y = "Probability of dispersal")

# final parameters
# define threshold for rate going to 0 at tot meters from focal tree
#d = 109.41
#d_t = 200
#b = 0.85

#d = 50
#b = 0.85
#d_t = 100

#d = 20
#b = 0.85
#d_t = 50
threshold = exp_power_kernel(d,b, d_t) ## otherwise, just keep d
threshold

# building matrices and plot
dispersal_matrices = build.reduced.dispersal.matrices(num_patch, distance_matrix, d, b, threshold)

dispersal_matrix <- dispersal_matrices$dispersal_matrix
connectivity_matrix <- dispersal_matrices$connectivity_matrix
rate_matrix <- dispersal_matrices$rate_matrix

dim(dispersal_matrix)
dim(connectivity_matrix)
dim(rate_matrix)

# check dispresal of one random focal patch
library(plot.matrix)
plot(matrix(dispersal_matrix[800,],ncol=n_cols,byrow= T))

# check distance against dispersal
distance_v <- distance_matrix[lower.tri(distance_matrix, diag = T)][order(row(distance_matrix)[lower.tri(row(distance_matrix), diag = T)])]
dispersal_v <- dispersal_matrix[lower.tri(dispersal_matrix, diag = T)][order(row(dispersal_matrix)[lower.tri(row(dispersal_matrix), diag = T)])]
plot(distance_v, dispersal_v)

      ## convert to raster and plot
      par(mfrow = c(1,3))
      disp_raster = rast(dispersal_matrix)
      plot(disp_raster$lyr.1, main = "Dispersal matrix")
      conn_raster = rast(connectivity_matrix)
      plot(conn_raster$lyr.1, main = "Connectivity matrix")
      rate_raster = rast(rate_matrix)
      plot(rate_raster$lyr.1, main = "Reduced matrix")

# (OPTIONAL) Write full dispersal matrix to file (if needed)
# write.dispersal.matrix(dispersal_matrix, paste0(input_path, "NEMO_pollen_dispersal_matrix_w.txt"), discard.NA = TRUE)

# write connectivity and rate matrix to files 
write.red.dispersal.matrix(connectivity_matrix, rate_matrix, 
                           paste0(input_path, "NEMO_pollen_connectivity_matrix_w_d", d, ".txt"), 
                           paste0(input_path, "NEMO_pollen_rate_matrix_w_d", d, ".txt"))



# -------------- FLAT dispersal

### create a flat dispresal matrix with constant value equal to the average dispersal of the exponential kernel

## remove the distance threshold (they will spread evenly throughout the space)
threshold = 0.0000005 
# mean dispersal probability  = 1 / number of patches
### give the num_patch of Waldi or Allewniller
mean_K = 1/num_patch

# build dispersal matrices using the flat kernel
build.reduced.dispersal.matrices.flat <- function(num_patch, distance_matrix, d, b, d_thresh = threshold) {
  stopifnot(dim(distance_matrix)[1] == num_patch)
  
  dispersal_matrix <- matrix(0, nrow = num_patch, ncol = num_patch)
  rate_matrix <- matrix(0, nrow = num_patch, ncol = num_patch)
  connectivity_matrix <- matrix(NA, nrow = num_patch, ncol = num_patch)
  
  
  for (i in 1:num_patch) {
    for (j in 1:num_patch) {
      # just constant value
      rate_matrix[i, j] <- mean_K
    }
  }
  
  # Define full dispersal matrix
  dispersal_matrix <- rate_matrix
  
  for (i in 1:num_patch) {
    # Get indices of the patches sorted by dispersal probability in descending order
    ord <- order(rate_matrix[i, ], decreasing = TRUE)
    
    # Re-order the values in row i of reduced rate matrix, connectivity matrix
    rate_matrix[i, ] <- rate_matrix[i, ord]
    connectivity_matrix[i, ] <- ord
    
    # Remove elements below the threshold by setting them to NA 
    to_remove <- which(rate_matrix[i, ] < d_thresh)
    rate_matrix[i, to_remove] <- NA
    connectivity_matrix[i, to_remove] <- NA
    
    # Normalize dispersal rates to sum to 1
    sum_d <- sum(rate_matrix[i, ], na.rm = TRUE)
    rate_matrix[i, !is.na(rate_matrix[i, ])] <- rate_matrix[i, !is.na(rate_matrix[i, ])] / sum_d
  }
  
  return(list(dispersal_matrix = dispersal_matrix, connectivity_matrix = connectivity_matrix, rate_matrix = rate_matrix))
}


dispersal_matrices_flat <- build.reduced.dispersal.matrices.flat(num_patch, distance_matrix, d, b, threshold)

dispersal_matrix <- dispersal_matrices_flat$dispersal_matrix
connectivity_matrix <- dispersal_matrices_flat$connectivity_matrix
rate_matrix <- dispersal_matrices_flat$rate_matrix

dim(dispersal_matrix)
dim(connectivity_matrix)
dim(rate_matrix)

# check that rate for one (random) row sums to 1
sum(rate_matrix[1500, ], na.rm = T)

# check dispresal of one random focal patch
grid_r <- terra::rast(ext = ext(grid), resolution = 4, crs = "EPSG:3035")
# grid dimensions 
n_rows = dim(grid_r)[1]
n_cols = dim(grid_r)[2]

library(plot.matrix)
plot(matrix(dispersal_matrix[,1500],ncol=n_cols,byrow= T))

      # check distance against dispersal
      distance_v <- distance_matrix[lower.tri(distance_matrix, diag = T)][order(row(distance_matrix)[lower.tri(row(distance_matrix), diag = T)])]
      dispersal_v <- dispersal_matrix[lower.tri(dispersal_matrix, diag = T)][order(row(dispersal_matrix)[lower.tri(row(dispersal_matrix), diag = T)])]
      plot(distance_v, dispersal_v)
      
      ## visualize the NEMO matrices
      par(mfrow = c(1,3))
      disp_raster = rast(dispersal_matrix)
      plot(disp_raster$lyr.1, main = "Dispersal matrix")
      conn_raster = rast(connectivity_matrix)
      plot(conn_raster$lyr.1, main = "Connectivity matrix")
      rate_raster = rast(rate_matrix)
      plot(rate_raster$lyr.1, main = "Reduced matrix")

# (OPTIONAL) Write full dispersal matrix to file (if needed)
# write.dispersal.matrix(dispersal_matrix, paste0(input_path, "NEMO_seed_dispersal_matrix_w.txt"), discard.NA = TRUE)


# write connectivity and rate matrix to files (both seed and pollen are the same, just change the name)
# WALDI #
write.red.dispersal.matrix(connectivity_matrix, rate_matrix,paste0(input_path, "NEMO_seed_connectivity_matrix_w_flat.txt"), paste0(input_path, "NEMO_seed_rate_matrix_w_flat.txt"))
file.copy(from =  paste0(input_path, "NEMO_seed_connectivity_matrix_w_flat.txt"), to = paste0(input_path, "NEMO_pollen_connectivity_matrix_w_flat.txt")) 
file.copy(from =paste0(input_path, "NEMO_seed_rate_matrix_w_flat.txt"),  to =paste0(input_path, "NEMO_pollen_rate_matrix_w_flat.txt"))

# ALLENWILLER #
write.red.dispersal.matrix(connectivity_matrix, rate_matrix,paste0(input_path, "NEMO_seed_connectivity_matrix_a_flat.txt"), paste0(input_path, "NEMO_seed_rate_matrix_a_flat.txt"))
file.copy(from =  paste0(input_path, "NEMO_seed_connectivity_matrix_a_flat.txt"), to = paste0(input_path, "NEMO_pollen_connectivity_matrix_a_flat.txt")) 
file.copy(from =paste0(input_path, "NEMO_seed_rate_matrix_a_flat.txt"),  to =paste0(input_path, "NEMO_pollen_rate_matrix_a_flat.txt"))

