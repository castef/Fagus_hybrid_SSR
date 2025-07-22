#### NEMO simulation - analysis of results of demographic test ######

library(ggplot2)
library(dplyr)
library(tidyr)
library(raster)
library(terra)
library(data.table)
library(viridis)
library(readr)
library(stringr)
library(ggh4x)
library(sf)
library(xlsx)

rm(list=ls())

# path where to move results for storage and analysis
store_path = "~/NEMO/NEMO_output"

# path to get other files
setwd(dir="~/Data_clean")

# WALDI #

# load simulation grid
grid <- vect("Waldi_grid_4m_5396patches.shp", crs = "EPSG:3035")
# load sampling scheme to extract the results 
sampling <- vect("Waldi_sampling_scheme.shp", crs = "EPSG:3035")
sampling = sampling[!is.na(geom(sampling)[,"x"]),]
# find the patches of the grid corresponding to the sampling
# patch where centroid intersects the sampling scheme
centroid_match <- terra::intersect(terra::centroids(grid), sampling)

# intersection between grid and sampling
area_match <- terra::intersect(grid, sampling)

# area of intersection for each cell
area_match$intersection_area <- terra::expanse(area_match)
# set threshold of intersection area to retain the cell (30%)
cell_area <- terra::expanse(grid)[1] # area of each cell (4 x 4 m)
area_thres_match <- area_match[area_match$intersection_area > (cell_area * 0.20), ]

# get patch ID from combining both approaches
all_match_patch <- unique(c(centroid_match$patch.ID, area_thres_match$patch.ID))

# check
plot(grid)
plot(sampling, add=TRUE)
matched_grid <- grid[grid$patch.ID %in% all_match_patch, ]
plot(matched_grid, add=TRUE, col='blue', alpha=0.5)

# compare the total sampling area in reality and the total sampling area in the simulation
sampling_area <- sum(terra::expanse(sampling))
matched_grid <- grid[grid$patch.ID %in% all_match_patch, ]
matched_area <- sum(terra::expanse(matched_grid))
sampling_area
matched_area

# load observed data of circle plots
obs_cp <- read.csv("Waldi_circleplot_info_withcoord.csv")
# load observed adult data
obs_ad <- read.csv("Waldi_adults_Coordinates_metric.csv")

# check (size of seedlings circles proportional to the observed number of seedlings and juvenils)
plot(grid)
plot(sampling, add = T)
plot(vect(obs_cp, geom =  c("x", "y"), crs = crs(grid)), cex = obs_cp$NClass1Class3*0.05, add = T)
plot(vect(obs_cp, geom =  c("x", "y"), crs = crs(grid)), cex = obs_cp$NClass4Class5*0.05, add = T, col = "red")
plot(vect(obs_ad, geom =  c("x", "y"), crs = crs(grid)), add = T, col = "blue")

# count number of seedlings and juveniles in the sampled area (in the circle plots)
obs_seedl_tot_wal = sum(obs_cp$NClass1Class3)
obs_juv_tot_wal = sum(obs_cp$NClass4Class5)
# count number of adults
obs_ad_tot_wal = nrow(obs_ad)

# scale the results on the sampling area
obs_seedl_scaled_wal = obs_seedl_tot_wal / sampling_area
obs_juv_scaled_wal = obs_juv_tot_wal / sampling_area
obs_ad_scaled_wal = obs_ad_tot_wal / sampling_area

# max number of seedlings and juveniles counted in the cp (for carrying capacity)
max(obs_cp$NClass1Class3)
max(obs_cp$NClass4Class5)

# store specific Waldi grid info for later
grid_wal <- grid
all_match_patch_wal <- all_match_patch
sampling_area_sim_wal <- matched_area


# ALLENWILLER #

# load simulation grid
grid <- vect("Allenwiller_grid_r4m_9844patches.shp")
# load sampling scheme to extract the results 
sampling <- vect("Allenwiller_sampling_scheme.shp")
sampling = sampling[!is.na(geom(sampling)[,"x"]),]

# find the patches of the grid corresponding to the sampling
# patch where centroid intersects the sampling scheme
centroid_match <- terra::intersect(terra::centroids(grid), sampling)
# intersection between grid and sampling
area_match <- terra::intersect(grid, sampling)
# area of intersection for each cell
area_match$intersection_area <- terra::expanse(area_match)
# set threshold of intersection area to retain the cell (30%)
cell_area <- terra::expanse(grid)[1] # area of each cell (4 x 4 m)
area_thres_match <- area_match[area_match$intersection_area > (cell_area * 0.20), ]
# get cell ID from combining both approaches
all_match_patch <- unique(c(centroid_match$patch.ID, area_thres_match$patch.ID))

# check
plot(grid)
plot(sampling, add=TRUE)
matched_grid <- grid[grid$patch.ID %in% all_match_patch, ]
plot(matched_grid, add=TRUE, col='blue', alpha=0.5)

# compare the total sampling area in reality and the total sampling area in the simulation
sampling_area <- sum(terra::expanse(sampling))
matched_grid <- grid[grid$patch.ID %in% all_match_patch, ]
matched_area <- sum(terra::expanse(matched_grid))
sampling_area
matched_area

# load observed data of circle plots
obs_cp <- read.csv("Allenwiller_circleplot_info_withcoord.csv")
# load observed adult data
obs_ad <- read.csv("Allenwiller_adults_Coordinates_metric.csv")

# check (size of seedlings circles proportional to the observed number of seedlings and juvenils)
plot(grid)
plot(sampling, add = T)
plot(vect(obs_cp, geom =  c("x", "y"), crs = crs(grid)), cex = obs_cp$CountedC1.new.*0.01, add = T)
plot(vect(obs_cp, geom =  c("x", "y"), crs = crs(grid)), cex = obs_cp$CountedNrC2.new.*0.01, add = T, col = "red")
plot(vect(obs_ad, geom =  c("x", "y"), crs = crs(grid)), add = T, col = "blue")

# count number of seedlings and juveniles in the sampled area (in the circle plots)
obs_seedl_tot_al = sum(obs_cp$CountedC1.new.)
obs_juv_tot_al = sum(obs_cp$CountedNrC2.new.)
# count number of adults
obs_ad_tot_al = nrow(obs_ad)

# scale the results on the sampling area
obs_seedl_scaled_al = obs_seedl_tot_al / sampling_area
obs_juv_scaled_al = obs_juv_tot_al / sampling_area
obs_ad_scaled_al = obs_ad_tot_al / sampling_area

# max number of seedlings and juveniles counted in the cp (for carrying capacity)
max(obs_cp$CountedC1.new.)
max(obs_cp$CountedNrC2.new.)

# store specific Allenwiller grid info for later
grid_al <- grid
all_match_patch_al <- all_match_patch
sampling_area_sim_al <- matched_area

# loading simulated data of demographic test
setwd(store_path)

# select the type of simulation
scenario <- "spatial" 
# patterns for reading files in demo folder
folder_pattern <- "_demo_test"
# folder_pattern <- "Allenwiller_demo_test"
base_pattern <- "_demo_spatial_k\\d+_a\\d+_s1\\d+_b\\d+\\.\\d+_results\\.txt$"
file_pattern <- paste0("^[A-Za-z]+", base_pattern)

# all files matching pattern
file_list <- list.files(path = store_path,pattern = file_pattern,recursive = TRUE, full.names = TRUE)
# keep only files in demo test folder
file_list <- file_list[grepl(folder_pattern, dirname(file_list))]

# to get only the patch for the specific grid 
grid_al_r <- terra::rast(ext = ext(grid_al), resolution = 4, crs = "EPSG:3035")
grid_wal_r <- terra::rast(ext = ext(grid_wal), resolution = 4, crs = "EPSG:3035")

grid_info <- list(
  "Allenwiller" = list(n_rows = dim(grid_al_r)[2], n_cols = dim(grid_al_r)[1], patch_ids = all_match_patch_al),
  "Waldi" = list(n_rows =  dim(grid_wal_r)[2], n_cols = dim(grid_wal_r)[1], patch_ids = all_match_patch_wal))


#  empty list to store processed data
data_list2 <- list()

# read and process files
for (i in seq_along(file_list)) {
  file <- file_list[i]
  
  tryCatch({
    data <- fread(file)
    
    # stand
    stand <- ifelse(grepl("Allenwiller", basename(file)), "Allenwiller", "Waldi")
    
    # patch.ID for the simulation sapmling (for each stand)
    sampling_patch <- grid_info[[stand]]$patch_ids
    
    # NB: old nemo version stage 0 is coded as "off.fem.p", while new versions are "a0.fem.p" 
    col_names <- names(data)
    col_off <- any(grepl("^off\\.fem\\.p", col_names))  # for old nemo versions
    col_a0 <- any(grepl("^a0\\.fem\\.p", col_names)) # for new nemo version
    stages <- if (col_off) {
      c("off", "a1", "a2", "a3")
    } else if (col_a0) {
      c("a0", "a1", "a2", "a3")
    } else {
      # otherwise use the original stages
      c("off", "a1", "a2", "a3")
    }
    # create col names
    generate_col_names <- function(stage, sex) {
      paste0(stage, ".", sex, ".p", sampling_patch)
    }
    
    # combine column names for female and male columns, also count all the adults
    cols_to_read <- c("replicate", "generation", "a3.tot",
                      unlist(lapply(stages, generate_col_names, sex = "fem")),
                      unlist(lapply(stages, generate_col_names, sex = "mal")))
    # read only the columns of interest (only specific patches to sample)
    cols_exist <- cols_to_read[cols_to_read %in% names(data)]
    data <- data[, ..cols_exist]
    
    # get metadata from file and folder name
    file_name <- gsub(".txt", "", basename(file))
    k <- as.numeric(sub(".*_k(\\d+)_.*", "\\1", basename(file)))
    a <- as.numeric(sub(".*_a(\\d+)_.*", "\\1", basename(file)))
    s1 <-  as.numeric(sub(".*_s1(\\d+)_.*", "\\1", basename(file)))/10
    b <- as.numeric(sub(".*_b([0-9]+\\.?[0-9]*e?-?[0-9]*).*", "\\1", basename(file))) 
    b <- as.numeric(format(b, scientific = FALSE))
    
    # add infos as columns
    data$file_name <- file_name
    data$stand <- stand
    data$k <- k
    data$a <- a
    data$s1 <- s1
    data$b <- as.numeric(b)
    
    # check that the values are stored properly from the filename
    cat(sprintf("Processing file %d/%d: %s (k=%s, a=%s, s1=%s,b=%s)\n",i, length(file_list), file, k, a, s1, b))

    # reshape data for each population stage
    stages_data_list <- lapply(stages, function(stage) {
      # Melt females
      fem_data <- melt(data, 
                       id.vars = c("stand", "replicate", "generation","a3.tot", "k", "a", "s1", "b"), 
                       measure.vars = patterns(paste0("^", stage, ".fem.p")), 
                       variable.name = "patch", 
                       value.name = "N_ind")
      
      # get cell ID from column name
      fem_data[, patch := as.numeric(gsub(paste0(stage, ".fem.p"), "", patch, fixed = TRUE))]
      
      # filter only patch corresponding to sampling patch
      fem_data <- fem_data[patch %in% sampling_patch]
      
      # add stage column ad grid infos
      fem_data$stage <- stage
      fem_data$n_rows <- grid_info[[stand]]$n_rows
      fem_data$n_cols <- grid_info[[stand]]$n_cols
      
      return(fem_data)
    })
    
    # combine the stage data for 1 file
    file_data <- rbindlist(stages_data_list)
    
    # store the processed file data with a unique identifier (to avoid overwriting if there are files with same name by mistake)
    data_list2[[paste0(file_name, "_", i)]] <- file_data
    gc()
  }, error = function(e) {
    # chek if there are errors
    cat("Error processing file:", file, "\n")
    print(e)
  })
}

# combine all files into 1 df
combined_fem_data <- rbindlist(data_list2)
rm(data, data_list2, file_data, stages_data_list)
gc()

# convert "off" stage into "a0" (new version)
combined_fem_data$stage[combined_fem_data$stage == "off"] <- "a0"

# remove b = 0.002 
combined_fem_data <- combined_fem_data[combined_fem_data$b != 0.002 , ]

# save results simulations combined 
setwd(dir="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Data_clean")
saveRDS(combined_fem_data, "NEMO/Demographic_test_results.RDS")

# tot individuals across the sampled area (summarize over patches)
combined_fem_data <- readRDS( "NEMO/Demographic_test_results.RDS")
summary_data <- combined_fem_data %>%
  group_by(replicate, generation, stand, stage, a3.tot, k, a, s1, b) %>%
  summarize(N_stage = sum(N_ind, na.rm = TRUE)) %>%
  ungroup()

# aggregate stage 0 and stage 1, as they are together in the observed data
# consider both off and a0 version
nemo_stage <- c("a0" = "Stage 0", "a1" = "Stage 1", "a2" = "Stage 2", "a3" = "Stage 3")
agg_stage <- c("a0" = "Stage 1", "a1" = "Stage 1", "a2" = "Stage 2", "a3" = "Stage 3")
summary_data$nemo_stage <- nemo_stage[summary_data$stage]
summary_data$agg_stage <- agg_stage[summary_data$stage]

# summary on aggreagate stages, but keep adults as they are
summary_data_agg = summary_data %>% 
  group_by(replicate,generation,stand,a3.tot,k, a, s1, b,agg_stage) %>%
  summarise(N_agg_stage=sum(N_stage))
summary_data_agg$N_agg_stage[summary_data_agg$agg_stage=="Stage 3"]=summary_data_agg$a3.tot[summary_data_agg$agg_stage=="Stage 3"]

sampling_area_sim_wal <- avg_sampling_area_per_stand %>% filter(stand == "Waldi") %>%  pull(avg_sampling_area)
sampling_area_sim_al <- avg_sampling_area_per_stand %>% filter(stand == "Allenwiller") %>%pull(avg_sampling_area)

# normalize numbers on the sampling area
summary_data_agg$N_agg_stage_scaled <- ifelse(summary_data_agg$stand == "Waldi", summary_data_agg$N_agg_stage/sampling_area_sim_wal, summary_data_agg$N_agg_stage/sampling_area_sim_al)
summary_data_agg$N_agg_stage_scaled <- ifelse(summary_data_agg$stand == "Waldi", summary_data_agg$N_agg_stage/sampling_area_sim_wal, summary_data_agg$N_agg_stage/sampling_area_sim_al)

#observed values dataframe
obs_val=data.frame(agg_stage=c("Stage 1","Stage 2","Stage 3"),
                   stand = c("Allenwiller", "Allenwiller", "Allenwiller", "Waldi","Waldi","Waldi" ),
                   N_obs = c(obs_seedl_tot_al, obs_juv_tot_al ,obs_ad_tot_al,obs_seedl_tot_wal, obs_juv_tot_wal,obs_ad_tot_wal),
                   N_obs_scaled=c(obs_seedl_scaled_al,obs_juv_scaled_al,obs_ad_scaled_al,obs_seedl_scaled_wal,obs_juv_scaled_wal,obs_ad_scaled_wal))


# chekc simulated data for each stand separately
# add source data type for legend
obs_val$data_source <- "Observed data"
summary_data_agg$data_source <- "Simulated data"

#3 check number of replicates: 
replicates_check <- summary_data_agg %>%
  group_by(stand,generation,agg_stage, a, s1, b, k) %>%
  summarise(num_replicates = n(), .groups = "drop") %>%
  arrange(desc(num_replicates))
replicates_check[replicates_check$num_replicates > 10, ]
replicates_check[replicates_check$num_replicates < 10, ]

pop_trend_w <- ggplot() +
  geom_line(data = subset(summary_data_agg, stand == "Waldi"), aes(x = generation, y = N_agg_stage_scaled,  group = interaction(agg_stage, replicate),color = agg_stage,linetype = data_source)) +
  geom_hline(data = subset(obs_val, stand == "Waldi"),aes(yintercept = N_obs_scaled,color = agg_stage,linetype = data_source)) +
  scale_linetype_manual(name = "Data type",
                        values = c("Observed data" = "dashed",
                                   "Simulated data" = "solid")) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", 
        strip.text = element_text(size = 10),
        legend.key = element_rect(colour = "white"),
        legend.key.spacing.y = unit(0.4, 'cm'), 
        strip.background = element_rect(colour = "black", fill = "white"), 
        plot.title = element_text(hjust = 0.5))+
  ggtitle("Wäldi")+
  
  facet_grid(b + a ~ k + s1,labeller = labeller(
    k = function(x) paste0("k = ", x),
    s1 = function(x) paste0("s1 = ", x),
    a = function(x) paste0("a = ", x),
    b = function(x) paste0("b = ", x)
  )) +
  labs(x = "Years",y = "Number of Individuals scaled for sampled area",color = "Nemo Stage") +
  scale_x_continuous(limits = c(0,150))+
  scale_y_continuous(labels = function(x) format(x, big.mark = "", scientific = FALSE)) +
  scale_y_sqrt() # square root transformation to better visualization
pop_trend_w

pop_trend_a <- ggplot() +
  geom_line(data = subset(summary_data_agg, stand == "Allenwiller"), aes(x = generation, y = N_agg_stage_scaled,  
                                                                         group = interaction(agg_stage, replicate),color = agg_stage,linetype = data_source)) +
  geom_hline(data = subset(obs_val, stand == "Allenwiller"),aes(yintercept = N_obs_scaled,color = agg_stage,linetype = data_source)) +
  scale_linetype_manual(name = "Data type",
                        values = c("Observed data" = "dashed",
                                   "Simulated data" = "solid")) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", 
        strip.text = element_text(size = 10),
        legend.key = element_rect(colour = "white"),
        legend.key.spacing.y = unit(0.4, 'cm'), 
        strip.background = element_rect(colour = "black", fill = "white"), 
        plot.title = element_text(hjust = 0.5))+
  ggtitle("Allenwiller")+
   facet_grid(b + a ~ k + s1,labeller = labeller(
    k = function(x) paste0("k = ", x),
    s1 = function(x) paste0("s1 = ", x),
    a = function(x) paste0("a = ", x),
    b = function(x) paste0("b = ", x))) +
  labs(x = "Years",y = "Number of Individuals scaled for sampled area",color = "Nemo Stage") +
  scale_x_continuous(limits = c(0,150))+
  scale_y_continuous(labels = function(x) format(x, big.mark = "", scientific = FALSE)) +
  scale_y_sqrt() # square root transformation to better visualization
pop_trend_a

FigS3 <- ggpubr::ggarrange(
  pop_trend_a, pop_trend_w, 
  nrow = 2, ncol = 1,
  legend = "right",
  common.legend = TRUE, 
  align = "v"
)  

FigS3

path="~/Hybridization/"

png(filename=paste0(path, "Manuscript/Figures/FigureS2_panel.png"), width =12000,height=10000, res=900)
plot(FigS3)
dev.off()


## get the best parameters for each stand
# merge simulated with observed
dim(summary_data_agg)
summary_data_agg=merge(summary_data_agg,obs_val,by=c("agg_stage", "stand"))
dim(summary_data_agg)

# get only year 70 (corresponding to 100 from planting adults approximately)
end_sim_obs=subset(summary_data_agg,generation==70)

# create combinations of params
tab_for_loop = unique(end_sim_obs[,c("b","a","k","s1", "stand")])

# chi square test to check the best model
end_sim_obs2=NULL
# for each stand
for(stand_name in unique(tab_for_loop$stand)) {
  tab_stand <- subset(tab_for_loop, stand == stand_name)
  end_sim_stand <- subset(end_sim_obs, stand == stand_name)

  for(i in 1:nrow(tab_stand)){
    bidon=subset(end_sim_stand,b==tab_stand$b[i] & a==tab_stand$a[i] & k==tab_stand$k[i] & s1==tab_stand$s1[i])
    tab_stand$chi[i]=sum(((bidon$N_agg_stage-bidon$N_obs)^2)/bidon$N_obs)
    tab_stand$chi_scaled[i]=sum(((bidon$N_agg_stage_scaled-bidon$N_obs_scaled)^2)/bidon$N_obs_scaled)
    end_sim_obs2=rbind(end_sim_obs2,bidon)
  }
  tab_for_loop[tab_for_loop$stand == stand_name, c("chi", "chi_scaled")] <- tab_stand[, c("chi", "chi_scaled")]
}


# WALDI #

tab_for_loop_wal=subset(tab_for_loop, stand == "Waldi")
tab_for_loop_wal = tab_for_loop_wal[order(tab_for_loop_wal$chi_scaled,decreasing=FALSE),]
head(tab_for_loop_wal)

# check hte real numbers of the best model
check <- subset(end_sim_obs2, k == tab_for_loop_wal$k[1] & s1 == tab_for_loop_wal$s1[1] & b == tab_for_loop_wal$b[1] &  a == tab_for_loop_wal$a[1])

wal_best_param <- ggplot(subset(check, stand == "Waldi"), aes(x = agg_stage, y = N_agg_stage, fill = agg_stage)) +
  geom_boxplot() +
  theme_bw() +
  labs(
    title = "",
    x = "",
    y = "Number of Individuals",
    fill = "Simulated data",
    color = "Observed data"
  ) +  facet_grid(b + a ~ k +  s1, labeller = labeller(
    k = function(x) paste0("k = ", x),
    s1 = function(x) paste0("s1 = ", x),
    a = function(x) paste0("a = ", x),
    b = function(x) paste0("b = ", x))) +
  geom_hline(data=subset(obs_val,stand == "Waldi"), aes(yintercept = N_obs,col = agg_stage), linetype = 2, linewidth = 1) +  
  scale_y_continuous(labels = function(x) format(x, big.mark = "", scientific = FALSE)) +
  theme_bw()
wal_best_param

# check the phenotype for the best parameters values (offspring only)
ggplot(subset(quanti_data, stand == "Waldi" & generation == 75 & stage == 1 & replicate == 2 & k == tab_for_loop_wal$k[1] & s1 == tab_for_loop_wal$s1[1] & b == tab_for_loop_wal$b[1] &
                a == tab_for_loop_wal$a[1]), aes(x = factor(P1))) +
  geom_histogram(stat="count", color = "black",boundary = 0) +
  facet_grid(b + a ~ k+ s1) +
  labs(
    title = "Species phenotype (-1 sylvatica, 1 orientalis)",
    x = "",
    y = "P1"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ALLENWILLER #
tab_for_loop_al=subset(tab_for_loop, stand == "Allenwiller")
tab_for_loop_al = tab_for_loop_al[order(tab_for_loop_al$chi_scaled,decreasing=FALSE),]
head(tab_for_loop_al)

# check hte real numbers of the best model
check <- subset(end_sim_obs2, k == tab_for_loop_al$k[1] & s1 == tab_for_loop_al$s1[1] & b == tab_for_loop_al$b[1] &  a == tab_for_loop_al$a[1])

al_best_param <- ggplot(subset(check, stand == "Allenwiller"), aes(x = agg_stage, y = N_agg_stage, fill = agg_stage)) +
  geom_boxplot() +
  theme_bw() +
  labs(
    title = "",
    x = "",
    y = "Number of Individuals",
    fill = "Simulated data",
    color = "Observed data"
  ) +
  facet_grid(b + a ~ k +  s1) +
  geom_hline(data=subset(obs_val,stand == "Allenwiller"), aes(yintercept = N_obs,col = agg_stage), linetype = 2, linewidth = 1) +  
  scale_y_continuous(labels = function(x) format(x, big.mark = "", scientific = FALSE)) +
  theme_bw()
al_best_param

              # check the phenotype for the best parameters values (offspring only) 
              ggplot(subset(quanti_data, stand == "Allenwiller" & generation == 70 & stage == 1 & replicate == 2 & k == tab_for_loop_al$k[1] & s1 == tab_for_loop_al$s1[1] & b == tab_for_loop_al$b[1] &
                              a == tab_for_loop_al$a[1]), aes(x = factor(P1))) +
                geom_histogram(stat="count", color = "black",boundary = 0) +
                facet_grid(b + a ~ k+ s1) +
                labs(
                  title = "Species phenotype (-1 sylvatica, 1 orientalis)",
                  x = "",
                  y = "P1"
                ) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
              
# merge with data to compare numbers (need to subset for the stand)
end_sim_obs2=merge(end_sim_obs2,tab_for_loop,by=c("b","s1","k","a", "stand"))
head(end_sim_obs2[order(end_sim_obs2$chi_scaled,decreasing=FALSE),])


fwrite(end_sim_obs2, "NEMO/NEMO_demographic_param_test.csv")


