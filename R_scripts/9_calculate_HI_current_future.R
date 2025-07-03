### HI calculation

library(ggplot2)
library(dplyr)
library(tidyr)
library(raster)
library(terra)
library(data.table)
library(readr)
library(stringr)
library(sf)
library(xlsx)
library(ggspatial)
library(ape)
library(ggh4x)
library(cowplot)
library(grid)
library(ggpubr)
library(gridExtra)


# path to get files
setwd(dir="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Data_clean")
#setwd(dir="C:/Users/Camilla/Dropbox (Old (1))/Dropbox/WSL_PhD/Projects/Hybridization/Data_clean")

## plotting theme
theme_p <- theme(
  axis.line = element_blank(),
  plot.background = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  plot.margin = margin(t = 15, r = 2, b = 2, l = 2),
  strip.background = element_blank(),
  strip.placement = "outside",
  strip.text.x = element_text(size = 7, margin = margin(b = 3)),
  strip.text.y.left = element_text(angle = 0, size = 7, margin = margin(r =3)),
  plot.title = element_text(size = 8, hjust = 0.5)
)

# --------- calculate observed HI individual and circle plot -------

# load trees info
al_info <- fread("NewHybrids/Allenwiller_indiv_info_class.csv")
wal_info <- fread("NewHybrids/Waldi_indiv_info_class.csv")

# classify individuals into the three HI categories as in the simulations (remove unassigned)
al_cp_HI <- al_info %>% mutate(HI = case_when(
  str_detect(NH_class_t80, "Pure") ~ 0, 
  NH_class_t80 == "F1" ~ 1,
  str_detect(NH_class_t80, "BC") ~ 0.5,
  NH_class_t80 == "F2" ~ 0.5,
  TRUE ~ NA_integer_
)) %>% filter(!is.na(HI))

wal_cp_HI <- wal_info %>% mutate(HI = case_when(
  str_detect(NH_class_t80_corrected, "Pure") ~ 0, 
  NH_class_t80_corrected == "F1" ~ 1,
  str_detect(NH_class_t80_corrected, "BC") ~ 0.5,
  NH_class_t80_corrected == "F2" ~ 0.5,
  TRUE ~ NA_integer_
))%>% filter(!is.na(HI))

# calculate proportions og HI classes into circle plots
al_cp_HI2 <-  al_cp_HI %>% 
  group_by(CirclePlotID, generation) %>%
  summarize(
    N_HI1 = sum(HI == 1, na.rm = TRUE),
    N_HI0.5 = sum(HI > 0 & HI < 1, na.rm = TRUE),
    N_HI0 = sum(HI == 0, na.rm = TRUE),
    prop_HI1 = sum(HI == 1, na.rm = TRUE) / n(),
    prop_HI0.5 = sum(HI > 0 & HI < 1, na.rm = TRUE) / n(),
    prop_HI0 = sum(HI == 0, na.rm = TRUE) / n(),
    .groups = "drop"
  ) %>%
  mutate(sum_check = prop_HI1 + prop_HI0.5 + prop_HI0)

wal_cp_HI2 <- wal_cp_HI %>%group_by(CirclePlotID, generation) %>%
  summarize(
    N_HI1 = sum(HI == 1),
    N_HI0.5 = sum(HI > 0 & HI < 1),  # Handles any HI strictly between 0 and 1
    N_HI0 = sum(HI == 0),
    prop_HI1 = sum(HI == 1) / n(),
    prop_HI0.5 = sum(HI > 0 & HI < 1) / n(),
    prop_HI0 = sum(HI == 0) / n(),
    .groups = "drop"
  ) %>%
  mutate(sum_check = prop_HI1 + prop_HI0.5 + prop_HI0) 


# merge and save for comparison with simulated
wal_cp_HI2$stand = "Waldi"
al_cp_HI2$stand = "Allenwiller"
observed_cp_HI = bind_rows(wal_cp_HI2, al_cp_HI2)
fwrite(observed_cp_HI, "NEMO/Observed_HI_circleplots.csv")

## check number of seedlings per circle plot against HI
wal_cp_HI
sampling <- st_read("Waldi_sampling_scheme.shp")
colnames(sampling)[2] <- "CirclePlotID"
colnames(sampling)[6] <- "N_seedlings"

indiv_HI_seedlings <- wal_cp_HI %>%
  filter(generation == "Offspring") %>%
  left_join(st_drop_geometry(sampling) %>% dplyr::select(CirclePlotID, N_seedlings), by = "CirclePlotID")

# HI vs seedling numer per circle plot
ggplot(indiv_HI_seedlings, aes(x = N_seedlings, y = HI)) +
  geom_jitter(width = 0.3, height = 0.02) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(
    x = "N seedlings per circle plot",
    y = "Individual HI value",
    title = "Waldi HI ~seedling number"
  ) +
  theme_minimal()

cor.test(HI_seedling_data$N_seedlings, HI_seedling_data$mean_HI)

al_cp_HI
sampling <- st_read("Allenwiller_sampling_scheme.shp")
colnames(sampling)[2] <- "CirclePlotID"
colnames(sampling)[5] <- "N_seedlings"

indiv_HI_seedlings <- al_cp_HI %>%
  filter(generation == "Offspring") %>%
  left_join(st_drop_geometry(sampling) %>% dplyr::select(CirclePlotID, N_seedlings), by = "CirclePlotID")

# HI vs seedling numer per circle plot
ggplot(indiv_HI_seedlings, aes(x = N_seedlings, y = HI)) +
  geom_jitter(width = 0.3, height = 0.02) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(
    x = "N seedlings per circle plot",
    y = "Individual HI value",
    title = "Allenwiller HI ~seedling number"
  ) +
  theme_minimal()

cor.test(HI_seedling_data$N_seedlings, HI_seedling_data$mean_HI)

              ##pie chart visualisations
              # check proportions of hybrids
              hist(observed_cp_HI$prop_HI0.5)
              hist(observed_cp_HI$prop_HI1)
              
              # build scatterpies
              library(scatterpie)
              
              # allenwiller
              sampling <- st_read("Allenwiller_sampling_scheme.shp")
              colnames(sampling)[2] <- "CirclePlotID"
              al_HI_shp <- sampling %>% left_join(al_cp_HI, by = "CirclePlotID")
              
              # centroids cirlce plots to place pie charts
              al_HI_shp_pie <- al_HI_shp %>%
                mutate(
                  centroid = st_centroid(geometry),
                  x = st_coordinates(centroid)[,1],
                  y = st_coordinates(centroid)[,2]
                )
              
              # Create a dataframe suitable for scatterpie
              al_HI_pie_data <- al_HI_shp_pie %>%
                st_drop_geometry() %>%
                dplyr::select(CirclePlotID, x, y, prop_HI0, prop_HI0.5, prop_HI1)
              
              base_map <- ggplot(al_HI_shp) +
                geom_sf(fill = NA, color = "black") +
                theme_void()
              
              pie_chart_map <- base_map +
                geom_scatterpie(
                  data = al_HI_pie_data,
                  aes(x = x, y = y, r = 2),  # Adjust radius as needed
                  cols = c("prop_HI0", "prop_HI0.5", "prop_HI1"),
                  color = "black",
                  alpha = 0.7
                ) +
                scale_fill_manual(
                  values = c("prop_HI0" = "#CCCCCC", "prop_HI0.5" = "#3CB371", "prop_HI1" = "#FF6347"),
                  labels = c("HI=0", "HI=0.5", "HI=1"),
                  name = "HI Value"
                ) 
              pie_chart_map
              
              
              # allenwiller
              sampling <- st_read("Waldi_sampling_scheme.shp")
              colnames(sampling)[2] <- "CirclePlotID"
              wal_HI_shp <- sampling %>% left_join(wal_cp_HI, by = "CirclePlotID")
              
              # Calculate centroids of each polygon for pie chart placement
              wal_HI_shp_pie <- wal_HI_shp %>%
                mutate(
                  centroid = st_centroid(geometry),
                  x = st_coordinates(centroid)[,1],
                  y = st_coordinates(centroid)[,2]
                )
              
              # Create a dataframe suitable for scatterpie
              wal_HI_pie_data <- wal_HI_shp_pie %>%
                st_drop_geometry() %>%
                dplyr::select(CirclePlotID, x, y, prop_HI0, prop_HI0.5, prop_HI1)
              
              # Create a spatial plot with pie charts
              # First, create the base map with circle plot boundaries
              base_map <- ggplot(wal_HI_shp) +
                geom_sf(fill = NA, color = "black") +
                theme_void()
              
              pie_chart_map <- base_map +
                geom_scatterpie(
                  data = wal_HI_pie_data,
                  aes(x = x, y = y, r = 4),  # Adjust radius as needed
                  cols = c("prop_HI0", "prop_HI0.5", "prop_HI1"),
                  color = "black",
                  alpha = 0.7
                ) +
                scale_fill_manual(
                  values = c("prop_HI0" = "#CCCCCC", "prop_HI0.5" = "#3CB371", "prop_HI1" = "#FF6347"),
                  labels = c("HI=0", "HI=0.5", "HI=1"),
                  name = "HI Value"
                ) 
              pie_chart_map


            ### observed Moran's I 
            
            # filter only the observed circle plots with data (some have no individuals)
            wal_HI <- wal_HI[!(is.na(wal_HI$HI)),]
            # get the cp ID for the simulated data
            wal_valid_cp <- wal_HI$CirclePlotID
            coords_indivs <- as.matrix(wal_HI[, c("x", "y")])
            # distance matrix coordinates individuals
            dist_matrix <- as.matrix(dist(coords_indivs))
            # take inverse of matrix values and set diagonal to zero (self-neighbors)
            inv_dist_matrix <- 1 / dist_matrix
            inv_dist_matrix[is.infinite(inv_dist_matrix)] <- 0 # for overlapping positions
            diag(inv_dist_matrix) <- 0 
            wal_HI <- wal_HI[!(is.na(wal_HI$HI)),]
            moran_result <- Moran.I(wal_HI$HI, inv_dist_matrix)
            
            moran_waldi <- data.frame(
              stand = "Waldi",
              morans_I = moran_result$observed,
              p_value = moran_result$p.value
            )
            
            # with spdep?
            
            
            ## repeat for allenwiller
            # filter only the observed circle plots with data (some have no individuals)
            al_HI <- al_HI[!(is.na(al_HI$HI)),]
            # get the cp ID for the simulated data
            al_valid_cp <- al_HI$CirclePlotID
            coords_indivs <- as.matrix(al_HI[, c("x", "y")])
            # distance matrix coordinates individuals
            dist_matrix <- as.matrix(dist(coords_indivs))
            # take inverse of matrix values and set diagonal to zero (self-neighbors)
            inv_dist_matrix <- 1 / dist_matrix
            inv_dist_matrix[is.infinite(inv_dist_matrix)] <- 0 # for overlapping positions
            diag(inv_dist_matrix) <- 0 
            al_HI <- al_HI[!(is.na(al_HI$HI)),]
            moran_result <- Moran.I(al_HI$HI, inv_dist_matrix)
            
            moran_allen <- data.frame(
              stand = "Allenwiller",
              morans_I = moran_result$observed,
              p_value = moran_result$p.value
            )
            
            
            moran_observed <- bind_rows(moran_waldi, moran_allen)
            moran_observed


# --------- calculate individual HI simulations (TO RUN IN THE CLUSTER) -----------

## load NEMO simulation results  (tr only with current data to speed up)
quanti_data_allen_current <- readRDS("NEMO/Quanti_data_Allenwiller_year2020.RDS")
quanti_data_waldi_current <- readRDS("Quanti_data_Waldi_year2020.RDS")
quanti_data_waldi_current <- quanti_data_waldi_current[, -c(24:31, 39:44)]

# col names for the loci
loci <- unique(gsub("(t1l[0-9]+)[xy]$", "\\1", grep("^t1l.*[xy]$", names(quanti_data_waldi_current), value = TRUE)))

quanti_data_waldi_current[, HI := sapply(1:nrow(quanti_data_waldi_current), function(i) {
  #count heterozygous loci for eahc individual
  hetero_count <- 0
  
  for (locus in loci) {
    x_col <- paste0(locus, "x")
    y_col <- paste0(locus, "y")
    
    # heterozygous if one locus is negative and other is positive
    x_val <- quanti_data_waldi_current[i, get(x_col)]
    y_val <- quanti_data_waldi_current[i, get(y_col)]
    
    if ((x_val > 0 && y_val < 0) || (x_val < 0 && y_val > 0)) {
      hetero_count <- hetero_count + 1
    }
  }
  
  # proportion of heterozygous loci
  return(hetero_count / length(loci))
})]

saveRDS(quanti_data_waldi_current, "Quanti_data_waldi_current_HI.RDS")

### same for allenwiller




# --------- Allenwiller remove buffering edge -------
## 2020
data <-  readRDS("NEMO/Quanti_data_allen_current_HI.RDS")
data <- data[, c("pop", "stage", "stand", "scenario", "dp", "replicate", "s3_age", "HI")]

grid <- vect("Allenwiller_grid_r4m_9844patches.shp")
sampling <- vect("Allenwiller_sampling_scheme.shp")
names(sampling)[2] <- "CirclePlotID"

# sampling plots
sampling_hull <- terra::convHull(sampling)  # convex hull around sampling points
buffer_sampling <- terra::buffer(sampling_hull, width = 10)  # buffer

# simulation triangle adults
coords <- rbind(
  c(4126700, 2840045),
  c(4126870, 2840200),
  c(4126705, 2840260),
  c(4126700, 2840045)  # Close the polygon
)
triangle <- vect(coords, type = "polygons", crs = "EPSG:3035")
triangle <- terra::buffer(triangle, width = 10) 

# combine geometries
triangle <- project(triangle, crs(buffer_sampling))  # just to be safe
combined_area <- rbind(buffer_sampling, triangle)
combined_area <- terra::aggregate(combined_area)
# intersect polygon with grid
grid_subset <- terra::intersect(grid, combined_area)

# filter data by patch IDs in the subset area
patchID_in_poly <- grid_subset$patch.ID
data_filtered <- data %>% filter(pop %in% patchID_in_poly)

plot(grid)
plot(combined_area, add = TRUE, col = "red")
plot(sampling, add = TRUE, col = "black")

saveRDS(data_filtered, "NEMO/Quanti_data_allen_current_HI_cut.RDS")

## 2120
data <-  readRDS("NEMO/Quanti_data_allen_future_HI.RDS")
data <- data[, c("pop", "stage", "stand", "scenario", "dp","replicate","s3_age","HI" )]
data_filtered <- data %>% filter(pop %in% patchID_in_poly)
saveRDS(data_filtered, "NEMO/Quanti_data_allen_future_HI_cut.RDS")

# --------- count HI per scenario -----
sim_HI_files_2020 <- c("NEMO/Quanti_data_waldi_current_HI.RDS", "NEMO/Quanti_data_allen_current_HI_cut.RDS")
sim_HI_files_2120 <- c("NEMO/Quanti_data_waldi_future_HI.RDS", "NEMO/Quanti_data_allen_future_HI_cut.RDS")

count_HI_scenario <- function(sim_HI_files) {
  sim_HI_data <- lapply(sim_HI_files, function(file) {
    readRDS(file) %>%
      mutate(scenario2 = paste0(scenario, "_a", s3_age)) %>%
      group_by(stand, scenario2, replicate,stage) %>%
      summarize(
        count_HI1 = sum(HI == 1),
        count_HI0.5 = sum(HI > 0 & HI < 1),
        count_HI0 = sum(HI == 0),
        prop_HI1 = sum(HI == 1) / n(),
        prop_HI0.5 = sum(HI > 0 & HI < 1) / n(),
        prop_HI0 = sum(HI == 0) / n(),
        .groups = "drop"
      ) %>%
      ungroup()
  })%>%
    bind_rows() %>%
    group_by(stand, scenario2, stage) %>%
    summarise(
      mean_count_HI1 = mean(count_HI1, na.rm = TRUE),
      sd_count_HI1 = sd(count_HI1, na.rm = TRUE),
      mean_count_HI0.5 = mean(count_HI0.5, na.rm = TRUE),
      sd_count_HI0.5 = sd(count_HI0.5, na.rm = TRUE),
      mean_count_HI0 = mean(count_HI0, na.rm = TRUE),
      sd_count_HI0 = sd(count_HI0, na.rm = TRUE),
      mean_prop_HI1 = mean(prop_HI1, na.rm = TRUE),
      sd_prop_HI1 = sd(prop_HI1, na.rm = TRUE),
      mean_prop_HI0.5 = mean(prop_HI0.5, na.rm = TRUE),
      sd_prop_HI0.5 = sd(prop_HI0.5, na.rm = TRUE),
      mean_prop_HI0 = mean(prop_HI0, na.rm = TRUE),
      sd_prop_HI0 = sd(prop_HI0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      total_count_mean = mean_count_HI1+ mean_count_HI0.5 + mean_count_HI0,
      total_prop_mean = mean_prop_HI1 + mean_prop_HI0.5 + mean_prop_HI0, 
      hybridization_rate = mean_prop_HI1 + mean_prop_HI0.5)
}

sim_HI_2020_summary <- count_HI_scenario(sim_HI_files_2020)
sim_HI_2120_summary <- count_HI_scenario(sim_HI_files_2120)

# visualise
long_data_2020 <- sim_HI_2020_summary %>%
  dplyr::select(stand, scenario2, stage,
         mean_prop_HI1, sd_prop_HI1,
         mean_prop_HI0.5, sd_prop_HI0.5,
         mean_prop_HI0, sd_prop_HI0) %>%
  pivot_longer(
    cols = starts_with("mean_prop_HI"),
    names_to = "HI",
    names_prefix = "mean_prop_",
    values_to = "mean_prop"
  ) %>%
  mutate(
    sd = case_when(
     HI == "HI1" ~ sd_prop_HI1,
     HI == "HI0.5" ~ sd_prop_HI0.5,
     HI == "HI0" ~ sd_prop_HI0
    )
  )

HI_comp_2020<- ggplot(subset(long_data_2020, HI == "HI1"| HI == "HI0.5"), aes(x = factor(scenario2), y = mean_prop, fill = HI)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_prop - sd, ymax = mean_prop + sd),
                position = position_dodge(width = 0.8),
                width = 0.2) +
  ggh4x::facet_grid2(stand ~ stage, labeller = labeller(stage =  c("1" = "Offspring","2" = "Juvenile","3" = "Adults" )) )+
  labs(x = "Stage", y = "Mean Proportion", fill = "HI", title = "Total proportions of HI1 and HI0.5 in 2020") +
  theme_bw() +
  theme(strip.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))
HI_comp_2020

long_data_2120 <- sim_HI_2120_summary %>%
  dplyr::select(stand, scenario2, stage,
                mean_prop_HI1, sd_prop_HI1,
                mean_prop_HI0.5, sd_prop_HI0.5,
                mean_prop_HI0, sd_prop_HI0) %>%
  pivot_longer(
    cols = starts_with("mean_prop_HI"),
    names_to = "HI",
    names_prefix = "mean_prop_",
    values_to = "mean_prop"
  ) %>%
  mutate(
    sd = case_when(
      HI == "HI1" ~ sd_prop_HI1,
      HI == "HI0.5" ~ sd_prop_HI0.5,
      HI == "HI0" ~ sd_prop_HI0
    )
  )

HI_comp_2120<- ggplot(subset(long_data_2120, HI == "HI1"| HI == "HI0.5"), aes(x = factor(scenario2), y = mean_prop, fill = HI)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_prop - sd, ymax = mean_prop + sd),
                position = position_dodge(width = 0.8),
                width = 0.2) +
  ggh4x::facet_grid2(stand ~ stage, labeller = labeller(stage = c("1" = "Offspring","2" = "Juvenile","3" = "Adults" )) )+
  labs(x = "Stage", y = "Mean Proportion", fill = "HI", title = "Total proportions of HI1 and HI0.5 in 2120") +
  theme_bw() +
  theme(strip.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))
HI_comp_2120

HI_comp <- ggarrange(HI_comp_2020, HI_comp_2120, ncol = 1, nrow = 2, common.legend = T)

# save in NEMO folder 
path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "NEMO/Figures/HI_scenario_comparison.png"), 
    width =15, height=17, units = "in", res=600)
plot(HI_comp)
dev.off()

# save for nex tplots
saveRDS(sim_HI_2020_summary, "NEMO/NEMO_HI_2020_summary.RDS")
saveRDS(sim_HI_2120_summary, "NEMO/NEMO_HI_2120_summary.RDS")


# --------- calculate patch HI and plot across grid (adults and offspring) -------
process_HI_data <- function(sim_HI_files, grid_files, stand_names, stage_filter) {
  sim_HI_data <- lapply(sim_HI_files, function(file) {
    readRDS(file) %>%
      filter(stage == stage_filter) %>%
      mutate(scenario2 = paste0(scenario, "_a", s3_age)) %>%
      group_by(stand, scenario2, replicate, pop) %>%
      summarize(
        prop_HI1 = sum(HI == 1) / n(),
        prop_HI0.5 = sum(HI > 0 & HI < 1) / n(),
        prop_HI0 = sum(HI == 0) / n(),
        .groups = "drop"
      ) %>%
      ungroup()
  }) %>%
    bind_rows() %>%
    group_by(stand, pop, scenario2) %>%
    summarise(
      mean_prop_HI1 = mean(prop_HI1, na.rm = TRUE),
      sd_prop_HI1 = sd(prop_HI1, na.rm = TRUE),
      mean_prop_HI0.5 = mean(prop_HI0.5, na.rm = TRUE),
      sd_prop_HI0.5 = sd(prop_HI0.5, na.rm = TRUE),
      mean_prop_HI0 = mean(prop_HI0, na.rm = TRUE),
      sd_prop_HI0 = sd(prop_HI0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(total_mean = mean_prop_HI1 + mean_prop_HI0.5 + mean_prop_HI0)
  
  # merge with grid
  grid_HI_data <- lapply(seq_along(grid_files), function(i) {
    grid <- read_sf(grid_files[i])
    stand_data <- filter(sim_HI_data, stand == stand_names[i])
    grid %>%
      left_join(stand_data, by = c("patch.ID" = "pop")) %>%
      filter(!is.na(scenario2))
  }) %>%
    bind_rows()
  
  grid_HI_data %>%
    mutate(
      Ar = as.integer(str_extract(scenario2, "_a\\d+$") %>% str_remove("_a")),
      dp = case_when(
        str_detect(scenario2, "dp0") ~ "Flat\ndispersal",
        str_detect(scenario2, "dp\\d+") ~ str_extract(scenario2, "dp\\d+") %>% str_remove("dp"),
        TRUE ~ NA_character_
      ),
      s = case_when(
        str_detect(scenario2, "s[0-9]\\.[0-9]") ~ as.numeric(str_extract(scenario2, "s[0-9]\\.[0-9]") %>% str_remove("s")),
        TRUE ~ 0
      ),
      dp = factor(dp, levels = c("Flat\ndispersal", "109", "50", "20")),
      s = factor(s, levels = sort(unique(na.omit(s))))
    ) %>%
    mutate(
      dominant_HI = case_when(
        mean_prop_HI1 >= mean_prop_HI0.5 & mean_prop_HI1 >= mean_prop_HI0 ~ "HI1",
        mean_prop_HI0.5 >= mean_prop_HI1 & mean_prop_HI0.5 >= mean_prop_HI0 ~ "HI0.5",
        TRUE ~ "HI0"
      ),
      dominance_strength = case_when(
        dominant_HI == "HI1" ~ mean_prop_HI1,
        dominant_HI == "HI0.5" ~ mean_prop_HI0.5,
        TRUE ~ mean_prop_HI0
      )
    )
}

# grids and stand names
grid_files <- c("Waldi_grid_4m_5396patches.shp", "Allenwiller_grid_r4m_9844patches.shp")
stand_names <- c("Waldi", "Allenwiller")

# process adult data
sim_HI_files_2020 <- c("NEMO/Quanti_data_waldi_current_HI.RDS", "NEMO/Quanti_data_allen_current_HI_cut.RDS")
sim_HI_files_2120 <- c("NEMO/Quanti_data_waldi_future_HI.RDS", "NEMO/Quanti_data_allen_future_HI_cut.RDS")

grid_adults_2020 <- process_HI_data(sim_HI_files_2020, grid_files, stand_names, stage_filter = 3)
grid_adults_2120 <- process_HI_data(sim_HI_files_2120, grid_files, stand_names, stage_filter = 3)

# process offspring data 
grid_off_2020 <- process_HI_data(sim_HI_files_2020, grid_files, stand_names, stage_filter = 1)
grid_off_2120 <- process_HI_data(sim_HI_files_2120, grid_files, stand_names, stage_filter = 1)



# save for figure
saveRDS(grid_adults_2020, "NEMO/HI_data_adults_2020.RDS")
saveRDS(grid_adults_2120, "NEMO/HI_data_adults_2120.RDS")
saveRDS(grid_off_2020, "NEMO/HI_data_offspring_2020.RDS")
saveRDS(grid_off_2120, "NEMO/HI_data_offspring_2120.RDS")

plot_HI <- function(data_subset, title, sampling_data = NULL) {
  p <- ggplot(data_subset) +
    geom_sf(aes(fill = dominant_HI, alpha = dominance_strength), color = NA) +
    scale_fill_manual(
      values = c("HI0" = "#FDE725FF", "HI0.5" = "#440154FF", "HI1" = "#25858EFF"),
      name = "HI",
      labels = c("HI = 0", "0 < HI < 1", "HI = 1")
    ) +
    scale_alpha_continuous(range = c(0.2, 1)) +
    ggh4x::facet_grid2(dp ~ s, labeller = label_value, switch = "both", render_empty = FALSE) +
    theme_minimal(base_size = 10) +
    guides(alpha = "none") +
    labs(title = title, x = "Selection strength (s)", y = "Dispersal distance (dp)") +
    theme_p
  
  if (!is.null(sampling_data)) {
    p <- p + geom_sf(data = sampling_data, shape = 21, color = "grey30", fill = NA, size = 0.5, stroke = 0.3, alpha = 0.6)
  }
  
  return(p)
}

## add sampling scheme to the offspring
samplings <- c("Waldi_sampling_scheme.shp", "Allenwiller_sampling_scheme.shp")
sampling_schemes <- list()
for (i in seq_along(samplings)) {
  sampling_schemes[[i]] <- read_sf(samplings[i]) %>%
    mutate(stand = stand_names[i])
  sampling_schemes[[i]] <- sampling_schemes[[i]][,-c(3:13)] 
}
sampling_all <- bind_rows(sampling_schemes)

# adults plots 2020
wal_adults_30 <- plot_HI(subset(grid_adults_2020, stand == "Waldi" & Ar == 30), "Ar = 30")
wal_adults_50 <- plot_HI(subset(grid_adults_2020, stand == "Waldi" & Ar == 50), "Ar = 50")
al_adults_30 <- plot_HI(subset(grid_adults_2020, stand == "Allenwiller" & Ar == 30), "Ar = 30")
al_adults_50 <- plot_HI(subset(grid_adults_2020, stand == "Allenwiller" & Ar == 50), "Ar = 50")

legend_HI <- cowplot::get_legend(
  wal_adults_30 +
    theme(
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 6), 
      legend.key.size = unit(0.3, "cm")
    ) +
    guides(fill = guide_legend(direction = "vertical"))
)
grid::grid.newpage()
grid::grid.draw(legend_HI)


Figure_adults_sim_HI_grids_2020 <- ggpubr::ggarrange(
  al_adults_30+ theme(axis.title.x = element_blank(), strip.text.x =element_blank()),
  al_adults_50 + theme(axis.title = element_blank(),strip.text.x =element_blank(), strip.text.y.left = element_blank()),
  wal_adults_30, 
  wal_adults_50 + theme(axis.title.y=element_blank(),  strip.text.y.left = element_blank() ),
  nrow = 2, ncol = 2,
  widths = c(1, 1), heights = c(1, 1),
  labels = c("Allenwiller", "", "Wäldi", ""),
  label.x = 0.1,  label.y = 1,
  font.label = list(size = 12, face = "bold"),
  align = "hv"
)  

Figure_adults_sim_HI_grids_2020
legend <- cowplot::get_legend(al_adults_30)

path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/Figure_adults_sim_HI_grids_2020.png"), 
    width =7.0866, height=7.0866, units = "in", res=600)
plot(Figure_adults_sim_HI_grids_2020)
dev.off()

# adults plots 2120
wal_adults_30 <- plot_HI(subset(grid_adults_2120, stand == "Waldi" & Ar == 30), "Ar = 30")
wal_adults_50 <- plot_HI(subset(grid_adults_2120, stand == "Waldi" & Ar == 50), "Ar = 50")
al_adults_30 <- plot_HI(subset(grid_adults_2120, stand == "Allenwiller" & Ar == 30), "Ar = 30")
al_adults_50 <- plot_HI(subset(grid_adults_2120, stand == "Allenwiller" & Ar == 50), "Ar = 50")

Figure_adults_sim_HI_grids_2120 <- ggpubr::ggarrange(
  al_adults_30+ theme(legend.position ="none", axis.title.x = element_blank(), strip.text.x =element_blank()),
  al_adults_50 + theme(legend.position ="none",axis.title = element_blank(),strip.text.x =element_blank(), strip.text.y.left = element_blank()),
  wal_adults_30 + theme(legend.position ="none"), 
  wal_adults_50 + theme(legend.position ="none",axis.title.y=element_blank(),  strip.text.y.left = element_blank() ),
  nrow = 2, ncol = 2,
  widths = c(1, 1), heights = c(1, 1),
  labels = c("Allenwiller", "", "Wäldi", ""),
  label.x = 0.1,  label.y = 1,
  font.label = list(size = 12, face = "bold"),
  align = "hv"
)  

Figure_adults_sim_HI_grids_2120
#legend <- cowplot::get_legend(al_adults_30)

path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/Figure_adults_sim_HI_grids_2120.png"), 
    width =7.0866, height=7.0866, units = "in", res=600)
plot(Figure_adults_sim_HI_grids_2120)
dev.off()

# offsrping plots 2020 
wal_off_30 <- plot_HI(subset(grid_off_2020, stand == "Waldi" & Ar == 30),"Ar = 30",sampling_data = sampling_all %>% filter(stand == "Waldi"))
wal_off_50 <- plot_HI(subset(grid_off_2020, stand == "Waldi" & Ar == 50), "Ar = 50",sampling_data = sampling_all %>% filter(stand == "Waldi"))
al_off_30 <- plot_HI(subset(grid_off_2020, stand == "Allenwiller" & Ar == 30), "Ar = 30",sampling_data = sampling_all %>% filter(stand == "Allenwiller"))
al_off_50 <- plot_HI(subset(grid_off_2020, stand == "Allenwiller" & Ar == 50), "Ar = 50",sampling_data = sampling_all %>% filter(stand == "Allenwiller"))

Figure_off_sim_HI_grids_2020 <- ggpubr::ggarrange(
  al_off_30+ theme(axis.title.x = element_blank(), strip.text.x =element_blank(),legend.position ="none"), 
  al_off_50+ theme(axis.title = element_blank(),strip.text.x =element_blank(), strip.text.y.left = element_blank(),legend.position ="none"),
  wal_off_30 + theme(legend.position ="none"),
  wal_off_50+ theme(axis.title.y=element_blank(),  strip.text.y.left = element_blank(),legend.position ="none" ),
  nrow = 2, ncol = 2,
  widths = c(1, 1), heights = c(1, 1),
  labels = c("Allenwiller", "", "Wäldi", ""),
  label.x = 0.1,  label.y = 1,
  font.label = list(size = 12, face = "bold"),
  align = "hv"
)  
  
Figure_off_sim_HI_grids_2020

path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/Figure_offs_sim_HI_grids_2020.png"), 
    width =7.0866, height=7.0866, units = "in", res=600)
plot(Figure_off_sim_HI_grids_2020)
dev.off()

# offsrping plots 2120
wal_off_30 <- plot_HI(subset(grid_2120, stand == "Waldi" & Ar == 30),"Ar = 30",sampling_data = sampling_all %>% filter(stand == "Waldi"))
wal_off_50 <- plot_HI(subset(grid_2120, stand == "Waldi" & Ar == 50), "Ar = 50",sampling_data = sampling_all %>% filter(stand == "Waldi"))
al_off_30 <- plot_HI(subset(grid_2120, stand == "Allenwiller" & Ar == 30), "Ar = 30",sampling_data = sampling_all %>% filter(stand == "Allenwiller"))
al_off_50 <- plot_HI(subset(grid_2120, stand == "Allenwiller" & Ar == 50), "Ar = 50",sampling_data = sampling_all %>% filter(stand == "Allenwiller"))

Figure_off_sim_HI_grids_2120 <- ggpubr::ggarrange(
  al_off_30+ theme(plot.margin = unit(c(1, 0, 1, 0), "mm"), axis.title.x = element_blank(), strip.text.x =element_blank(),legend.position ="none"), 
  al_off_50+ theme(plot.margin = unit(c(1, 0, 1, 0), "mm"),axis.title = element_blank(),strip.text.x =element_blank(), strip.text.y.left = element_blank(),legend.position ="none"),
  wal_off_30 + theme(plot.margin = unit(c(1, 0, 1, 0), "mm"),legend.position ="none"),
  wal_off_50+ theme(plot.margin = unit(c(1, 0, 1, 0), "mm"),axis.title.y=element_blank(),  strip.text.y.left = element_blank(),legend.position ="none" ),
  nrow = 2, ncol = 2,
  widths = c(1, 1), heights = c(1, 1),
  labels = c("Allenwiller", "", "Wäldi", ""),
  label.x = 0.1,  label.y = 1,
  font.label = list(size = 12, face = "bold"),
  align = "hv"
)  

Figure_off_sim_HI_grids_2120
path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/Figure_offs_sim_HI_grids_2120.png"), 
    width =7.0866, height=7.0866, units = "in", res=600)
plot(Figure_off_sim_HI_grids_2120)
dev.off()

# --------- calculate simulated HI per circle plot 2020 ----------
grids <- c("Waldi_grid_4m_5396patches.shp","Allenwiller_grid_r4m_9844patches.shp" )
samplings <- c("Waldi_sampling_scheme.shp", "Allenwiller_sampling_scheme.shp")
sim_HI_files <- c("NEMO/Quanti_data_waldi_current_HI.RDS", "NEMO/Quanti_data_allen_current_HI_cut.RDS")

sim_HI_data <- list()
# again dplyr not working
for (i in seq_along(sim_HI_files)) {
  sim_HI <- readRDS(sim_HI_files[i]) %>%
    filter(stage == 1) %>%
    mutate(
      scenario2 = paste0(scenario, "_a", s3_age),
      HI_class = case_when(
        HI == 1 ~ "HI1",
        HI == 0 ~ "HI0",
        HI > 0 & HI < 1 ~ "HI0.5",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(HI_class)) %>%
    group_by(stand, scenario2, replicate, pop, HI_class) %>%
    dplyr::summarise(n = n(), .groups = "drop")%>%
    tidyr::pivot_wider(
      names_from = HI_class,
      values_from = n,
      values_fill = list(n = 0)
    ) %>%
    dplyr::rename(
      n_HI1 = HI1,
      n_HI0.5 = HI0.5,
      n_HI0 = HI0
    ) %>%
    mutate(n_total = n_HI1 + n_HI0.5 + n_HI0)
  
  sim_HI_data[[i]] <- sim_HI
}

# check if there are empty patches 
check <-sim_HI_data[[1]]
hist(check$n_total)
subset(sim_HI_data[[1]], n_HI1== 0)
subset(sim_HI_data[[2]], n_HI1 == 0)

# calculate the circle plot HI as sum of HI of patches overlapping with cp, weigthed by the respective overallping area
result_list <- list()

for (i in seq_along(grids)) {
  grid <- vect(grids[i])
  sampling <- vect(samplings[i])
  names(sampling)[2]<- "CirclePlotID"
  sampling <- sampling[!is.na(geom(sampling)[, "x"]),]
  grid_HI <- sim_HI_data[[i]]
  # get combinations to process
  combinations <- grid_HI %>%
    dplyr::select(stand, scenario2, replicate) %>%
    distinct()
  # filter combination to process
  for (j in seq_len(nrow(combinations))) {
    combo <- combinations[j, ]
    combo_HI <- grid_HI %>%
      filter(stand == combo$stand, scenario2 == combo$scenario2, replicate == combo$replicate)
    # attach counts to the grid
    grid$n_HI1 <- combo_HI$n_HI1[match(grid$patch.ID, combo_HI$pop)]
    grid$n_HI0.5 <- combo_HI$n_HI0.5[match(grid$patch.ID, combo_HI$pop)]
    grid$n_HI0 <- combo_HI$n_HI0[match(grid$patch.ID, combo_HI$pop)]
    grid$n_total <- combo_HI$n_total[match(grid$patch.ID, combo_HI$pop)]
    # area of intersection between the sampled circle plots and the grid
    area_match <- terra::intersect(grid, sampling)
    area_match$intersection_area <- terra::expanse(area_match)
    
    # weigth the contribution of each intersecting patch on its intersection area with the circle plot
    area_match1 <- area_match %>%
      as.data.frame() %>%
      group_by(CirclePlotID) %>%
      mutate(weight = intersection_area / sum(intersection_area)) %>%
      summarise(
        N_HI1 = round(sum(n_HI1 * weight, na.rm = TRUE)),
        N_HI0.5 = round(sum(n_HI0.5 * weight, na.rm = TRUE)),
        N_HI0 = round(sum(n_HI0 * weight, na.rm = TRUE)),
        stand = combo$stand,
        scenario2 = combo$scenario2,
        replicate = combo$replicate,
        .groups = "drop"
      ) %>%
      mutate(
        N_total = N_HI1 + N_HI0.5 + N_HI0,
        prop_HI1 = N_HI1 / N_total,
        prop_HI0.5 = N_HI0.5 / N_total,
        prop_HI0 = N_HI0 / N_total,
        sum_check = prop_HI1 + prop_HI0.5 + prop_HI0 
      )    
    result_list[[paste0("Stand_", i, "_Comb_", j)]] <- area_match1
  }
}

sim_HI_cp <- bind_rows(result_list)

# empty circle plots -- convert NA to 0
subset(sim_HI_cp, N_total == 0)
sim_HI_cp[is.na(sim_HI_cp$N_total), ]
sim_HI_cp_clean <- sim_HI_cp %>%
  mutate(
    prop_HI1    = if_else(N_total == 0, 0, prop_HI1),
    prop_HI0.5  = if_else(N_total == 0, 0, prop_HI0.5),
    prop_HI0    = if_else(N_total == 0, 0, prop_HI0),
    sum_check   = if_else(N_total == 0, 0, sum_check)
  )

sim_HI_cp_clean[is.na(sim_HI_cp_clean$prop_HI1), ]
sim_HI_cp_clean[is.na(sim_HI_cp_clean$sum_check), ]

saveRDS(sim_HI_cp_clean, "NEMO/Simulated_HI_counts_circleplots_2020.RDS")

# --------- plot difference P_hi(SIM - OBS) in the circle plots -----------
sim_HI <- readRDS("NEMO/Simulated_HI_counts_circleplots_2020.RDS")
obs_HI <- fread("NEMO/Observed_HI_circleplots.csv")

#avergae across replicates 
mean_sim_HI <- sim_HI %>% group_by(scenario2, stand, CirclePlotID) %>%
  summarize(prop_HI1 = mean(prop_HI1, na.rm = T), 
            prop_HI0.5 = mean(prop_HI0.5, na.rm = T),
            prop_HI0 = mean(prop_HI0, na.rm = T))

# merge and calculate difference
obs_HI <- obs_HI[obs_HI$generation=="Offspring"]
merged_HI <- obs_HI %>%  inner_join(mean_sim_HI, by = c("CirclePlotID", "stand"), suffix = c("_obs", "_sim"))

merged_HI <- merged_HI %>% mutate(diff_prop_HI1 = prop_HI1_sim - prop_HI1_obs,
                                  diff_prop_HI0.5 = prop_HI0.5_sim - prop_HI0.5_obs,
                                  diff_prop_HI0 = prop_HI0_sim - prop_HI0_obs)

# attach values to the sampling shapefile
samplings <- c("Waldi_sampling_scheme.shp", "Allenwiller_sampling_scheme.shp")
stand_names <- c("Waldi", "Allenwiller")

HI_data <- list()
for (i in seq_along(samplings)) {
  sampling <- read_sf(samplings[i])
  colnames(sampling)[2]<- "CirclePlotID"
  sampling <- sampling[,-c(3:13)] 
  stand_data <- merged_HI %>% filter(stand == stand_names[i])
  HI <- sampling %>% left_join(stand_data, by = "CirclePlotID") %>%
    filter(!is.na(scenario2))
  
  HI_data[[i]] <- HI
}

HI_data_grid <- bind_rows(HI_data)

HI_data_grid <- HI_data_grid %>%
  mutate(
    Ar = as.integer(str_extract(scenario2, "_a\\d+$") %>% str_remove("_a")),
    dp = case_when(
      str_detect(scenario2, "dp0") ~ "Flat\ndispersal", 
      str_detect(scenario2, "dp\\d+") ~ str_extract(scenario2, "dp\\d+") %>% str_remove("dp"),
      TRUE ~ NA_character_
    ),
    s = case_when(
      str_detect(scenario2, "s[0-9]\\.[0-9]") ~ as.numeric(str_extract(scenario2, "s[0-9]\\.[0-9]") %>% str_remove("s")),
      TRUE ~ 0
    ),
    dp = factor(dp, levels = c("Flat\ndispersal", "109", "50", "20")),
    s = factor(s, levels = sort(unique(na.omit(s))))  
  )

# read full files for background
sim_HI_files_2020 <- c("NEMO/Quanti_data_waldi_current_HI.RDS", "NEMO/Quanti_data_allen_current_HI_cut.RDS")
stand_names <- c("Waldi", "Allenwiller")
grid_files <- c("Waldi_grid_4m_5396patches.shp", "Allenwiller_grid_r4m_9844patches.shp")

plot_HI_diff <- function(data_subset, title, diff_prop, HI_value, zoom = TRUE) {
  stand_name <- unique(data_subset$stand)
  
  data_subset <- data_subset %>%
    mutate(
      dp_label = as.character(dp),
      dp_label = ifelse(dp_label == "1000", "Flat\ndispersal", dp_label),
      dp_label = factor(dp_label, levels = c("Flat\ndispersal", "20", "50", "109"))
    )
  
  base_plot <- ggplot() +
    geom_sf(
      data = background_grid %>%
        filter(stand == stand_name, has_data == TRUE),
      fill = "grey90", color = NA
    ) +
    geom_sf(data = data_subset, aes(fill = {{diff_prop}}), color = NA) +
    scale_fill_gradient2(
      name = bquote(p[.(HI_value)]^sim - p[.(HI_value)]^obs),
      low = "#00365F",
      mid = "#FEFDBF",
      high = "#680A1C",
      midpoint = 0,
      limits = c(-1, 1),
      breaks = c(-1, -0.5, 0, 0.5, 1),
      oob = scales::rescale_none
    ) +
    ggh4x::facet_grid2(dp ~ s, labeller = label_value, switch = "both", render_empty = FALSE) +
    theme_minimal(base_size = 10) +
    guides(alpha = "none") +
    labs(
      title = title,
      x = "Selection strength (s)",
      y = "Dispersal distance (dp)"
    ) +
    theme_p
  
  # Zoom in for Allenwiller only
  if (stand_name == "Allenwiller" && zoom) {
    data_bbox <- sf::st_bbox(data_subset)
    # little padding around
    padding_x <- 20
    padding_y <- 20
    
    base_plot <- base_plot +
      coord_sf(
        xlim = c(data_bbox["xmin"] - padding_x, data_bbox["xmax"] + padding_x),
        ylim = c(data_bbox["ymin"] - padding_y, data_bbox["ymax"] + padding_y),
        expand = FALSE
      )
  } else {
    # coord_sf without zooming for Waldi to avoid distortion
    base_plot <- base_plot +
      coord_sf(expand = FALSE)
  }
  
  return(base_plot)
}

wal_cp_30 <- plot_HI_diff(subset(HI_data_grid, stand == "Waldi" & Ar == 30), "Ar = 30", diff_prop_HI1, 1, zoom = T)
wal_cp_50 <- plot_HI_diff(subset(HI_data_grid, stand == "Waldi" & Ar == 50), "Ar = 50", diff_prop_HI1, 1, zoom = T)
al_cp_30 <- plot_HI_diff(subset(HI_data_grid, stand == "Allenwiller" & Ar == 30), "Ar = 30", diff_prop_HI1, 1, zoom = T)
al_cp_50 <- plot_HI_diff(subset(HI_data_grid, stand == "Allenwiller" & Ar == 50), "Ar = 50", diff_prop_HI1, 1, zoom = T)

HI_diff_plot <- ggpubr::ggarrange(
  al_cp_30+ theme(axis.title.x = element_blank(), strip.text.x =element_blank(),legend.position ="none"), 
  al_cp_50+ theme(axis.title = element_blank(),strip.text.x =element_blank(), strip.text.y.left = element_blank(),legend.position ="none"),
  wal_cp_30 + theme(legend.position ="none"),
  wal_cp_50+ theme(axis.title.y=element_blank(),  strip.text.y.left = element_blank(),legend.position ="none" ),
  nrow = 2, ncol = 2,
  common.legend = T, legend = "right",
  widths = c(1, 1), heights = c(1, 1),
  labels = c("Allenwiller", "", "Wäldi", ""),
  label.x = 0,  label.y = 1.02,
  font.label = list(size = 12, face = "bold"),
  align = "hv")
HI_diff_plot

path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/Figure_HI1_diff.png"), 
    width =7.0866, height=7.0866, units = "in", res=600)
plot(HI_diff_plot)
dev.off()

## repeat for advanced generation hybrids
wal_cp_30 <- plot_HI_diff(subset(HI_data_grid, stand == "Waldi" & Ar == 30), "Ar = 30", diff_prop_HI0.5, 0.5)
wal_cp_50 <- plot_HI_diff(subset(HI_data_grid, stand == "Waldi" & Ar == 50), "Ar = 50", diff_prop_HI0.5, 0.5)
al_cp_30 <- plot_HI_diff(subset(HI_data_grid, stand == "Allenwiller" & Ar == 30), "Ar = 30", diff_prop_HI0.5, 0.5)
al_cp_50 <- plot_HI_diff(subset(HI_data_grid, stand == "Allenwiller" & Ar == 50), "Ar = 50", diff_prop_HI0.5, 0.5)

HI_diff_plot <- ggpubr::ggarrange(
  al_cp_30+ theme(axis.title.x = element_blank(), strip.text.x =element_blank()), 
  al_cp_50+ theme(axis.title = element_blank(),strip.text.x =element_blank(), strip.text.y.left = element_blank()),
  wal_cp_30,
  wal_cp_50+ theme(axis.title.y=element_blank(),  strip.text.y.left = element_blank() ),
  nrow = 2, ncol = 2,
  widths = c(1, 1), heights = c(1, 1),
  labels = c("Allenwiller", "", "Wäldi", ""),
  label.x = 0,  label.y = 1,
  font.label = list(size = 12, face = "bold"),
  common.legend = TRUE,legend = "right",
  align = "hv" )

HI_diff_plot

path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/Figure_HI0.5_diff.png"), 
    width =7.0866, height=7.0866, units = "in", res=600)
plot(HI_diff_plot)
dev.off()

# --------- calculate and plot SSE circle plots ----------
sim_HI_cp <- readRDS("NEMO/Simulated_HI_counts_circleplots_2020.RDS")
obs_HI_cp <- fread("NEMO/Observed_HI_circleplots.csv")

obs_props <- obs_HI_cp %>%
 dplyr::select(CirclePlotID, obs_prop_HI1 = prop_HI1, obs_prop_HI0.5 = prop_HI0.5, obs_prop_HI0 = prop_HI0)

merged_prop <- sim_HI_cp %>%
  left_join(obs_props, by = "CirclePlotID") %>%
  filter(!is.na(obs_prop_HI1) & !is.na(prop_HI1))  # ensure both observed and simulated are present

# calculate SSE per circle plot accounting for the three HI categories
SSE_HI <- merged_prop %>%
  rowwise() %>%
  mutate(SSE = sum((c(prop_HI1, prop_HI0.5, prop_HI0) - c(obs_prop_HI1, obs_prop_HI0.5, obs_prop_HI0))^2) ) %>%
  ungroup()
hist(SSE_HI$SSE)

saveRDS(SSE_HI,"NEMO/NEMO_SSE_HI_circleplots_2020.RDS")

# mean SSE across replicates and per scenario
SSE_HI <- readRDS("NEMO/NEMO_SSE_HI_circleplots_2020.RDS")
SSE_summary <- SSE_HI %>%
  group_by(CirclePlotID, stand, scenario2) %>%
  mutate(mean_cp_SSE = mean(SSE, na.rm = T), # mean across replicates
         sd_cp_SSE = sd(SSE, na.rm = T)) %>% ungroup() %>%
  group_by(stand, scenario2) %>%
  mutate(mean_scenario_SSE = mean(SSE, na.rm = T),
          sd_scenario_SSE = sd(SSE, na.rm = T))


# attach values to the sampling shapefile
samplings <- c("Waldi_sampling_scheme.shp", "Allenwiller_sampling_scheme.shp")
stand_names <- c("Waldi", "Allenwiller")

SSE_data <- list()
for (i in seq_along(samplings)) {
  sampling <- read_sf(samplings[i])
  colnames(sampling)[2]<- "CirclePlotID"
  sampling <- sampling[,-c(3:13)] 
  stand_data <- SSE_summary %>% filter(stand == stand_names[i])
  SSE <- sampling %>% left_join(stand_data, by = "CirclePlotID") %>%
    filter(!is.na(scenario2))
  
  SSE_data[[i]] <- SSE
}

SSE_results_grid <- bind_rows(SSE_data)

SSE_results_grid <- SSE_results_grid %>%
  mutate(
    Ar = as.integer(str_extract(scenario2, "_a\\d+$") %>% str_remove("_a")),
    dp = case_when(
      str_detect(scenario2, "dp0") ~ "Flat\ndispersal", 
      str_detect(scenario2, "dp\\d+") ~ str_extract(scenario2, "dp\\d+") %>% str_remove("dp"),
      TRUE ~ NA_character_
    ),
    s = case_when(
      str_detect(scenario2, "s[0-9]\\.[0-9]") ~ as.numeric(str_extract(scenario2, "s[0-9]\\.[0-9]") %>% str_remove("s")),
      TRUE ~ 0
    ),
    dp = factor(dp, levels = c("Flat\ndispersal", "109", "50", "20")),
    s = factor(s, levels = sort(unique(na.omit(s))))  
  )

# read full files for background
sim_HI_files_2020 <- c("NEMO/Quanti_data_waldi_current_HI.RDS", "NEMO/Quanti_data_allen_current_HI_cut.RDS")
stand_names <- c("Waldi", "Allenwiller")
grid_files <- c("Waldi_grid_4m_5396patches.shp", "Allenwiller_grid_r4m_9844patches.shp")

grid_HI_data <- lapply(seq_along(sim_HI_files_2020), function(i) {
  file <- readRDS(sim_HI_files_2020[i])
  grid <- read_sf(grid_files[i])
  stand_data <- file %>% filter(stand == stand_names[i])
  filled_patches <- unique(stand_data$pop)
  grid$has_data <- ifelse(grid$patch.ID %in% filled_patches, TRUE, NA)
  grid$stand <- stand_names[i]
  grid
})

background_grid <- bind_rows(grid_HI_data)

## plot function
hist(SSE_results_grid$SSE)
plot_SSE <- function(data_subset, title, zoom = TRUE) {
  stand_name <- unique(data_subset$stand)
  
  data_subset <- data_subset %>%
    mutate(
      dp_label = as.character(dp),
      dp_label = ifelse(dp_label == "1000", "Flat\ndispersal", dp_label),
      dp_label = factor(dp_label, levels = c("Flat\ndispersal", "20", "50", "109"))
    )
  
  base_plot <- ggplot() +
    geom_sf(data = background_grid %>%filter(stand == stand_name, has_data == TRUE),fill = "grey90", color = NA ) +
    geom_sf(data = data_subset, aes(fill = mean_cp_SSE), color = NA) +
    #scale_fill_viridis_c(option = "viridis", name = "SSE") +
    scale_fill_gradient2(
      name =  "SSE",
      low = "#1B9E77",
      mid = "#EACB2B",
      high = "#D95F02",
      midpoint = 1,
      limits = c(0, 2),
      breaks = c( 0,0.5, 1,1.5,2 ),
      oob = scales::rescale_none
    ) +
    #scale_fill_viridis_c(option = "plasma", direction = -1, name = "SSE")+
    ggh4x::facet_grid2(dp ~ s, labeller = label_value, switch = "both", render_empty = FALSE) +
    geom_text(data = data_subset, aes( x = Inf, y = -Inf,label = paste0(round(mean_scenario_SSE, 2))),inherit.aes = FALSE, 
              size = 2, hjust =1.1, vjust = -0.5) +
    theme_minimal(base_size = 10) +
    guides(alpha = "none") +
    labs(
      title = title,
      x = "Selection strength (s)",
      y = "Dispersal distance (dp)"
    ) +
    theme_p
  
  # Zoom in for Allenwiller only
  if (stand_name == "Allenwiller"&&zoom) {
    data_bbox <- sf::st_bbox(data_subset)
    # little padding around
    padding_x <- 20
    padding_y <- 20
    
    base_plot <- base_plot +
      coord_sf(
        xlim = c(data_bbox["xmin"] - padding_x, data_bbox["xmax"] + padding_x),
        ylim = c(data_bbox["ymin"] - padding_y, data_bbox["ymax"] + padding_y),
        expand = FALSE
      )
  } else {
    # coord_sf without zooming for Waldi to avoid distortion
    base_plot <- base_plot +
      coord_sf(expand = FALSE)
  }
  
  return(base_plot)
}

wal_cp_30 <- plot_SSE(subset(SSE_results_grid, stand == "Waldi" & Ar == 30), "Ar = 30",  zoom = T)
wal_cp_50 <- plot_SSE(subset(SSE_results_grid, stand == "Waldi" & Ar == 50), "Ar = 50", zoom = T)
al_cp_30 <- plot_SSE(subset(SSE_results_grid, stand == "Allenwiller" & Ar == 30), "Ar = 30", zoom = T)
al_cp_50 <- plot_SSE(subset(SSE_results_grid, stand == "Allenwiller" & Ar == 50), "Ar = 50", zoom = T)

SSE_plot <- ggpubr::ggarrange(
  al_cp_30+ theme(axis.title.x = element_blank(), strip.text.x =element_blank()), 
  al_cp_50+ theme(axis.title = element_blank(),strip.text.x =element_blank(), strip.text.y.left = element_blank()),
  wal_cp_30,
  wal_cp_50+ theme(axis.title.y=element_blank(),  strip.text.y.left = element_blank() ),
  nrow = 2, ncol = 2,
  widths = c(1, 1), heights = c(1, 1),
  labels = c("Allenwiller", "", "Wäldi", ""),
  label.x = 0,  label.y = 1.02,
  font.label = list(size = 12, face = "bold"),
  common.legend = TRUE,legend = "right",
  align = "hv" 
)  
SSE_plot

path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/Figure_SSE_grid.png"), 
    width =7.0866, height=7.0866, units = "in", res=600)
plot(SSE_plot)
dev.off()

