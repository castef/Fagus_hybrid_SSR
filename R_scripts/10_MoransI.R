
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

# path to get files
setwd(dir="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Data_clean")

# --------- Observed Moran's I Calculation ---------
obs_HI_cp <- fread("NEMO/Observed_HI_circleplots.csv")
wal_cp_HI <- subset(obs_HI_cp, stand == "Waldi")
wal_valid_cp <- wal_cp_HI$CirclePlotID

sampling_shp <- st_read("Waldi_sampling_scheme.shp")
colnames(sampling_shp)[2]<- "CirclePlotID"
sampling_shp <- sampling_shp[sampling_shp$CirclePlotID %in% wal_valid_cp,]

# attach HI value to the shapefile
merged_obs <- sampling_shp %>%
  left_join(wal_cp_HI, by = c("CirclePlotID")) %>% filter(!is.na(prop_HI1)) 

centroids <- st_centroid(merged_obs)
coords <- st_coordinates(centroids)

# distance matrix
dist_matrix <- as.matrix(dist(coords))
inv_dist_matrix <- 1 / dist_matrix
diag(inv_dist_matrix) <- 0
inv_dist_matrix[is.infinite(inv_dist_matrix)] <- 0

# calculate Moran's I on the proportion of HI
moran_result <- Moran.I(merged_obs$prop_HI1, inv_dist_matrix)

moran_waldi <- data.frame(
  stand = "Waldi",
  morans_I = moran_result$observed,
  p_value = moran_result$p.value
)


## repeat for allenwiller
al_cp_HI <- subset(obs_HI_cp, stand == "Allenwiller")
# get the cp ID for the simulated data
al_valid_cp <- al_cp_HI$CirclePlotID

sampling_shp <- st_read("Allenwiller_sampling_scheme.shp")
colnames(sampling_shp)[2]<- "CirclePlotID"
sampling_shp <- sampling_shp[sampling_shp$CirclePlotID %in% al_valid_cp,]

# attach HI value to the shapefile
merged_obs <- sampling_shp %>%
  left_join(al_cp_HI, by = c("CirclePlotID")) %>% filter(!is.na(prop_HI1)) 

centroids <- st_centroid(merged_obs)
coords <- st_coordinates(centroids)

# distance matrix
dist_matrix <- as.matrix(dist(coords))
inv_dist_matrix <- 1 / dist_matrix
diag(inv_dist_matrix) <- 0
inv_dist_matrix[is.infinite(inv_dist_matrix)] <- 0

# calculate Moran's I
moran_result <- Moran.I(merged_obs$prop_HI1, inv_dist_matrix)

moran_allen <- data.frame(
  stand = "Allenwiller",
  morans_I = moran_result$observed,
  p_value = moran_result$p.value
)

moran_observed <- bind_rows(moran_waldi, moran_allen)
moran_observed


#-------- Simulated Moran's I Calculation ----------
sim_HI_cp <- readRDS("NEMO/Simulated_HI_counts_circleplots_2020.RDS")
stands <- c("Waldi", "Allenwiller")
samplings <- c("Waldi_sampling_scheme.shp", "Allenwiller_sampling_scheme.shp")

moran_res <- list()

#calculate moran's I for the 2 stands
for (i in seq_along(stands)) {
  stand_name <- stands[i]
  sampling_file <- samplings[i]
  
  # filter circle plots based on those where there are observed data
  valid_cp <- if (stand_name == "Waldi") wal_valid_cp else al_valid_cp
  stand_HI_cp <- sim_HI_cp %>%filter(stand == stand_name, CirclePlotID %in% valid_cp)
  
  sampling_shp <- st_read(sampling_file)
  colnames(sampling_shp)[2]<- "CirclePlotID"
  
  # filter the shapefile for the valid circle plots
  sampling_shp <- sampling_shp %>%filter(CirclePlotID %in% stand_HI_cp$CirclePlotID)
  
  # centroids and distance matrix
  centroids <- st_centroid(sampling_shp)
  coords <- st_coordinates(centroids)
  dist_matrix <- as.matrix(dist(coords))
  
  # inverse distance matrix
  inv_dist_matrix <- 1 / dist_matrix
  diag(inv_dist_matrix) <- 0
  
  #unique combinations of scenario and replicate for the current stand
  combinations <- stand_HI_cp %>% distinct(stand, scenario2, replicate)
  
  # Moran's I for each combination
  for (j in seq_len(nrow(combinations))) {
    combo <- combinations[j, ]
    combo_HI <- stand_HI_cp %>%
      filter(
        scenario2 == combo$scenario2,
        replicate == combo$replicate
      )
    

    # match CirclePlotIDs and subset inv_dist_matrix
    ids <- combo_HI$CirclePlotID
    row_inds <- match(ids, sampling_shp$CirclePlotID)
    
    # check for valid ID s
    if (any(is.na(row_inds))) {
      cat("Missing IDs in distance matrix for combination:", j, "\n")
      next
    }
    
    sub_inv_dist_matrix <- inv_dist_matrix[row_inds, row_inds]
    if (nrow(sub_inv_dist_matrix) != nrow(combo_HI)) {
      cat("Mismatch in matrix and data length for combination:", j, "\n")
      next
    }
    
    moran_result <- tryCatch({
      Moran.I(combo_HI$prop_HI1, sub_inv_dist_matrix)
    }, error = function(e) {
      combo_str <- paste(stand_name, combo$scenario2, combo$replicate, sep = "_")
      cat("Error in Moran's I calculation for combination:",j, combo_str, "\n")
      return(data.frame(observed = NA, p.value = NA))
    })
    
    # Store the results
    moran_res[[paste0("Sim_", stand_name, "_", j)]] <- data.frame(
      stand = combo$stand,
      scenario2 = combo$scenario2,
      replicate = combo$replicate,
      morans_I = moran_result$observed,
      p_value = moran_result$p.value
    )
  }
}

moran_simulated <- bind_rows(moran_res)

# check rep with NULL morans I
moran_simulated[is.na(moran_simulated$morans_I), ]

# summarise across replicates
moran_simulated <- moran_simulated %>%
  group_by(stand, scenario2) %>%
  mutate(mean_morans_I = mean(morans_I, na.rm = T), 
         sd_morans_I = sd(morans_I, na.rm = T))

#check highest and lowest values
summary(subset(moran_simulated, stand == "Waldi"))
summary(subset(moran_simulated, stand == "Allenwiller"))

# --------------------- SSE and Moran's plot -----

# combine datasets
SSE_HI <- readRDS("NEMO/NEMO_SSE_HI_circleplots_2020.RDS")

# mean SSE across replicates and per scenario
SSE_summary <- SSE_HI %>%
  group_by(stand, scenario2, replicate) %>%
  mutate(mean_rep_SSE = mean(SSE, na.rm = T), # mean across circle plots
         sd_rep_SSE = sd(SSE, na.rm = T)) %>% ungroup() %>%
  group_by(stand, scenario2) %>% # mean across replicates
  mutate(mean_scenario_SSE = mean(SSE, na.rm = T),
         sd_scenario_SSE = sd(SSE, na.rm = T))

# combine values per replicate and across replicate
combined_df <- bind_rows(
  moran_simulated %>%
    dplyr::select(stand, scenario2, replicate, value = morans_I, mean_value = mean_morans_I, sd_value = sd_morans_I) %>%
    mutate(metric = "Moran's I"),
  SSE_summary %>%
    dplyr::select(stand, scenario2, replicate, value = mean_rep_SSE , mean_value = mean_scenario_SSE, sd_value = sd_scenario_SSE) %>%
    mutate(metric = "SSE")
)

# add observed morans I
observed_moran <- moran_observed %>%
  dplyr::select(stand, morans_I) %>%
  rename(observed_value = morans_I)

# merge the observed and add the parameters (for the colors)
combined_df <- combined_df %>%
  left_join(observed_moran, by = c("stand")) %>%
  mutate(observed_value = ifelse(metric == "Moran's I", observed_value, NA),
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

# add the lable for  axis
combined_df$scenario_label <- NA
combined_df$scenario_label[grepl("^flatdispersal", combined_df$scenario2)] <- paste0("No space, Ar = ",gsub(".*_a([0-9]+)$", "\\1",combined_df$scenario2[grepl("^flatdispersal", combined_df$scenario2)]))
combined_df$scenario_label[ grepl("^spatial_", combined_df$scenario2)] <- paste0("Space (dp = ", gsub(".*_dp([0-9]+)_.*", "\\1",combined_df$scenario2[ grepl("^spatial_", combined_df$scenario2)]),"), Ar = ",  gsub(".*_a([0-9]+)$", "\\1", combined_df$scenario2[ grepl("^spatial_", combined_df$scenario2)]))
combined_df$scenario_label[grepl("^spatialselection", combined_df$scenario2)] <- paste0("Space (dp = ", gsub(".*_dp([0-9]+)_.*", "\\1",combined_df$scenario2[grepl("^spatialselection", combined_df$scenario2)]),") + Selection (s = ",  gsub(".*_s([0-9]\\.[0-9])_.*", "\\1", combined_df$scenario2[grepl("^spatialselection", combined_df$scenario2)]),"), Ar = ",gsub(".*_a([0-9]+)$", "\\1",combined_df$scenario2[grepl("^spatialselection", combined_df$scenario2)]))


# final plot
combined_plot <- ggplot(combined_df, aes(x = scenario_label, y = value)) +
  #geom_point(alpha = 0.6, size = 1) +                 # all replicate points
  geom_point(aes(y = mean_value, color = dp), size = 2) +  # mean value
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value, , color = dp), 
                 width = 0.2) +
  #scale_color_discrete(low = "gold", high = "purple")+
  geom_hline(aes(yintercept = observed_value), linetype = "dashed", color = "black") + 
  facet_grid(metric ~ stand, scales = "free_y") +     
  theme_bw() + guides(color = "none")+
  labs(
    x = "Simulation Scenario",
    y = "Value"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20, unit = "pt"),
    axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(size = 10, margin = margin(b = 5)),
    strip.text.y.left = element_text(angle = 0, size = 10, margin = margin(r = 5)),
    plot.title = element_text(size = 10, hjust = 0.5)
  )

combined_plot

path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/FigureS9_NEMO_HI_SSE.png"), width =13000,height=8000, res=1200)
plot(combined_plot)
dev.off()


