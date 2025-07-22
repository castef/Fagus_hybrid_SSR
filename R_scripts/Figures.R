
# FIGURES


# max size per figure: 10MB
# dimensions: 80 - 180 mm (7.0866 inches) width
# resolution: 300 - 600 dpi (res)
# all figures need to be named with the figure number
# use  theme_bw(base_size = 7) to standardise text size

# plot.title(size = 8)

# check which packages are not needed
install.packages("forcats") # to merge plots
install.packages("ggpubr")
install.packages("xlsx")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("grid")
install.packages("plyr")
install.packages("stringr")
install.packages("scales")
install.packages("reshape2")
install.packages("tidyr")
install.packages("ggnewscale") # to erase the fill and color scale in the same ggplot
install.packages("scatterpie") # to build pie chart
install.packages("basemaps")   # for background map
install.packages("terra")      # to import raster data
install.packages("tidyterra")  # to use raster file as background in ggplot
install.packages("cowplot")    # to create inset into the maps
install.packages("ggspatial")  # to add the scale to the maps

library(xlsx)
library(plyr)
library(dplyr)
library(ggplot2)
library(grid)
library(stringr)
library(scales)
library(reshape2)
library(tidyr)
library(ggpubr)
library(forcats)
library(ggnewscale)
library(scatterpie)
library(basemaps)
library(terra)
library(tidyterra)
library(cowplot)
library(ggspatial)
library(sf)
#library(rgeos)
library(data.table)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

path="~/Hybridization/"

setwd(dir="~/Data_clean")

# -------- Colour vectors ------
ColourVector = c( "#ffcc00","#482173FF", "#25858EFF", "#85D54AFF", "orange","#9C99C8",  "gray8")
ColourVector1 <- c("N_P1" = "#ffcc00", "N_P2" = "#482173FF", "N_F1" = "#25858EFF", "N_F2" = "#85D54AFF", "N_BC1" = "orange", "N_BC2" = "#9C99C8", "N_unassigned" = "black")
ColourVector2 = c("Pure1" = "#ffcc00", "Pure2" = "#482173FF", "F1" = "#25858EFF", "F2" = "#85D54AFF", "BC1" = "orange","BC2" ="#9C99C8", "unassigned" = "black")
ColourVector3 = c("Pure sylvatica" = "#ffcc00", "Pure hohenackeriana" = "#482173FF", "F1" = "#25858EFF", "F2" = "#85D54AFF", "BC sylvatica" = "orange","BC hohenackeriana" ="#9C99C8", "unassigned" = "gray8")

# -------- just legend ------

## with unassigned
ColourVector3 <- c("Pure sylvatica" = "#ffcc00",
                   "Pure hohenackeriana" = "#482173FF",
                   "F1" = "#25858EFF",
                   "F2" = "#85D54AFF",
                   "BC sylvatica" = "orange",
                   "BC hohenackeriana" = "#9C99C8", 
                   "Unassigned" = "black")

# dummy dataframe
dummy_df <- data.frame(group = factor(names(ColourVector3), levels = names(ColourVector3)))

# dummy plot
dummy_plot <- ggplot(dummy_df, aes(x = group, fill = group)) +
  geom_bar() +
  scale_fill_manual(values = ColourVector3) +
  guides(fill = guide_legend(ncol = 1, title = "Genotype class")) + 
  theme_bw(base_size = 7) +
  theme(legend.position = "right",
        legend.key = element_rect(fill = "black", colour = NA) ,
        strip.background = element_blank(), 
        legend.key.spacing.y = unit(0.3, 'cm'), 
        legend.text = element_text(size = 6),      # Smaller legend text
        legend.title = element_text(size = 7),     # Smaller legend title
        legend.key.size = unit(0.3, "cm"),         # Smaller legend keys
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.spacing.y = unit(0.05, "cm"),        # Less spacing
        legend.margin = margin(2, 2, 2, 2))

legend_with_unassigned <- get_legend(dummy_plot)
grid.newpage()
grid.draw(legend_with_unassigned)

# without unassigned

## with unassigned
ColourVector4 <- c("Pure sylvatica" = "#ffcc00",
                   "Pure hohenackeriana" = "#482173FF",
                   "F1" = "#25858EFF",
                   "F2" = "#85D54AFF",
                   "BC sylvatica" = "orange",
                   "BC hohenackeriana" = "#9C99C8" )

# dummy dataframe
dummy_df2 <- data.frame(group = factor(names(ColourVector4), levels = names(ColourVector4)))

dummy_plot2 <- ggplot(dummy_df2, aes(x = group, fill = group)) +
  geom_bar() +
  scale_fill_manual(values = ColourVector3) +
  guides(fill = guide_legend(ncol = 1, title = "Genotype class")) + 
  theme_bw(base_size = 7) +
  theme(legend.position = "right",
        legend.key = element_rect(fill = "black", colour = NA) ,
        strip.background = element_blank(), 
        legend.key.spacing.y = unit(0.1, 'cm'), 
        legend.text = element_text(size = 6),      # Smaller legend text
        legend.title = element_text(size = 7),     # Smaller legend title
        legend.key.size = unit(0.3, "cm"),         # Smaller legend keys
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.spacing.y = unit(0.01, "cm"),        # Less spacing
        legend.margin = margin(2, 2, 2, 2) )

legend_no_unassigned <- get_legend(dummy_plot2)
grid.newpage()
grid.draw(legend_no_unassigned)

# -------- Figure 1. Sampling scheme -------

# Waldi - load adults, trasnects and sampling circles
sampling_v <- st_read("Waldi_sampling_scheme.shp")
trans <- st_read("Waldi_transects.shp")
obs <- read.csv("Newhybrids/Waldi_indiv_info_class.csv")
obs <- obs[!is.na(obs$x),]
ndvi1=rast(paste0(path,"Data/NDVI_classification/Waldi_NDVI_50.tif"))
ndvi1=project(ndvi1, "EPSG:3035")

waldi_ss <- ggplot() +
  geom_spatraster(data = ndvi1, alpha = 0.8)+scale_fill_gradient(high="#006400",low = "white", expand = c(0,0))+
  geom_point(data = subset(obs,generation == "Adult"), aes(x = x, y = y),size = 0.1, fill ="black", color = "black", stroke = 1, shape = 21) +
  geom_sf(data = sampling_v, mapping= aes(fill= "Sampling scheme"), color = "NA", fill = "white", alpha = 0.7) +
  geom_sf(data = trans,  mapping= aes(fill= "transects"), color = "black", fill = "black", alpha = 0.5) +
  coord_sf(crs = "epsg:3035", xlim = c(min(obs$x, na.rm = T)-17, max(obs$x, na.rm = T)+28), 
           ylim = c(min(obs$y, na.rm = T)-5,max(obs$y, na.rm = T)+3))+
  theme_void() +
  theme(legend.position = "none",
    axis.title = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
    axis.text = element_blank(),
    plot.title = element_text(hjust = 0.5,vjust = 1,  size = 8, face = "bold", margin = margin(b = 0))
  ) +
  annotation_scale(location = "bl", unit_category= "metric",  width_hint = 0.5, style = "ticks", line_width = 1.2)+
  ggtitle("Wäldi (CH)")

waldi_ss          

# allenwiller
sampling_v <- st_read("Allenwiller_sampling_scheme.shp")
trans <- st_read("Allenwiller_transects.shp")
obs <- read.csv("Newhybrids/Allenwiller_indiv_info_class_subset.csv")
obs <- obs[!is.na(obs$x),]
ndvi=rast(paste0(path,"Data/NDVI_classification/All_NDVI_image.tif"))
ndvi=project(ndvi, "EPSG:3035")

# triangle
coords <- matrix(c(
  4126715, 2840045,
  4126800, 2840120,
  4126710, 2840152,
  4126715, 2840045  # Closing the polygon
), ncol = 2, byrow = TRUE)

polygon_sf <- st_sf(geometry = st_sfc(st_polygon(list(coords))),crs = 3035)

allen_ss <- ggplot() +
  geom_spatraster(data = ndvi, alpha = 0.8)+scale_fill_gradient(high="#006400",low = "white", expand = c(0,0))+
  geom_sf(data = polygon_sf, fill = NA, linetype = "dashed", linewidth = 0.5, col = "white")+
  
  geom_point(data = subset(obs,generation == "Adult"), aes(x = x, y = y), fill ="black", color = "black",size = 0.1, alpha = 0.8,stroke = 1, shape = 21) +
  geom_sf(data = sampling_v, 
          mapping= aes(fill= "Sampling scheme"), color = "NA", fill = "white", alpha = 0.7) +
  geom_sf(data = trans,  mapping= aes(fill= "transects"), color = "black", fill = "black", alpha = 0.5) +
  coord_sf(crs = "epsg:3035", xlim = c(min(obs$x, na.rm = T)-20, max(obs$x, na.rm = T)+20), 
           ylim = c(min(obs$y, na.rm = T)-20,max(obs$y, na.rm = T)+5))+
  theme_void() +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_text(hjust = 0.5,vjust = 1, size = 8, face = "bold", margin = margin(b = 0))
  ) +
  annotation_scale(location = "bl",unit_category= "metric",  width_hint = 0.5,  style = "ticks", line_width = 1.2)+
  ggtitle("Allenwiller (FR)") 

allen_ss          

          ## get the simple map
          world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf", countries = NULL)
          europe <- world[world$continent %in% c("Europe"), ]
          base_map <- ggplot() +
            geom_sf(data = europe, fill = "grey", color = "white", linewidth = 0.8) +  
            geom_point(aes(x = c(4126800, 4252950), y = c(2840149, 2724958)), size = 1.3, color = "black") + 
            coord_sf(crs = 3035, xlim = c(3600000, 4909000), ylim = c(2000000, 3500000)) +
            theme_void()
          base_map
          
          ## arrows
          base_map <- base_map +
            geom_segment(aes(x = 4126800, y = 2840149, xend = 3600000, yend = 2840149),arrow = arrow(length = unit(0.1, "inches")), color = "black") +  
            geom_segment(aes(x = 4252950, y = 2724958,xend = 4909000, yend = 2724958),arrow = arrow(length = unit(0.1, "inches")), color = "black")


### new map version in WGS 84
install.packages("basemaps")
library(basemaps)

ext <- st_bbox(c(xmin = 0, xmax = 20, ymin =40, ymax = 55), crs = st_crs(4326))

# view all available maps
get_maptypes()

# set defaults for the basemap
set_defaults(map_service = "esri", map_type = "world_physical_map")

# load and return basemap map as class of choice, e.g. as image using magick:
#basemap_magick(ext)
#basemap_ggplot(ext)
#basemap_terra(ext)

base_map <- ggplot() + 
  basemap_gglayer(ext,alpha = 0.7) + 
  theme_void()+
  scale_fill_identity() + 
  coord_sf()


base_map <- base_map +
  geom_point(aes(x = 819750.5828826678, y = 6214421.276157718),size = 2)+
  geom_point(aes(x =1012455.8008027739, y = 6046261.621746812),size = 2)+
  geom_segment(
    aes(x = 819750.5828826678, y = 6214421.276157718, xend = 29422.421664471054, yend = 6214421.276157718),
    arrow = arrow(length = unit(0.1, "inches")), linewidth = 0.7,   
    color = "black"
  ) +  
  geom_segment(
    aes(x =1012455.8008027739, y = 6046261.621746812, xend = 2081500.3829862862, yend = 6046261.621746812),
    arrow = arrow(length = unit(0.1, "inches")), linewidth = 0.7,   
    color = "black"
  )
base_map

final_plot <- (allen_ss | base_map | waldi_ss ) +
  plot_layout(guides = "collect", widths = c(0.94, 0.9,1)) &  
  theme(plot.margin = margin(5,0,5)) 

final_plot
gc()

png(filename=paste0(path, "Manuscript/Figures/Figure1_sampling_subset.png"),  
    width =7, height=3.3, units = "in", res=600)
plot(final_plot)
dev.off()

ggsave(paste0(path, "Manuscript/Figures/Figure1_sampling_subset.pdf"),  plot = final_plot, 
       width =7, height=3.3, units = "in")

# -------- Figure 2. NH results plots and maps ------
 ## Allenwiller SUBSET##

# importing dataset
path.hold = paste0(path,"NewHybrids/Results/Allenwiller/")
# importing output newhybrid and individuals ID and merging
nh_output = read.table(paste0(path.hold, "/Allenwiller_NH_all.txt_PofZ.txt"), header = TRUE)[,-1]
Indiv_name = read.table(paste0(path.hold, "/Allenwiller_all_individuals.txt"))
nh_output$IndivName = as.factor(Indiv_name$V1)
colnames(nh_output) = c("IndivName", "Pure2", "Pure1","F1", "F2", "BC2", "BC1") # assigning names with Pure1 = sylvatica, Pure2 = hohenackeriana (reversed!!!)
subset <- fread("NewHybrids/Allenwiller_indiv_info_class_subset.csv")
nh_output <- nh_output[nh_output$IndivName %in% subset$SampleID,]
  
# grouping individuals according to the genotype class with highest posterior probability
nh_output$class_max = colnames(nh_output[,2:7])[apply(nh_output[,2:7],1,which.max)]
nh_output$class_max_value = apply((nh_output[2:7]), 1, max)

# reordering individuals according to the class with max value and that value (decreasing order)
nh_output$class_max <- as.character(nh_output$class_max)
nh_output$class_max <- factor(nh_output$class_max, ordered=T, levels=c("Pure1","F1","F2", "BC1","Pure2"))
nh_output_reord = nh_output[order(nh_output$class_max, -nh_output$class_max_value),]

# data stacked by individual and reordering
NH_melt = reshape2::melt(data = nh_output_reord[, 1:7], id.vars = "IndivName") 
colnames(NH_melt) = c("IndivName", "PopProb", "CumProb")
NH_melt$IndivName= factor(NH_melt$IndivName, levels=as.character(nh_output_reord$IndivName)) 

# add information of ageclass 
age_class = read.csv("NewHybrids/Allenwiller_indiv_info_class.csv")
colnames(age_class)[2] <- "IndivName"
NH_melt = merge(NH_melt, age_class[, c(2,11)], by = "IndivName")

# renaming for the plot
NH_melt$PopProb2 = as.character(NH_melt$PopProb)
NH_melt$PopProb2[NH_melt$PopProb=="Pure1"] = "Pure sylvatica"
NH_melt$PopProb2[NH_melt$PopProb=="Pure2"] = "Pure hohenackeriana"
NH_melt$PopProb2[NH_melt$PopProb=="F1"] = "F1"
NH_melt$PopProb2[NH_melt$PopProb=="F2"] = "F2"
NH_melt$PopProb2[NH_melt$PopProb=="BC1"] = "Backcross sylvatica"
NH_melt$PopProb2[NH_melt$PopProb=="BC2"] = "Backcross hohenackeriana"
NH_melt$PopProb2 = factor(NH_melt$PopProb2,levels = c("Pure sylvatica","Pure hohenackeriana","F1","F2","Backcross sylvatica", "Backcross hohenackeriana"))

al_nh_plot_subset = ggplot(NH_melt, aes(x = IndivName, y=as.numeric(CumProb), fill=PopProb2))+ # CumProb changed to numeric
  geom_bar(stat="identity", position = "stack", width = 1) + 
  scale_fill_manual(values=ColourVector)+
  scale_colour_manual(values=ColourVector, drop = FALSE)+ylab("Cumulative Probability")+xlab("Individual") +
  scale_y_continuous(limits = c(0, 1.05), expand=c(0, 0)) + 
  theme(panel.background = element_rect(fill = "white", colour = "white"), 
        panel.spacing = unit(0.1, "cm"),
        plot.margin = margin(t = 8, r = 0, b = 3, l = 1),
        plot.background = element_rect(colour = "white"),
        plot.title = element_text(hjust = 0.5, vjust = 1, size = 8, face = "bold"),
        strip.background = element_blank(), 
        strip.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 7, vjust = -3), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7), 
        axis.title.y = element_text(size = 7),
        ) + 
  ggtitle("Allenwiller")+
  guides(colour = "none",fill="none") +
  ggh4x::facet_wrap2(~Ageclass_name, drop = T, scales = "free_x", strip.position = "right", ncol = 1,
                     labeller = labeller(Ageclass_name = c("Seedling" = "Seedlings", "Sapling" = "Saplings", "Juvenile" = "Juveniles", "Adult" = "Adults")))
al_nh_plot_subset

## Waldi ##

# importing dataset
path.hold = paste0(path,"NewHybrids/Results/Waldi/")
# importing output newhybrid and individuals ID and merging
nh_output1 = read.table(paste0(path.hold, "/Waldi_NH_all.txt_PofZ.txt"), header = TRUE)[,-1]
Indiv_name1 = read.table(paste0(path.hold, "/Waldi_NH_all_individuals.txt"))
nh_output1$IndivName = as.factor(Indiv_name1$V1) # make the individuals factors to block out the data
colnames(nh_output1) = c("IndivName", "Pure1", "Pure2","F1", "F2", "BC1", "BC2") # assigning names with Pure1 = sylvatica, Pure2 = hohenackeriana (reversed!!!)
nh_output1 = nh_output1[order(nh_output1$IndivName),] 

# data stacked by individual and reordering
NH_melt1 = reshape2::melt(data = nh_output1, id.vars = "IndivName") 
colnames(NH_melt1) = c("IndivName", "PopProb", "CumProb")
NH_melt1$PopProb = factor(NH_melt1$PopProb, levels = c("Pure1", "Pure2", "F1", "F2", "BC1", "BC2")) 

# grouping individuals according to the genotype class with highest posterior probability
nh_output1$class_max = colnames(nh_output1[,2:7])[apply(nh_output1[,2:7],1,which.max)]
nh_output1$class_max_value = apply((nh_output1[2:7]), 1, max)

# reordering individuals according to the class with max value and that value (decreasing order)
nh_output1$class_max <- as.character(nh_output1$class_max)
nh_output1$class_max <- factor(nh_output1$class_max, ordered=T, levels=c("Pure1","F1","F2", "BC1","Pure2"))
nh_output_reord1 = nh_output1[order(nh_output1$class_max, -nh_output1$class_max_value),]

# data stacked by individual and reordering
NH_melt1 = reshape2::melt(data = nh_output_reord1[, 1:7], id.vars = "IndivName") 
colnames(NH_melt1) = c("IndivName", "PopProb", "CumProb")
NH_melt1$IndivName= factor(NH_melt1$IndivName, levels=as.character(nh_output_reord1$IndivName)) 
# add information of generation 
age_class1 = read.csv("NewHybrids/Waldi_indiv_info_class.csv")
colnames(age_class1)[2] <- "IndivName"
NH_melt1 = merge(NH_melt1, age_class1[, c(2,11)], by = "IndivName")

# renaming for the plot
NH_melt1$PopProb2 = as.character(NH_melt1$PopProb)
NH_melt1$PopProb2[NH_melt1$PopProb=="Pure1"] = "Pure sylvatica"
NH_melt1$PopProb2[NH_melt1$PopProb=="Pure2"] = "Pure hohenackeriana"
NH_melt1$PopProb2[NH_melt1$PopProb=="F1"] = "F1"
NH_melt1$PopProb2[NH_melt1$PopProb=="F2"] = "F2"
NH_melt1$PopProb2[NH_melt1$PopProb=="BC1"] = "Backcross sylvatica"
NH_melt1$PopProb2[NH_melt1$PopProb=="BC2"] = "Backcross hohenackeriana"
NH_melt1$PopProb2 = factor(NH_melt1$PopProb2,levels = c("Pure sylvatica","Pure hohenackeriana","F1","F2","Backcross sylvatica", "Backcross hohenackeriana"))

# plot
wal_nh_plot = ggplot(NH_melt1, aes(x = IndivName, y=as.numeric(CumProb), fill=PopProb2))+ # CumProb changed to numeric
  geom_bar(stat="identity", position = "stack", width = 1) + 
  scale_fill_manual(values=ColourVector)+
  scale_colour_manual(values=ColourVector)+ylab("Cumulative Probability")+xlab("Individual") +
  scale_y_continuous(limits = c(0, 1.05), expand=c(0, 0)) + 
  theme(panel.background = element_rect(fill = "white", colour = "white"), 
        panel.spacing = unit(0.1, "cm"),
        plot.margin = margin(t = 8, r = 0, b = 3, l = 1),
        plot.background = element_rect(colour = "white"),
        plot.title = element_text(hjust = 0.5, vjust = 1, size = 8, face = "bold"),
        strip.background = element_blank(), 
        strip.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 7, vjust = -3), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7), 
        axis.title.y = element_text(size = 7),
        )+ 
  ggtitle("Wäldi")+
  guides(colour = "none",fill="none") +
  ggh4x::facet_wrap2(~Ageclass_name, drop = T, scales = "free_x", strip.position = "right", ncol = 1,
                     labeller = labeller(Ageclass_name = c("Seedling" = "Seedlings", "Sapling" = "Saplings", "Juvenile" = "Juveniles", "Adult" = "Adults")))

wal_nh_plot

# ------ maps of the genotype distribution 


# importing circle plot data, tree info, ndvi 

## Allenwiller SUBSET ## 
all_tree = fread(paste0(path, "Data_clean/NewHybrids/Allenwiller_indiv_info_class_subset.csv"))
all_cp_info = fread(paste0(path, "Data_clean/Allenwiller_circleplot_info_withcoord.csv"))
ndvi=rast(paste0(path,"Data/NDVI_classification/All_NDVI_image.tif"))
ndvi=project(ndvi, "EPSG:3035")
# remove unassigned
all_tree <- all_tree[all_tree$NH_class_t80 != "unassigned"]
# attach mother tree coords to cp
all_cp = merge(all_tree, all_cp_info[, c(2,16,17)], by = "CirclePlotID", all.x = T)

mt_cp <- all_cp %>%
  filter(generation == "Offspring", !is.na(x_mothertree), !is.na(y_mothertree)) %>%
  group_by(x = x_mothertree, y = y_mothertree) %>%
  summarise(
    N_tot = n(),
    N_P1 = sum(NH_class_t80_2 == "Pure sylvatica", na.rm = TRUE),
    N_P2 = sum(NH_class_t80_2 == "Pure orientalis", na.rm = TRUE),
    N_F1 = sum(NH_class_t80_2 == "F1", na.rm = TRUE), 
    N_BC1 = sum(NH_class_t80_2 == "BC sylvatica", na.rm = TRUE),
    N_BC2 = sum(NH_class_t80_2 == "BC orientalis", na.rm = TRUE),
    N_F2 = sum(NH_class_t80_2 == "F2", na.rm = TRUE),
    #N_unassigned = sum(NH_class_t80 == "unassigned", na.rm = TRUE),
    .groups = "drop"
  )

# count offpsirng genotype but for the plots wihtout mother tree
no_mt_cp <- all_cp%>%
  filter(generation == "Offspring", is.na(x_mothertree) | is.na(y_mothertree)) %>%
  group_by(x, y) %>%
  summarise(
    N_tot = n(),
    N_P1 = sum(NH_class_t80_2 == "Pure sylvatica", na.rm = TRUE),
    N_P2 = sum(NH_class_t80_2 == "Pure orientalis", na.rm = TRUE),
    N_F1 = sum(NH_class_t80_2 == "F1", na.rm = TRUE), 
    N_BC1 = sum(NH_class_t80_2 == "BC sylvatica", na.rm = TRUE),
    N_BC2 = sum(NH_class_t80_2 == "BC orientalis", na.rm = TRUE),
    N_F2 = sum(NH_class_t80_2 == "F2", na.rm = TRUE),
    #N_unassigned = sum(NH_class_t80 == "unassigned", na.rm = TRUE),
    .groups = "drop"
  )
all_cp_grouped <- bind_rows(mt_cp, no_mt_cp)

## find the pie charts that are overlaping and move slightly one coords to improve visibility
all_cp_grouped[which.min(all_cp_grouped$y), "x"] <- all_cp_grouped[which.min(all_cp_grouped$y), "x"] - 15
all_cp_grouped[which(all_cp_grouped$y == 2840074.77452002 ), "y"] <- all_cp_grouped[which(all_cp_grouped$y == 2840074.77452002), "y"] + 10

## convert orientalis --> hohenackeriana
all_tree$NH_class_t80_2 <- gsub("orientalis", "hohenackeriana", all_tree$NH_class_t80_2)


## Waldi ##
# importing tree data info
wal_tree = fread(paste0(path, "Data_clean/NewHybrids/Waldi_indiv_info_class.csv"))
wal_cp_info = fread(paste0(path, "Data_clean/Waldi_circleplot_info_withcoord.csv"))
ndvi1=rast(paste0(path,"Data/NDVI_classification/Waldi_NDVI_50.tif"))
ndvi1=project(ndvi1, "EPSG:3035")

# remove unassigned
wal_tree <- wal_tree[wal_tree$NH_class_t80_corrected != "unassigned"]

# attach mother tree coords to cp
wal_cp = merge(wal_tree, wal_cp_info[, c(2,16,17)], by = "CirclePlotID")

# count and group offspring data for cp around mother tree
mt_cp <- wal_cp %>%
  filter(generation != "Adult", !is.na(x_mothertree), !is.na(y_mothertree)) %>%
  group_by(x = x_mothertree, y = y_mothertree) %>%
  summarise(
    N_tot = n(),
    N_P1 = sum(NH_class_t80_corrected == "Pure sylvatica", na.rm = TRUE),
    N_P2 = sum(NH_class_t80_corrected == "Pure orientalis", na.rm = TRUE),
    N_F1 = sum(NH_class_t80_corrected == "F1", na.rm = TRUE), 
    N_BC1 = sum(NH_class_t80_corrected == "BC sylvatica", na.rm = TRUE),
    N_BC2 = sum(NH_class_t80_corrected == "BC orientalis", na.rm = TRUE),
    N_F2 = sum(NH_class_t80_corrected == "F2", na.rm = TRUE),
    #N_unassigned = sum(NH_class_t80_corrected == "unassigned", na.rm = TRUE),
    .groups = "drop"
  )

# count offpsirng genotype but for the BA plots (they have no mother tree!!)
B_R_cp <- wal_cp %>%
  filter(generation != "Adult", is.na(x_mothertree) | is.na(y_mothertree)) %>%
  group_by(x, y) %>%
  summarise(
    N_tot = n(),
    N_P1 = sum(NH_class_t80_corrected == "Pure sylvatica", na.rm = TRUE),
    N_P2 = sum(NH_class_t80_corrected == "Pure orientalis", na.rm = TRUE),
    N_F1 = sum(NH_class_t80_corrected == "F1", na.rm = TRUE), 
    N_BC1 = sum(NH_class_t80_corrected == "BC sylvatica", na.rm = TRUE),
    N_BC2 = sum(NH_class_t80_corrected == "BC orientalis", na.rm = TRUE),
    N_F2 = sum(NH_class_t80_corrected == "F2", na.rm = TRUE),
    #N_unassigned = sum(NH_class_t80_corrected == "unassigned", na.rm = TRUE),
    .groups = "drop"
  )

# merge
wal_cp_grouped <- bind_rows(mt_cp, B_R_cp)
wal_tree$NH_class_t80_corrected <- gsub("orientalis", "hohenackeriana", wal_tree$NH_class_t80_corrected)


### maps

# allenwiller subset (different coordinates)
all_map_subset= ggplot()+
  geom_spatraster(data = ndvi)+
  scale_fill_gradient(high="gray50",low = "white", expand = c(0,0))+
  labs(fill = "Ndvi")+
  guides(fill = "none")+
  ggnewscale::new_scale("fill")+
  new_scale("fill")+
  geom_scatterpie(data = all_cp_grouped, mapping = aes(x = x, y = y, r=8),colour="black",alpha = 0.8,cols = c("N_P1", "N_P2", "N_F1", "N_F2" ,"N_BC1" ,"N_BC2"))+
  scale_fill_manual(values = ColourVector1,labels=c("Pure sylvatica", "Pure hohenackeriana", "F1","F2", "BC slvatica", "unassigned"))+
  geom_point(data=subset(all_tree,generation == "Adult"), mapping = aes(x=x,y=y,color=NH_class_t80_2),size=0.8, alpha = 0.8)+
  scale_color_manual(values = ColourVector3, expand = c(0, 0))+
  theme_void()+ guides(fill = "none")+
  coord_sf(
    xlim = c(min(all_tree$x, na.rm = T)-10, max(all_tree$x, na.rm = T)+10),
    ylim = c(min(all_tree$y, na.rm = T)-20, max(all_tree$y, na.rm = T)+6)
  )+
  theme( plot.margin = margin(t = 2, r = 1, b = 0, l = 1)) + 
  guides(colour = "none",fill="none") +
  annotation_scale(location = "bl", width_hint = 0.5, style = "ticks", line_col = "black", text_col = "black",line_width = 1.2, text_cex = 0.5)
all_map_subset

## waldi map
wal_map = ggplot()+
  geom_spatraster(data = ndvi1)+
  scale_fill_gradient(high="gray50",low = "white", expand = c(0,0))+
  guides(colour = "none", fill = "none")+
  new_scale("fill")+ 
  geom_scatterpie(data = wal_cp_grouped,mapping = aes(x = x, y = y, r=8),colour="black",alpha = 0.8,cols = c("N_P1", "N_P2", "N_F1"))+
  scale_fill_manual(values = c("N_P1" = "#ffcc00", "N_P2" = "#482173FF", "N_F1" = "#25858EFF"), labels=c("Pure sylvatica", "Pure hohenackeriana", "F1","unassigned"))+
  geom_point(data=subset(wal_tree,generation == "Adult"),aes(x=x,y=y,color=NH_class_t80_corrected),size=0.8, alpha = 0.8)+
  scale_color_manual(values = ColourVector3,labels=c("Pure sylvatica", "Pure hohenackeriana", "F1"),expand = c(0, 0))+
  theme_void()+ guides(colour = "none", fill = "none")+
  theme( plot.margin =  margin(t = 2, r = 1, b = 0, l = 1))+
  coord_sf(expand = F, 
           xlim= c(min(wal_tree$x, na.rm = T) - 10, max(wal_tree$x, na.rm = T) + 15), 
           ylim= c(min(wal_tree$y, na.rm = T) - 10, max(wal_tree$y, na.rm = T) + 15)) +
  annotation_scale(location = "bl", width_hint = 0.5, style = "ticks", line_col = "black", 
                   text_col = "black", line_width = 1.2, text_cex = 0.5)

wal_map
gc()

# subset version
Fig_NH_subset <- ggpubr::ggarrange(
  al_nh_plot_subset, wal_nh_plot, all_map_subset, wal_map,
  nrow = 2, ncol = 2,
  labels = c("(a)","","(b)",""), font.label = list(size = 8) ,
  align = "hv",
  widths = c(3,3,1,3,3,1),
  heights = c(1.4, 2)  #
)  

Fig_NH_subset <- ggarrange(Fig_NH_subset, legend_no_unassigned, ncol = 2, widths = c(5,1))
Fig_NH_subset


png(filename=paste0(path, "Manuscript/Figures/Figure2_NH.png"), 
    width =7, height=6.5, units = "in", res=600)
plot(Fig_NH_subset)
dev.off()

ggsave(paste0(path, "Manuscript/Figures/Figure2_NH.pdf"),  plot = Fig_NH_subset, 
       width =7, height=6.5, units = "in")

# -------- Figure 3. phenology ------

# custom colors for subspecies
mycols <- c("orientalis" = "#482173FF", "sylvatica" = "#ffcc00", "hybrid" = "#25858EFF")
colSet <- scale_colour_manual(values = mycols, labels = c("sylvatica"="European beech", "orientalis"="Caucasian beech"))
fillSet <- scale_fill_manual(values = mycols)

# GDD and CD annual trends per site
t_dat = fread("Phenology/T_GDD_CD_all_sites.csv")

GDD_CD_sites = ggplot(subset(t_dat,site =="Waldi"), aes(x=date2 )) +
  geom_area( aes(y=CD*10),position = "identity",linewidth = 0.9, alpha = 0.3) + # multiply by 10 to get the same range than GDD
  scale_y_continuous(name = "Growing degree days", sec.axis = sec_axis(~./10, name="Accumulated chilling days"))+
  geom_line( aes(y=Cum_GDD), linewidth = 0.5) + theme_bw()+
  labs(x = "Year")+
  ggtitle("    ")+  theme_bw(base_size = 7) +
  scale_color_gradient(low = "yellow", high = "red")+
  facet_wrap(~site, ncol = 1)+
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom", 
        strip.text = element_blank(),
        legend.key = element_rect(colour = "white"),
        legend.key.spacing.y = unit(0.4, 'cm'), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        plot.title = element_text(hjust = 0.5))
GDD_CD_sites

# plot for year variation
all_dat <- fread("Phenology/Spring_phenology_all_sites_final.csv")
all_dat <- all_dat[site != "Leiselheim"]
all_dat_long <-  all_dat %>% pivot_longer(cols = c("Cum_GDD","CD", "doy"),names_to = "Measure", values_to = "value" )

# wilcox test for each observation and 
wilcox_df <- all_dat_long %>%
  filter(site == "Waldi", species %in% c("orientalis", "sylvatica")) %>%  
  group_by(year, Measure, value) %>%
  #filter(n_distinct(species) >= 2) %>%  
  summarise(
    p = tryCatch(wilcox.test(average.stage ~ species)$p.value, error = function(e) NA),
    .groups = "drop"
  ) %>%
  mutate(
    p_signif = cut(p,
                   breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, 1),
                   labels = c("***", "**", "*", ".", "ns")),
    label_y = 4.15  
  )

# mean and SE for the various measures
doy_summary <- all_dat %>%
  filter(site == "Waldi") %>%
  group_by(year, species, doy) %>%
  summarise(
    x_var = unique(doy),
    mean_stage = mean(average.stage, na.rm = TRUE),
    se_stage = sd(average.stage, na.rm = TRUE) / sqrt(n()),
    Measure = "doy",
    .groups = "drop"
  )

gdd_summary <- all_dat %>%
  filter(site == "Waldi") %>%
  group_by(year, species, Cum_GDD) %>%
  summarise(
    x_var = unique(Cum_GDD),
    mean_stage = mean(average.stage, na.rm = TRUE),
    se_stage = sd(average.stage, na.rm = TRUE) / sqrt(n()),
    Measure = "Cum_GDD",
    .groups = "drop"
  )

cd_summary <- all_dat %>%
  filter(site == "Waldi") %>%
  group_by(year, species,CD) %>%
  summarise(
    x_var = unique(CD),
    mean_stage = mean(average.stage, na.rm = TRUE),
    se_stage = sd(average.stage, na.rm = TRUE) / sqrt(n()),
    Measure = "CD",
    .groups = "drop"
  )

summary<- bind_rows(doy_summary, gdd_summary, cd_summary)
# remove hybrids 
summary <- subset(summary, species != "hybrid")

# plot only gdd and doy
waldi_pheno <- ggplot(subset(summary,Measure !="CD"),aes(x = x_var, y = mean_stage, col = species, fill = species)) +
  geom_errorbar(aes(ymin = mean_stage - se_stage, ymax = mean_stage + se_stage),
                width = 0.3, size = 0.5) +
  geom_point(shape = 21, size = 0.7) +
  geom_line(aes(group = species), size = 0.5) +
  geom_text(data = subset(wilcox_df,Measure!="CD"),
            aes(x = value, y = label_y, label = p_signif),
            inherit.aes = FALSE,
            size = 2) +
  ggh4x::facet_grid2(year ~ Measure, scales = "free_x",
             labeller = labeller(Measure = c(doy = "Day of the year (DOY)", Cum_GDD = "Growing degree days (GDD)")), switch = "x") +
  colSet + fillSet + theme_bw(base_size = 7) +labs(color = "Species")+guides(fill = "none")+
  ylab("Average phenological stage") + xlab("") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom", 
        strip.text = element_text(size = 7, angle= 0),
        legend.key = element_rect(colour = "white"),
        legend.key.spacing.y = unit(0.4, 'cm'), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        plot.title = element_text(hjust = 0.5), 
        strip.placement = "outside",
        strip.switch.pad.grid = unit(0.1, "cm")) +
  ggtitle("Wäldi")
waldi_pheno

# single boxplot allenwiller
all_boxplot <- ggplot(data = subset(all_dat,site == "Allenwiller"), aes(x = as.factor(doy), y = average.stage, fill = species, group = species)) +
  #geom_jitter( aes(color = species), size = 0.8)+
  #geom_point(position=position_jitterdodge( jitter.width = 0.1, jitter.height = 0), shape = 21, alpha = 0.9, size = 2) + 
  #stat_smooth(aes(group=species, fill = species, col = species))+ 
  geom_boxplot(alpha = 0.8, col = "black", linewidth =0.3, outlier.size = 0.7, width = 0.5) +
  fillSet + colSet + guides(color = "none")+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  #geom_text(data = subset(avg_labels, site == "Allenwiller"), aes(x = Inf, y = -Inf, label = label), inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 3.5)
  stat_compare_means(
    method = "wilcox.test",
    aes(group = species),
    label = "p.signif",
    label.y = wilcox_df$label_y - 0.1,     # same as other plot
    size = 2.5
  )+
  theme_bw(base_size = 7) +
  ylab("Average phenological stage") +labs(color = "Species", fill = "Species")+
  xlab("Day of the year (DOY)\n     \n    ") +
  facet_grid( ~year) +
  ggtitle("Allenwiller")+
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom", 
        strip.text =element_blank(),
        legend.key = element_rect(colour = "white"),
        legend.key.spacing.y = unit(0.4, 'cm'), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        plot.title = element_text(hjust = 0.5))
all_boxplot

Fig_pheno <- ggpubr::ggarrange(
  waldi_pheno , 
  ggarrange(GDD_CD_sites,all_boxplot + theme(legend.position = "none"),
           nrow = 2, labels = c("(b)", "(c)"),  font.label=list(size=8),align = "hv"), ncol = 2, labels = c("(a)", ""), font.label=list(size=8), 
  widths = c(3,2),  align = "hv", common.legend = T, legend = "bottom")  

Fig_pheno

path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/Figure3_phenology.png"), 
    width =6.2, height=4.8, units = "in", res=600)
plot(Fig_pheno)
dev.off()

ggsave(paste0(path, "Manuscript/Figures/Figure3_phenology.pdf"),  plot = Fig_pheno, 
       width =6.2, height=4.8,units = "in")

# -------- Figure 4 and S4. NMpi2 selection gradients + genealogies ---------

### panel  selection gradients
path = "C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Parentage analysis/NMpi2/Camilla_data/"
sel_grad = read.xlsx(paste0(path, "/NMpi2_overall_results.xlsx"), sheetName = "Selection gradients")

sel_grad$Site <- gsub("Waldi", "Wäldi", sel_grad$Site)
sel_grad$Fem_male <- gsub("Female", "Mother tree (Seed donor)", sel_grad$Fem_male)
sel_grad$Fem_male <- gsub("Male", "Father tree (Pollen donor)", sel_grad$Fem_male)

## effect of selection gradients on reproductinve success
sel_grad$labels2 = factor(sel_grad$Sel_gradient,
                          levels = c("DBH", "Competition_coeff", "Psyl", "Psyl biparental"),
                          labels = c("Tree size", "Competition", "P(syl)", "Assortative\nmating"))

sel_grad_clean <- subset(sel_grad,  !(labels2 == "Assortative\nmating" & Fem_male == "Mother tree (Seed donor)"))
sel_grad_clean$labels2 <- factor(sel_grad_clean$labels2,levels = c("Tree size", "Competition", "P(syl)", "Assortative\nmating"))

sel_grad_plot = ggplot(sel_grad_clean,mapping = aes(x = as.numeric(Effect), y = labels2, color=Fem_male)) + 
  geom_point(position = position_dodge(width = 0.5), size = 1.3) +
  geom_errorbarh(aes(xmin = as.numeric(Effect) - SE*1.96, xmax = as.numeric(Effect) + SE*1.96), 
                 height = 0.3, size = 0.6, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("black", "grey50"))+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.6) +
  facet_grid( ~ Site) +
  labs(y = "Phenotypic character", x = "Standardized effect on reproductive success", color = "Parent type") +
  xlim(c(-1.9, 1.9)) +
  guides(fill = "none")+
  theme_bw(base_size = 7) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0.5, "cm"),
          strip.background = element_blank(),
          legend.key = element_rect(colour = "white"),
          legend.key.spacing.y = unit(0.4, 'cm'), 
          strip.text.x = element_text(hjust = 0.5,vjust = 1, size = 8, face = "bold", margin = margin(b = 0)),
          axis.line.x.bottom = element_line(color = "black")
  )
sel_grad_plot

# genealogies -- data prep
par_results = read.csv("Mating_model/Parentage_results.csv")
par_results <- par_results %>%mutate(across(where(is.character), ~ gsub("orientalis", "hohenackeriana", .)))
par_results <- par_results %>%mutate(across(where(is.character), ~ gsub("Waldi", "Wäldi", .)))

filtered <- par_results %>% filter(Pr1 > 0.8)
filtered <- par_results %>% filter(Prog_NH != "unassigned")

# classify offspring into no parents, 2 parents and both parents assigned
both_parents <- filtered %>%filter(!is.na(Mo1_ID) & !is.na(Fa1_ID))
one_parent <- filtered %>% filter(xor(is.na(Mo1_ID), is.na(Fa1_ID)))
no_parents <- filtered %>%filter(is.na(Mo1_ID) & is.na(Fa1_ID))

combined <- bind_rows(
  both_parents %>% mutate(group = "Both parents"),
  one_parent %>% mutate(group = "One parent"),
  no_parents %>% mutate(group = "No parents")
)
combined <- combined[!(combined$Prog_NH == "unassigned"),]
combined$group <- factor(combined$group, levels = c("Both parents", "One parent", "No parents"))
combined$Prog_NH <- factor(combined$Prog_NH,levels = c("Pure sylvatica","Pure hohenackeriana","F1","F2","BC sylvatica", "BC hohenackeriana"))
combined$Prog_NH <- factor(combined$Prog_NH,levels = names(ColourVector3))

  # Figure S4. NMpi2 genealogies : genotpye of single parent of each offspring
one_parent$Mo1_NH <- factor(one_parent$Mo1_NH, levels = c("Pure sylvatica", "Pure hohenackeriana", "F1", "F2", "BC sylvatica", "BC hohenackeriana", "unassigned"))


  single_parents <- ggplot(one_parent, aes(x = Prog_NH, fill = Mo1_NH)) +
    geom_bar(position = "fill", col = "black", width = 0.3, linewidth = 0.3,) +
    scale_fill_manual(values = ColourVector3, 
                      labels = c("Pure\nsylvatica", "Pure\nhohenack.", "F1",  "Unassigned") ) +
    scale_x_discrete(labels = c("Pure\nsylvatica", "Pure\nhohenack.", "F1", "F2","BC sylvatica", "BC\nhohenackeriana", "Unassigned"))+
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "Offspring genotype", y = "Percentage", fill = "Parent NH Class") +
    theme_bw(base_size = 7)+ 
    ggh4x::facet_wrap2(~Test, labeller = labeller(Test = c("Allenwiller"="Allenwiller","Waldi"="Wäldi")))+
    theme(panel.background = element_rect(fill = "white", colour = "black"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0.5, "cm"),
          strip.background = element_blank(),
          legend.key = element_rect(colour = "white"),
          legend.key.spacing.y = unit(0.4, 'cm'), 
          #strip.text.x = element_blank(),
          axis.line.x.bottom = element_line(color = "black")
    )
  single_parents

path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/FigureS4_NMpi2_single_parents.png"), 
    width =7,height=4, res=600, units = "in")
plot(single_parents)
dev.off()

ggsave(paste0(path, "Manuscript/Figures/FigureS4_NMpi2_single_parents.pdf"),  plot = single_parents, 
       width =7,height=4,units = "in")


## Figure 4: offspring genotype x genealogy 
offspring_genealogy <- ggplot(combined, aes(x = Prog_NH, fill = group)) +
  geom_bar(position = "fill", col = "black", width = 0.5, linewidth = 0.3) +
  scale_fill_manual(values = c("white", "grey", "black")) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_discrete(labels = c("Pure\nsylvatica", "Pure\nhohenack.", "F1", "F2", "BC sylvatica" ))+
  labs(x = "Offspring genotype", y = "Percentage", fill = "N. of assigned\nparents") +
  theme_bw(base_size = 7) + facet_wrap(~Test)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.5, "cm"),
        strip.background = element_blank(),
        legend.key = element_rect(colour = "white"),
        legend.key.spacing.y = unit(0.4, 'cm'), 
        strip.text.x = element_blank(),
        axis.line.x.bottom = element_line(color = "black")
  )
offspring_genealogy

#  spatial relationship of F1 with their parents 
# euclidean distance between offspring and hohenackeriana and sylavtica parent
f1_dist <- combined %>%
  filter(Prog_NH == "F1") %>%
  mutate(
    hohenackeriana_X = ifelse(Mo1_NH == "Pure hohenackeriana", Mo1_X,ifelse(Fa1_NH == "Pure hohenackeriana", Fa1_X, NA)),
    hohenackeriana_Y = ifelse(Mo1_NH == "Pure hohenackeriana", Mo1_Y,ifelse(Fa1_NH == "Pure hohenackeriana", Fa1_Y, NA)),
    
    sylvatica_X = ifelse(Mo1_NH == "Pure sylvatica", Mo1_X, ifelse(Fa1_NH == "Pure sylvatica", Fa1_X, NA)),
    sylvatica_Y = ifelse(Mo1_NH == "Pure sylvatica", Mo1_Y, ifelse(Fa1_NH == "Pure sylvatica", Fa1_Y, NA)),
    
    dist_to_hohenackeriana = sqrt((Prog_X - hohenackeriana_X)^2 + (Prog_Y - hohenackeriana_Y)^2),
    dist_to_sylvatica  = sqrt((Prog_X - sylvatica_X)^2 + (Prog_Y - sylvatica_Y)^2)
  )

dist_compare <- f1_dist %>%
  dplyr::select(Test, dist_to_hohenackeriana, dist_to_sylvatica) %>%
  pivot_longer(cols = starts_with("dist_to_"), names_to = "Parent_NH", values_to = "Distance") %>%
  mutate(Parent_Type = recode(Parent_NH, dist_to_hohenackeriana = "Pure hohenackeriana",dist_to_sylvatica = "Pure sylvatica"))

dist_compare$Parent_Type <- factor(dist_compare$Parent_Type, levels = c("Pure sylvatica", "Pure hohenackeriana" ))

## check mean and max distnace of the two parents genotypes: 
dist_compare %>%
  group_by(Test, Parent_NH) %>%
  summarise(
    mean_distance = mean(Distance, na.rm = TRUE),
    max_distance = max(Distance, na.rm = TRUE),
    .groups = "drop"
  )

ColourVector3 = c("Pure sylvatica" = "#ffcc00", "Pure hohenackeriana" = "#482173FF", "F1" = "#25858EFF", "F2" = "#85D54AFF", "BC sylvatica" = "orange","BC hohenackeriana" ="#9C99C8", "unassigned" = "gray8")

# remove outlier for beeter visualization
dist_compare <- dist_compare[which(dist_compare$Distance <=200),]

distance_plot <- ggplot(dist_compare, aes(x = Parent_Type, y = Distance, fill = Parent_Type, colour = Parent_Type)) +
  geom_boxplot(alpha = 0.8, col = "black",width = 0.3, linewidth = 0.3,outlier.size = 0.5) +
  scale_fill_manual(values = ColourVector3)+  
  scale_color_manual(values = ColourVector3)+  
  facet_wrap(~Test) +
  labs(fill = "Parent genotype", x = "Parent genotype", y = "Distance (m) from F1 offspring") +
  theme_bw(base_size = 7)+
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.5, "cm"),
        strip.background = element_blank(),
        legend.key = element_rect(colour = "white"),
        legend.key.spacing.y = unit(0.4, 'cm'), 
        strip.text.x = element_blank(),
        axis.line.x.bottom = element_line(color = "black")
  )

distance_plot

## combine plots
nmpi2_results <- offspring_genealogy / distance_plot / sel_grad_plot  +
  plot_layout(ncol = 1, heights = c(1, 1, 1)) & 
  plot_annotation(
    tag_levels = "a",        # lower case letters
    tag_prefix = "(", 
    tag_suffix = ")"
  ) & 
  theme(plot.tag = element_text(size = 8, face = "bold"))

path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/Figure4_NMpi2.png"), 
    width =7, height=6, units = "in", res=600)
plot(nmpi2_results)
dev.off()

ggsave(paste0(path, "Manuscript/Figures/Figure4_NMpi2.pdf"),  plot = nmpi2_results, 
       width =7, height=6, units = "in")

# -------- Figure 5 NEMO HI patterns Ar 30   ----

## NEMO starting scenario
all_sim <- fread(paste0(path, "Data_clean/NEMO/Allenwiller_simulated_spatial.txt"))
all_sim$stand = "Allenwiller"
wal_sim <- fread(paste0(path, "Data_clean/NEMO/Waldi_simulated_spatial.txt"))
wal_sim$stand = "Waldi"

all_sim$species = ifelse(grepl("ori" , all_sim$SampleID), "Pure hohenackeriana", "Pure sylvatica")
wal_sim$species = ifelse(grepl("ori" , wal_sim$SampleID), "Pure hohenackeriana", "Pure sylvatica")

# load simulation grid
grid_wal <- vect(paste0(path,"Data_clean/Waldi_grid_4m_5396patches.shp"), crs = "EPSG:3035")
grid_al <- vect(paste0(path,"Data_clean/Allenwiller_grid_r4m_9844patches.shp"), crs = "EPSG:3035")

all_grid <- ggplot(grid_al) +
  geom_spatvector(fill = "white", col = "grey70")+theme_void()+
  geom_point(all_sim,mapping =  aes(x = x, y = y, color = species), size = 0.8)+
  scale_color_manual(values= ColourVector3) + 
  guides(col = "none")+
  theme( plot.title = element_text(hjust = 0.5,vjust = 1, size = 8, face = "bold", margin = margin(b = 0)) )+ 
  labs(col = "")
all_grid

wal_grid <- ggplot(grid_wal) +
  geom_spatvector(fill = "white", col = "grey70")+theme_void()+
  geom_point(wal_sim,mapping =  aes(x = x, y = y, color = species), size = 0.8)+
  scale_color_manual(values= ColourVector3) + 
  guides(col = "none")+
  theme( plot.title = element_text(hjust = 0.5,vjust = 1, size = 8, face = "bold", margin = margin(b = 0)) )+ 
  labs(col = "")
wal_grid

gc()

## data prep
setwd("C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Data_clean/")
grid_adults_2020 <- readRDS("NEMO/HI_data_adults_2020.RDS")
grid_adults_2120 <- readRDS("NEMO/HI_data_adults_2120.RDS")
grid_juv_2020 <- readRDS("NEMO/HI_data_juveniles_2020.RDS")
grid_juv_2120 <- readRDS("NEMO/HI_data_juveniles_2020.RDS")
grid_seed_2020 <- readRDS("NEMO/HI_data_seedlings_2020.RDS")
grid_seed_2120 <- readRDS("NEMO/HI_data_seedlings_2120.RDS")

      # read files with hybridization rate 
      # per circle plot? Simulated_HI_counts_circleplots_2020.RDS
      sim_HI_2020_summary <- readRDS("NEMO/NEMO_HI_2020_summary.RDS")
      sim_HI_2120_summary <- readRDS("NEMO/NEMO_HI_2120_summary.RDS")
              
      # add hybridization rate to grids
      #grid_adults_2020 <- grid_adults_2020 %>% left_join(sim_HI_2020_summary[which(sim_HI_2020_summary$stage == 3), c("stand", "scenario2", "hyb_rate_scenario_stage")],  by = c("stand", "scenario2"))
      #grid_adults_2120 <- grid_adults_2120 %>% left_join(sim_HI_2120_summary[which(sim_HI_2120_summary$stage == 3), c("stand", "scenario2", "hyb_rate_scenario_stage")],  by = c("stand", "scenario2"))
      grid_juv_2020 <- grid_juv_2020 %>% left_join(sim_HI_2020_summary[which(sim_HI_2020_summary$stage == 2), c("stand", "scenario2", "hyb_rate_scenario_stage")],  by = c("stand", "scenario2"))
      grid_juv_2120 <- grid_juv_2120 %>% left_join(sim_HI_2120_summary[which(sim_HI_2120_summary$stage == 2), c("stand", "scenario2", "hyb_rate_scenario_stage")],  by = c("stand", "scenario2"))
      grid_seed_2020 <- grid_seed_2020 %>% left_join(sim_HI_2020_summary[which(sim_HI_2020_summary$stage == 1), c("stand", "scenario2", "hyb_rate_scenario_stage")],  by = c("stand", "scenario2"))
      grid_seed_2120 <- grid_seed_2120 %>% left_join(sim_HI_2120_summary[which(sim_HI_2120_summary$stage == 1), c("stand", "scenario2", "hyb_rate_scenario_stage")],  by = c("stand", "scenario2"))
 
      # change column names (hybridization rate)
      colnames(grid_juv_2020)[17] <- "hybridization_rate"
      colnames(grid_juv_2120)[17] <- "hybridization_rate"
      colnames(grid_seed_2020)[17] <- "hybridization_rate"
      colnames(grid_seed_2120)[17] <- "hybridization_rate"
      
## plotting theme
theme_p <- theme(
  axis.line = element_blank(),
  plot.background = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  plot.margin = margin(t = 8, r = 0, b = 1, l = 1),
  strip.background = element_blank(),
  strip.placement = "outside",
  strip.text.x = element_text(size = 7, margin = margin(b = 1)),
  strip.text.y.left = element_text(angle = 0, size = 7, margin = margin(r =1)),
  plot.title = element_text(size = 7, hjust = 0.5)
)

# function
plot_HI <- function(data_subset, title, sampling_data = NULL) {
  stand_name <- unique(data_subset$stand)
  p <- ggplot(data_subset) +
    geom_sf(aes(fill = dominant_HI, alpha = dominance_strength), color = NA) +
    scale_fill_manual(name = "HI category",
      values = c("HI0" = "#ffcc00", "HI0.5" = "#440154FF", "HI1" = "#25858EFF"),
      labels = c("HI = 0", "0 < HI < 1", "HI = 1")
    ) +
    scale_alpha_continuous(range = c(0, 1)) +
    ggh4x::facet_grid2(dp ~ s, labeller = label_value, switch = "both", render_empty = FALSE) +
    theme_minimal(base_size = 7) +
    geom_text(aes(x = Inf, y = -Inf, label = paste0(round(hybridization_rate, 2))),hjust = 1.1, vjust = -0.5, size = 1.8, inherit.aes = FALSE)+
    labs(title = title, x = "Selection strength (s)", y = "Dispersal distance (dp)", alpha = "Proportion of\ndominant\nHI category") +
    theme_p
  
  if (!is.null(sampling_data)) {
    p <- p + geom_sf(data = sampling_data, shape = 21, color = "grey30", fill = NA, size = 0.5, stroke = 0.3, alpha = 0.6)
  }
  
  return(p)
}

    
    ## add sampling scheme to the offspring (not needed)
    samplings <- c("Waldi_sampling_scheme.shp", "Allenwiller_sampling_scheme.shp")
    stand_names <- c("Waldi", "Allenwiller")
    
    sampling_schemes <- list()
    for (i in seq_along(samplings)) {
      sampling_schemes[[i]] <- read_sf(samplings[i]) %>%
        mutate(stand = stand_names[i])
      sampling_schemes[[i]] <- sampling_schemes[[i]][,-c(3:13)] 
    }
    sampling_all <- bind_rows(sampling_schemes)

# adults 2020 
wal_adults_30 <- plot_HI(subset(grid_adults_2020, stand == "Waldi" & Ar == 30), "")
al_adults_30 <- plot_HI(subset(grid_adults_2020, stand == "Allenwiller" & Ar == 30), "")
wal_adults_50 <- plot_HI(subset(grid_adults_2020, stand == "Waldi" & Ar == 50), "")
al_adults_50 <- plot_HI(subset(grid_adults_2020, stand == "Allenwiller" & Ar == 50), "")

# juveniles 2020  -- removed sampling plots
wal_juv_30 <- plot_HI(subset(grid_juv_2020, stand == "Waldi" & Ar == 30),"")
al_juv_30 <- plot_HI(subset(grid_juv_2020, stand == "Allenwiller" & Ar == 30), "")
wal_juv_50 <- plot_HI(subset(grid_juv_2020, stand == "Waldi" & Ar == 50),"")
al_juv_50 <- plot_HI(subset(grid_juv_2020, stand == "Allenwiller" & Ar == 50), "")

# seedling 2020 
wal_seed_30 <- plot_HI(subset(grid_seed_2020, stand == "Waldi" & Ar == 30),"")
al_seed_30 <- plot_HI(subset(grid_seed_2020, stand == "Allenwiller" & Ar == 30), "")
wal_seed_50 <- plot_HI(subset(grid_seed_2020, stand == "Waldi" & Ar == 50),"")
al_seed_50 <- plot_HI(subset(grid_seed_2020, stand == "Allenwiller" & Ar == 50), "")

# adults plots 2120
wal_adults_30_f <- plot_HI(subset(grid_adults_2120, stand == "Waldi" & Ar == 30), " ")
wal_adults_50_f <- plot_HI(subset(grid_adults_2120, stand == "Waldi" & Ar == 50), " ")
al_adults_30_f <- plot_HI(subset(grid_adults_2120, stand == "Allenwiller" & Ar == 30), " ")
al_adults_50_f <- plot_HI(subset(grid_adults_2120, stand == "Allenwiller" & Ar == 50), " ")

### juvlings plots 2120
wal_juv_30_f <- plot_HI(subset(grid_juv_2120, stand == "Waldi" & Ar == 30)," ")
wal_juv_50_f <- plot_HI(subset(grid_juv_2120, stand == "Waldi" & Ar == 50), "")
al_juv_30_f <- plot_HI(subset(grid_juv_2120, stand == "Allenwiller" & Ar == 30), " ")
al_juv_50_f <- plot_HI(subset(grid_juv_2120, stand == "Allenwiller" & Ar == 50), " ")

### seedlings plots 2120
wal_seed_30_f <- plot_HI(subset(grid_seed_2120, stand == "Waldi" & Ar == 30)," ")
wal_seed_50_f <- plot_HI(subset(grid_seed_2120, stand == "Waldi" & Ar == 50), " ")
al_seed_30_f <- plot_HI(subset(grid_seed_2120, stand == "Allenwiller" & Ar == 30), " ")
al_seed_50_f <- plot_HI(subset(grid_seed_2120, stand == "Allenwiller" & Ar == 50), " ")

# get legends
# fill legend
legend_fill <- cowplot::get_legend(
  wal_adults_30 +
    guides(
      fill = guide_legend(direction = "vertical"),
      alpha = "none"  # disable alpha legend
    ) +
    theme(
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7), 
      legend.key.size = unit(0.6, "cm")
    )
)

# alpha legend
legend_alpha <- cowplot::get_legend(
  wal_adults_30 +
    guides(
      fill = "none",  # disable fill legend
      alpha = guide_legend(direction = "vertical")
    ) +
    theme(
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7), 
      legend.key.size = unit(0.6, "cm")
    )
)

legend_fill <- ggdraw(legend_fill) + theme(plot.margin = margin(0, 0, 0, 0))
legend_alpha <- ggdraw(legend_alpha) + theme(plot.margin = margin(0, 0, 0, 0))

grid::grid.newpage()
gridExtra::grid.arrange(legend_fill)
gridExtra::grid.arrange(legend_alpha)

### composite figure allenwiller ar = 30
figure_hi_allenwiller_ar30 <- ggpubr::ggarrange(
  all_grid + theme( plot.margin = margin(0, 0, 0, 0)),
  al_adults_30+ theme( plot.margin = margin(0, 0, 0, 0),axis.title.x = element_blank(), strip.text.x =element_blank(),legend.position ="none"), 
  al_adults_30_f+ theme( plot.margin = margin(0, 0, 0, 0),axis.title.x = element_blank(), axis.title.y = element_blank(),strip.text.x =element_blank(),legend.position ="none"), 
  legend_fill,
  al_juv_30+ theme( axis.title.x = element_blank(), strip.text.x =element_blank(),legend.position ="none"), 
  al_juv_30_f+ theme( axis.title.x = element_blank(), axis.title.y = element_blank(),strip.text.x =element_blank(),legend.position ="none"), 
  legend_alpha,
  al_seed_30+ theme(strip.text.x =element_text(size = 7),legend.position ="none"), 
  al_seed_30_f+ theme( strip.text.x =element_text(size = 7),axis.title.y = element_blank(),legend.position ="none"), 
  
  nrow = 3, ncol = 3,
  widths = c(1,2,2), 
  align = "hv"
)  

top_labels <- plot_grid(
  ggdraw() + draw_label("1920", hjust = 0.5, size = 10, fontface = "bold"),
  ggdraw() + draw_label("2020", hjust = 0.2, size = 10, fontface = "bold"),
  ggdraw() + draw_label("2120", hjust =0.2, size = 10, fontface = "bold"),
  ncol = 3,
  rel_widths = c(1,2,2) 
)

row_labels <- plot_grid(
  ggdraw() + draw_label("Adults", angle = 270, hjust = 0.5, size = 10),
  ggdraw() + draw_label("Juveniles", angle = 270, hjust = 0.5, size = 10),
  ggdraw() + draw_label("Seedlings", angle = 270, hjust = 0.5, size = 10),
  ncol = 1,
  rel_heights = c(1, 1, 1)  
)

figure_hi_allenwiller_ar30_1 <- plot_grid(
  top_labels,figure_hi_allenwiller_ar30 ,
  nrow = 2,align = "v",
  rel_heights = c(0.02,1), rel_widths =  c(2, 1)
)

figure_hi_allenwiller_ar30_2 <- plot_grid(
  figure_hi_allenwiller_ar30_1,row_labels,
  ncol = 2,
  rel_widths = c( 1,0.02)  
)


path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/Figure5a_allenwiller_ar30.png"), 
    width =7, height=7, units = "in", res=600)
plot(figure_hi_allenwiller_ar30_2)
dev.off()

ggsave(paste0(path, "Manuscript/Figures/Figure5a_allenwiller_ar30.pdf"),  plot = figure_hi_allenwiller_ar30, 
       width =7,height=7,units = "in")

### composite figure waldi ar = 30
figure_hi_waldi_ar30 <- ggpubr::ggarrange(
  wal_grid, 
  wal_adults_30+ theme( axis.title.x = element_blank(), strip.text.x =element_blank(),legend.position ="none"), 
  wal_adults_30_f+ theme( axis.title.x = element_blank(), axis.title.y = element_blank(),strip.text.x =element_blank(),legend.position ="none"), 
  legend_fill,
  wal_juv_30+ theme( axis.title.x = element_blank(), strip.text.x =element_blank(),legend.position ="none"), 
  wal_juv_30_f+ theme( axis.title.x = element_blank(), axis.title.y = element_blank(),strip.text.x =element_blank(),legend.position ="none"), 
  legend_alpha,
  wal_seed_30+ theme(strip.text.x =element_text(size = 7),legend.position ="none"), 
  wal_seed_30_f+ theme( strip.text.x =element_text(size = 7),axis.title.y = element_blank(),legend.position ="none"), 
  
  nrow = 3, ncol = 3,
  widths = c(1,2,2), 
  align = "h"
)  

figure_hi_waldi_ar30_1 <- plot_grid(
  top_labels,figure_hi_waldi_ar30 ,
  nrow = 2,align = "v",
  rel_heights = c(0.02,1), rel_widths =  c(2, 1)
)

figure_hi_waldi_ar30_2 <- plot_grid(
  figure_hi_waldi_ar30_1,row_labels,
  ncol = 2,
  rel_widths = c( 1,0.02)  
)


path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/Figure5b_waldi_ar30.png"), 
    width =7, height=7, units = "in", res=600)
plot(figure_hi_waldi_ar30_2)
dev.off()

ggsave(paste0(path, "Manuscript/Figures/Figure5b_waldi_ar30.pdf"),  plot = figure_hi_waldi_ar30_2, 
       width =7,height=7,units = "in")

# -------- Figure 6. Figure Moran's I  -----
# run script 10 first

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

# add the lable for axis
combined_df$scenario_label <- NA
combined_df$scenario_label[grepl("^flatdispersal", combined_df$scenario2)] <- paste0("Flat dispersal + Ar = ",gsub(".*_a([0-9]+)$", "\\1",combined_df$scenario2[grepl("^flatdispersal", combined_df$scenario2)]))
combined_df$scenario_label[ grepl("^spatial_", combined_df$scenario2)] <- paste0("dp = ", gsub(".*_dp([0-9]+)_.*", "\\1",combined_df$scenario2[ grepl("^spatial_", combined_df$scenario2)])," + Ar = ",  gsub(".*_a([0-9]+)$", "\\1", combined_df$scenario2[ grepl("^spatial_", combined_df$scenario2)]))
combined_df$scenario_label[grepl("^spatialselection", combined_df$scenario2)] <- paste0("dp = ", gsub(".*_dp([0-9]+)_.*", "\\1",combined_df$scenario2[grepl("^spatialselection", combined_df$scenario2)])," + s = ",  gsub(".*_s([0-9]\\.[0-9])_.*", "\\1", combined_df$scenario2[grepl("^spatialselection", combined_df$scenario2)])," + Ar = ",gsub(".*_a([0-9]+)$", "\\1",combined_df$scenario2[grepl("^spatialselection", combined_df$scenario2)]))

# sort the scenarios with increasing dp order - gradient color
levels(unique(factor(combined_df$scenario_label)))
combined_df$scenario_label <- factor(combined_df$scenario_label, levels=c("dp = 20 + Ar = 30", 
                                                                          "dp = 20 + Ar = 50", 
                                                                          "dp = 20 + s = 0.3 + Ar = 30", 
                                                                          "dp = 20 + s = 0.3 + Ar = 50" , 
                                                                          "dp = 20 + s = 0.5 + Ar = 30", 
                                                                          "dp = 20 + s = 0.5 + Ar = 50", 
                                                                          "dp = 20 + s = 0.8 + Ar = 30" , 
                                                                          "dp = 20 + s = 0.8 + Ar = 50", 
                                                                          "dp = 50 + Ar = 30", 
                                                                          "dp = 50 + Ar = 50", 
                                                                          "dp = 50 + s = 0.3 + Ar = 30", 
                                                                          "dp = 50 + s = 0.3 + Ar = 50" , 
                                                                          "dp = 50 + s = 0.5 + Ar = 30", 
                                                                          "dp = 50 + s = 0.5 + Ar = 50", 
                                                                          "dp = 50 + s = 0.8 + Ar = 30" , 
                                                                          "dp = 50 + s = 0.8 + Ar = 50",
                                                                          "dp = 109 + Ar = 30", 
                                                                          "dp = 109 + Ar = 50", 
                                                                          "dp = 109 + s = 0.3 + Ar = 30", 
                                                                          "dp = 109 + s = 0.3 + Ar = 50" , 
                                                                          "dp = 109 + s = 0.5 + Ar = 30", 
                                                                          "dp = 109 + s = 0.5 + Ar = 50", 
                                                                          "dp = 109 + s = 0.8 + Ar = 30" , 
                                                                          "dp = 109 + s = 0.8 + Ar = 50", 
                                                                          "Flat dispersal + Ar = 30", 
                                                                          "Flat dispersal + Ar = 50"))

combined_df <- combined_df %>%
  mutate(
    s_numeric = as.numeric(as.character(s)),
    dp_s_combo = paste0("dp", dp, "_s", s_numeric)
  )


color_map <- c(
  # Flat dispersal (gray)
  "dpFlat\ndispersal_s0"   = "#999999",
  "dpFlat\ndispersal_s0.3" = "#999999",
  "dpFlat\ndispersal_s0.5" = "#999999",
  "dpFlat\ndispersal_s0.8" = "#999999",
  
  # dp = 109 
  "dp109_s0"   = "#FBC15E",
  "dp109_s0.3" = "#F5B042",
  "dp109_s0.5" = "#E69F00",
  "dp109_s0.8" = "#C28500",
  
  # dp = 50 
  "dp50_s0"   = "#B6DFF6",
  "dp50_s0.3" = "#8DCEF1",
  "dp50_s0.5" = "#56B4E9",
  "dp50_s0.8" = "#3091C7",
  
  # dp = 20 
  "dp20_s0"   = "#A1E8DA",
  "dp20_s0.3" = "#5ED3BD",
  "dp20_s0.5" = "#009E73",
  "dp20_s0.8" = "#007256"
)

# my legend (simplified)
## colors of legend
vector_dp <- c(
  "Flat dispersal" = "#999999",   # gray
  "dp = 109"       = "#E69F00",   # orange
  "dp = 50"        = "#56B4E9",   # sky blue
  "dp = 20"        = "#009E73"    # bluish green
)

dummy_df <- data.frame(
  group = factor(names(vector_dp), levels = names(vector_dp)),
  x = 1:4,
  y = rep(1, 4),
  ymin = rep(0.7, 4),
  ymax = rep(1.3, 4)
)

# fake plot fir the legend
dummy_plot <- ggplot(dummy_df, aes(x = x, y = y)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax, color = group), width = 0.2) +
  geom_point(aes(color = group), size = 2) +
  scale_color_manual(values = vector_dp, name = "Dispersal distance (m)") +
  guides(color = guide_legend(override.aes = list(
    shape = 16, size = 2,
    linetype = 1, stroke = 0.5
  ), ncol = 1)) +
  theme_void(base_size = 7) +
  theme(
    legend.position = "right",
    legend.key = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.key.spacing.y = unit(0.3, 'cm'),
    legend.spacing.y = unit(0.05, "cm"),
    legend.margin = margin(2, 2, 2, 2)
  )


legend_dp <- cowplot::get_legend(dummy_plot)
grid.newpage()
grid.draw(legend_dp)

combined_plot <- ggplot(combined_df, aes(x = scenario_label, y = value)) +
  geom_point(aes(y = mean_value, color = dp_s_combo), size = 2) +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value, color = dp_s_combo), width = 0.2) +
  geom_hline(aes(yintercept = observed_value), linetype = "dashed", color = "black") +
  ggh4x::facet_grid2(metric ~ stand, labeller = labeller(stand = c("Allenwiller"="Allenwiller", "Waldi"="Wäldi")), scales = "free_y") +
  scale_color_manual(values = color_map, name = "Dispersal × Selection") +
  labs(x = "Simulation Scenario", y = "Value") +
  theme_bw(base_size = 7) +
  guides(fill = "none", color = "none")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20, unit = "pt"),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 9, margin = margin(b = 5), face = "bold"),
    strip.text.y.left = element_text(angle = 0, size = 9, margin = margin(r = 5)),
    axis.line = element_line(color = 'black'),
    plot.title = element_text(size = 10, hjust = 0.5)
  )

combined_plot <- ggarrange(combined_plot, legend_dp, ncol = 2, widths = c(6,1), align = "hv")

## separate SSE and Morans I
moran_plot <- ggplot(subset(combined_df,metric=="Moran's I"), aes(x = scenario_label, y = value)) +
  geom_point(aes(y = mean_value, color = dp_s_combo), size = 1.7) +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value, color = dp_s_combo), width = 0.2) +
  geom_hline(aes(yintercept = observed_value), linetype = "dashed", color = "black") +
  ggh4x::facet_grid2(metric ~ stand, labeller = labeller(stand = c("Allenwiller"="Allenwiller", "Waldi"="Wäldi")), scales = "free_y") +
  scale_color_manual(values = color_map, name = "Dispersal × Selection") +
  labs(x = "Simulation scenario", y = "Moran's I") +
  theme_bw(base_size = 7) +
  guides(fill = "none", color = "none")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20, unit = "pt"),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 9, margin = margin(b = 5), face = "bold"),
    strip.text.y.right =  element_blank(),
    axis.line = element_line(color = 'black'),
    plot.title = element_text(size = 10, hjust = 0.5)
  )

moran_plot <- ggarrange(moran_plot, legend_dp, ncol = 2, widths = c(6,1), align = "hv")

path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/Figure6_moransI.png"), 
    width =7,height=4, res=600, units = "in")
plot(moran_plot)
dev.off()

ggsave(paste0(path, "Manuscript/Figures/Figure6_moransI.pdf"),  plot = moran_plot, 
       width =7,height=4,units = "in")


# -------- Figure 7. Figure proportions sim obs ----
sim_HI <- readRDS("NEMO/Simulated_HI_counts_circleplots_2020.RDS")
obs_HI <- fread("NEMO/Observed_HI_circleplots.csv")

sim_HI$generation = ifelse(sim_HI$stage == 1, "Offspring", "Juvenile")

sim_summary <- sim_HI %>%
  group_by(stand, generation, scenario2, replicate) %>%
  summarise(
    mean_HI0 = mean(prop_HI0, na.rm = TRUE),
    mean_HI1 = mean(prop_HI1, na.rm = TRUE),
    mean_HI0.5 = mean(prop_HI0.5, na.rm = TRUE),
    .groups = "drop"
  )


obs_summary <- obs_HI %>%
  filter(generation %in% c("Offspring", "Juvenile")) %>%
  group_by(stand, generation) %>%
  summarise(
    mean_HI0 = mean(prop_HI0, na.rm = TRUE),
    mean_HI1 = mean(prop_HI1, na.rm = TRUE),
    mean_HI0.5 = mean(prop_HI0.5, na.rm = TRUE),
    .groups = "drop"
  )


sim_long <- sim_summary %>%
  pivot_longer(cols = starts_with("mean_HI"), names_to = "HI_type", values_to = "mean_prop")

# exclude scenario with flat dispersal from the comparison
sim_long_filtered <- sim_long %>%filter(!str_starts(scenario2, "flatdispersal"))
sim_long_filtered[2,5] <- lapply(sim_long_filtered[2,5], factor) 
obs_long <- obs_summary %>% pivot_longer(cols = starts_with("mean_HI"), names_to = "HI_type", values_to = "mean_prop")
obs_long[2,3] <- lapply(obs_long[2,3], factor) 

sim_long_filtered$generation <- factor(sim_long_filtered$generation , levels=c("Offspring", "Juvenile"))
obs_long$generation <- factor(obs_long$generation , levels=c("Offspring", "Juvenile"))
sim_long_filtered$HI_type <- factor(sim_long_filtered$HI_type , levels= c("mean_HI0",  "mean_HI1" , "mean_HI0.5"))
obs_long$HI_type <- factor(obs_long$HI_type , levels= c("mean_HI0",  "mean_HI1" , "mean_HI0.5"))

boxplot_sim_obs <- ggplot() +
  geom_boxplot(sim_long_filtered, position=position_dodge(width=0.8),
               mapping = aes(x = generation, y = mean_prop, fill = HI_type),
               outlier.shape = NA, alpha = 0.5, width = 0.6) +
  #geom_jitter(aes(color = generation), width = 0.15, size = 0.8, alpha = 0.6) +
  geom_point(data = obs_long, position=position_dodge(width=0.8), aes(x = generation, y = mean_prop,  group = HI_type, col = HI_type),size = 1.8) +
  scale_x_discrete(labels = c("Seedlings\n+ Saplings", "Juveniles"))+
  scale_fill_manual(name = "HI category",
                    values = c("mean_HI0" = "#ffcc00",  "mean_HI1" = "#25858EFF", "mean_HI0.5" = "#440154FF"),
                    labels = c("HI = 0",  "HI = 1","0 < HI < 1"))+
  scale_color_manual(name = "HI category",
                     values = c("mean_HI0" = "#ffcc00", "mean_HI1" = "#25858EFF",  "mean_HI0.5" = "#440154FF"),
                     labels = c("HI = 0",  "HI = 1", "0 < HI < 1"))+
  
  ggh4x::facet_grid2( ~ stand, labeller = labeller(stand = c("Allenwiller"="Allenwiller", "Waldi"="Wäldi"))) +
  labs(y = "Mean proportion", x = NULL, col = NULL) +
  theme_bw(base_size = 7) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.position = "none", 
        strip.text.x = element_text(size = 9, margin = margin(b = 5), face = "bold"),
        legend.key = element_rect(colour = "white"),
        legend.key.spacing.y = unit(0.4, 'cm'), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        plot.title = element_text(hjust = 0.5))

boxplot_sim_obs

## summary
sim_long_filtered %>%
  group_by(stand, generation, HI_type) %>%
  summarise(mean_HI_prop = mean(mean_prop), .groups = "drop")


path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/Figure7_boxplots.png"), 
    width =7, height=3, units = "in", res=600)
plot(boxplot_sim_obs)
dev.off()


ggsave(paste0(path, "Manuscript/Figures/Figure7_boxplots.pdf"),  plot = boxplot_sim_obs, 
       width =7,height=3,units = "in")


# -------- Figure S1. Kriging seedling density -----
setwd(path)
# read kriging data
all_kr <- fread("Data_clean/Allenwiller_kriging_results.csv")
wal_kr <- fread("Data_clean/Waldi_kriging_results.csv")

all_obs <- st_read("Data_clean/Allenwiller_sampling_scheme.shp")
wal_obs <- st_read("Data_clean/Waldi_sampling_scheme.shp")

obs_al <- ggplot() +
  geom_sf(data = all_obs, mapping= aes(fill = CountedC1.),  alpha = 0.8)+
  scale_fill_viridis_c(option = "D", direction = -1) +
  labs(title = "Observed", fill = "Number of seedlings") + 
  theme_void() +
  theme(legend.position = "right",
        axis.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 1,  size = 8, face = "bold", margin = margin(1,1,1,1)), 
        plot.margin = margin(t = 2, r = 5, b = 2, l = 5)
  ) +
  annotation_scale(location = "bl", unit_category= "metric",  width_hint = 0.5, style = "ticks", line_width = 1.2)+
  ggtitle("Allenwiller (Observed)") 

obs_wal <- ggplot() +
  geom_sf(data = wal_obs, mapping= aes(fill = NClass1Cla),  alpha = 0.8)+
  scale_fill_viridis_c(option = "D", direction = -1) +
  labs(title = "Observed", fill = "Number of seedlings") + 
  theme_void() +
  theme(legend.position = "right",
        axis.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 1,  size = 8, face = "bold", margin = margin(1,1,1,1)), 
        plot.margin = margin(t = 2, r = 5, b = 2, l = 5)
  ) +
  annotation_scale(location = "bl", unit_category= "metric",  width_hint = 0.5, style = "ticks", line_width = 1.2)+
  ggtitle("Wäldi (Observed)") 

# get the extent for the kriging plot
bbox_obs <- st_bbox(all_obs)
bbox_obs2 <- st_bbox(wal_obs)

# predicted
krig_al <- ggplot(all_kr, aes(x = x_centroid, y = y_centroid, fill = predicted)) +
  geom_tile() + 
  geom_sf(data = all_obs, mapping= aes(x = x, y = y), fill = NA, col = "black")+
  coord_sf(xlim = c(bbox_obs["xmin"], bbox_obs["xmax"]),
           ylim = c(bbox_obs["ymin"], bbox_obs["ymax"])) +
  scale_fill_viridis_c(option = "D", direction = -1) +
  labs(title = "Observed", fill = "Number of seedlings") + 
  theme_void() +
  theme(legend.position = "right",
        axis.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 1,  size = 8, face = "bold", margin = margin(1,1,1,1)), 
        plot.margin = margin(t = 2, r = 5, b = 2, l = 5)
  ) +
  annotation_scale(location = "bl", unit_category= "metric",  width_hint = 0.5, style = "ticks", line_width = 1.2)+
  ggtitle("Allenwiller (Predicted)") 

krig_wal <- ggplot() +
  geom_tile(wal_kr, mapping = aes(x = x_centroid, y = y_centroid, fill = predicted)) +  
  geom_sf(data = wal_obs,   col = "black", fill = "NA")+
  coord_sf(xlim = c(bbox_obs2["xmin"], bbox_obs2["xmax"]),
           ylim = c(bbox_obs2["ymin"], bbox_obs2["ymax"])) +
  scale_fill_viridis_c(option = "D", direction = -1) +
  labs(title = "Observed", fill = "Number of seedlings") + 
  theme_void() +
  theme(legend.position = "right",
        axis.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 1,  size = 8, face = "bold", margin = margin(1,1,1,1)),
        plot.margin = margin(t = 2, r = 5, b = 2, l = 5)
  ) +
  annotation_scale(location = "bl", unit_category= "metric",  width_hint = 0.5, style = "ticks", line_width = 1.2)+
  ggtitle("Wäldi (Predicted)") 


krig_plot <- ggarrange(obs_al, krig_al, obs_wal, krig_wal, ncol = 2, nrow = 2, common.legend = T, legend = "bottom", align = "hv")
png(filename=paste0(path, "Manuscript/Figures/FigureS1_kriging.png"),  
    width =5, height=6, units = "in", res=600)
plot(krig_plot)
dev.off()

ggsave(paste0(path, "Manuscript/Figures/FigureS1_kriging.pdf"),  plot = krig_plot, 
       width =5, height=6,units = "in")

# -------- Figure S2. informative score loci --------

# take allele frequencies from Supp material Kurz et al 2023
# for each locus, calculating (|peak size Syl - peak size Ori|)
# dividing for the total bin size (e.g., the difference between the smallest and largest allele) to normalize on the tot variabiliy of the locus
informative_score <- c(
  casolfagus_29 = (138 - 134) / (150 - 130),
  concat14 = (197 - 195) / (211 - 177),
  csolfagus_05 = (197 - 195) / (182 - 160),
  csolfagus_06 = (217 - 217) / (237 - 203),
  csolfagus_19 = (173 - 173) / (189 - 159),
  csolfagus_31 = (127 - 113) / (131 - 95),
  DE576 = (227 - 221) / (233 - 209),
  DUKCT = (89 - 81) / (97 - 77),
  DZ447 = (191 - 191) / (243 - 183),
  EEU75 = (97 - 91) / (119 - 87),
  EJV8T = (153 - 147) / (163 - 145),
  EMILY = (148 - 136) / (168 - 136),
  ERHBI = (167 - 161) / (183 - 153),
  FS1_15 = (115 - 95) / (139 - 83),
  sfc_0036 = (100 - 96) / (118 - 90),
  sfc_1143 = (129 - 121) / (147 - 99)
)

# convert to df
informative_score <- data.frame(Locus = names(informative_score),Score = as.numeric(informative_score))
informative_score <- informative_score[order(informative_score$Score, decreasing = TRUE), ]

# most informative 
most_inf <- informative_score$Locus[1:8]
less_inf <- informative_score$Locus[9:16]
ssr_panel <- ggplot(informative_score, aes(x = reorder(Locus, Score), y = Score)) +
  geom_bar(stat = "identity", fill = "black") +
  coord_flip() +  # Flip coordinates for easier reading
  labs(
    x = "Locus",
    y = "Informative Score"
  ) +
  theme_bw(base_size = 10) + 
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom", 
        strip.text = element_blank(),
        legend.key = element_rect(colour = "white"),
        legend.key.spacing.y = unit(0.4, 'cm'), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        plot.title = element_text(hjust = 0.5))
ssr_panel

png(filename=paste0(path, "Manuscript/Figures/FigureS2_panel_informativescore.png"), 
    width =5, height=4, units = "in", res=600)
plot(ssr_panel)
dev.off()

ggsave(paste0(path, "Manuscript/Figures/FigureS2_panel_informativescore.pdf"),  plot = ssr_panel, 
       width =5, height=4,units = "in")



# -------- Figure S3. SSR Panel analysis -----------
rm(list=ls())
path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"

accuracy_res <- fread(paste0(path, "Data_clean/Hybriddetective/Accuracy_results.csv"))
efficiency_res <- fread(paste0(path, "Data_clean/Hybriddetective/Efficiency_results.csv"))
power_res <- fread(paste0(path, "Data_clean/Hybriddetective/Power_results.csv"))

accuracy_res <- accuracy_res %>%
  rename(Value = Accuracy)

efficiency_res <- efficiency_res %>%
  rename(Value = Efficiency)

power_res <- power_res %>%
  rename(Value = Power)

# Add a Measure column
accuracy_res$Measure <- "Accuracy"
efficiency_res$Measure <- "Efficiency"
power_res$Measure <- "Power"

# Combine datasets
all_res <- bind_rows(accuracy_res, efficiency_res, power_res)

# correcting names
all_res$panel <- gsub("informative", "diagnostic", all_res$panel)
all_res$panel <- factor(all_res$panel,levels = c("all loci", "8 most diagnostic loci", "8 less diagnostic loci"), 
                        labels = c("All loci", "Eight most diagnostic loci", "Eight less diagnostic loci"))
all_res$colour <- gsub("Pure1", "P1", all_res$colour)
all_res$colour <- gsub("Pure2", "P2", all_res$colour)


# Genotype_class (harmonized)
all_res$Genotype_class <- factor(all_res$colour, 
                                 levels = c("P1", "P2", "F1", "F2", "BC1", "BC2"),
                                 labels = c("Pure sylvatica", "Pure hohenackeriana", "F1", "F2", 
                                            "BC sylvatica", "BC hohenackeriana"))


ColourVector = c("#ffcc00", "#482173FF", "#25858EFF", "#85D54AFF", "orange", "#9C99C8")

# allenwiler
all_plot <- ggplot(subset(all_res,stand == "Allenwiller"), aes(x = PostProb, y = Value)) +
  geom_ribbon(aes(x = PostProb, ymin = sd.lower, ymax = sd.upper,  fill = Genotype_class,group = Genotype_class), alpha = 0.2) +
  geom_line(aes(colour = Genotype_class), lwd = 0.5) +
  geom_vline(xintercept = 0.8, linetype = 2, lwd = 0.5) +
  facet_grid(Measure ~   panel) +
  scale_colour_manual(values = ColourVector) +
  scale_fill_manual(values = ColourVector) +
  labs(x = "Posterior probability threshold", 
       y = expression("Metric "%+-%"sd"), 
       col = "", fill = "") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 7)+
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", 
        legend.key = element_rect(colour = "white"),
        legend.text = element_text(size = 5),      # Smaller legend text
        legend.title = element_text(size = 7),     # Smaller legend title
        legend.key.size = unit(0.3, "cm"),         # Smaller legend keys
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.spacing.y = unit(0.05, "cm"),        # Less spacing
        legend.key.spacing.y = unit(0.05, 'cm'), 
        legend.margin = margin(2, 2, 2, 2),
        strip.background = element_rect(colour = "white", fill = "white"), 
        plot.title = element_text(hjust = 0.5,vjust = 1, size = 8, face = "bold", margin = margin(b = 0)))+
  guides(fill = "none",color=guide_legend(ncol=1, title = "Genotype class"))+
  ggtitle("Allenwiller")
all_plot

# waldi
wal_plot <- ggplot(subset(all_res,stand == "Waldi"), aes(x = PostProb, y = Value)) +
  geom_ribbon(aes(x = PostProb, ymin = sd.lower, ymax = sd.upper,  fill = Genotype_class,group = Genotype_class), alpha = 0.2) +
  geom_line(aes(colour = Genotype_class), lwd = 0.5) +
  geom_vline(xintercept = 0.8, linetype = 2, lwd = 0.5) +
  facet_grid(Measure ~   panel) +
  scale_colour_manual(values = ColourVector) +
  scale_fill_manual(values = ColourVector) +
  labs(x = "Posterior probability threshold", 
       y = expression("Metric "%+-%"sd"), 
       col = "", fill = "") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 7)+
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", 
        legend.key = element_rect(colour = "white"),
        legend.text = element_text(size = 5),      # Smaller legend text
        legend.title = element_text(size = 7),     # Smaller legend title
        legend.key.size = unit(0.3, "cm"),         # Smaller legend keys
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.spacing.y = unit(0.05, "cm"),        # Less spacing
        legend.key.spacing.y = unit(0.05, 'cm'), 
        legend.margin = margin(2, 2, 2, 2),
        strip.background = element_rect(colour = "white", fill = "white"), 
        plot.title = element_text(hjust = 0.5,vjust = 1, size = 8, face = "bold", margin = margin(b = 0)))+
  guides(fill = "none",color=guide_legend(ncol=1, title = "Genotype class"))+
  ggtitle("Wäldi")
wal_plot


panel_power <- ggpubr::ggarrange(
  all_plot, wal_plot, 
  nrow = 2, ncol = 1,
  legend = "right",
  common.legend = TRUE, 
  align = "v"  # increase space between the rows
)  

panel_power

png(filename=paste0(path, "Manuscript/Figures/FigureS3_panel_power.png"), 
    width =7,height=7,units = "in",  res=600)
plot(panel_power)
dev.off()

ggsave(paste0(path, "Manuscript/Figures/FigureS3_panel_power.pdf"),  plot = panel_power, 
       width =7,height=7,units = "in")

# -------- Figure S5. NEMO result demographic test -------

# NEMO demographic parameters test results
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

# calculate area of sampling
patch_area_m2 <- 16
avg_sampling_area_per_stand <- combined_fem_data %>%
  distinct(stand, replicate, patch) %>%
  group_by(stand, replicate) %>%
  summarise(n_patches = n(), .groups = "drop") %>%
  mutate(sampling_area = n_patches * patch_area_m2) %>%
  group_by(stand) %>%
  summarise(avg_sampling_area = mean(sampling_area))

sampling_area_sim_wal <- avg_sampling_area_per_stand %>% filter(stand == "Waldi") %>%  pull(avg_sampling_area)
sampling_area_sim_al <- avg_sampling_area_per_stand %>% filter(stand == "Allenwiller") %>%pull(avg_sampling_area)

# normalize numbers on the sampling area 
summary_data_agg$N_agg_stage_scaled <- ifelse(summary_data_agg$stand == "Waldi", summary_data_agg$N_agg_stage/sampling_area_sim_wal, summary_data_agg$N_agg_stage/sampling_area_sim_al)
summary_data_agg$N_agg_stage_scaled <- ifelse(summary_data_agg$stand == "Waldi", summary_data_agg$N_agg_stage/sampling_area_sim_wal, summary_data_agg$N_agg_stage/sampling_area_sim_al)

#observed values dataframe (load from script Number 8)

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

vector_age <- c("Stage 1" = "#009E73",
               "Stage 2" = "#56B4E9",
               "Stage 3" = "#E69F00")

pop_trend_w <- ggplot() +
  geom_line(data = subset(summary_data_agg, stand == "Waldi"), aes(x = generation, y = N_agg_stage_scaled,  group = interaction(agg_stage, replicate),color = agg_stage,linetype = data_source), linewidth= 0.2) +
  geom_hline(data = subset(obs_val, stand == "Waldi"),aes(yintercept = N_obs_scaled,color = agg_stage,linetype = data_source), linewidth= 0.3) +
  scale_linetype_manual(name = "Data type",values = c("Observed data" = "dashed","Simulated data" = "solid")) +
  scale_color_manual(name = "Nemo-age age stage",values = vector_age) +
  theme_bw(base_size = 7)+
  ggtitle("Wäldi")+
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", 
        legend.key = element_rect(colour = "white"),
        legend.text = element_text(size = 5),      
        legend.title = element_text(size = 5),    
        legend.key.size = unit(0.3, "cm"),         
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.spacing.y = unit(0.05, "cm"),       
        legend.key.spacing.y = unit(0.05, 'cm'), 
        legend.margin = margin(2, 2, 2, 2),
        strip.background = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size = 5),
        plot.title = element_text(hjust = 0.5)
  )+
  facet_grid(b + a ~ k + s1,labeller = labeller(
    k = function(x) paste0("K = ", x),
    s1 = function(x) paste0("S1 = ", x),
    a = function(x) paste0("Ar = ", x),
    b = function(x) paste0("b = ", x))) +
  labs(x = "Years",y = "Number of individuals scaled for sampled area",color = "Nemo-age age stage") +
  scale_x_continuous(limits = c(0,150))+
  scale_y_continuous(labels = function(x) format(x, big.mark = "", scientific = FALSE)) +
  scale_y_sqrt() # square root transformation to better visualization
pop_trend_w

pop_trend_a <- ggplot() +
  geom_line(data = subset(summary_data_agg, stand == "Allenwiller"), aes(x = generation, y = N_agg_stage_scaled, group = interaction(agg_stage, replicate),color = agg_stage,linetype = data_source), linewidth= 0.2) +
  geom_hline(data = subset(obs_val, stand == "Allenwiller"),aes(yintercept = N_obs_scaled,color = agg_stage,linetype = data_source), linewidth= 0.3) +
  scale_linetype_manual(name = "Data type",values = c("Observed data" = "dashed","Simulated data" = "solid")) +
  scale_color_manual(name = "Nemo-age age stage",values = vector_age) +
  theme_bw(base_size = 7)+ggtitle("Allenwiller")+
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 5),
        strip.background = element_rect(colour = "black", fill = "white"),
        legend.position = "right", 
        legend.key = element_rect(colour = "white"),
        legend.text = element_text(size = 5),      # Smaller legend text
        legend.title = element_text(size = 7),     # Smaller legend title
        legend.key.size = unit(0.3, "cm"),         # Smaller legend keys
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.spacing.y = unit(0.05, "cm"),        # Less spacing
        legend.key.spacing.y = unit(0.05, 'cm'), 
        legend.margin = margin(2, 2, 2, 2),
        plot.title = element_text(hjust = 0.5)
  )+
  facet_grid(b + a ~ k + s1,labeller = labeller(
    k = function(x) paste0("K = ", x),
    s1 = function(x) paste0("S1 = ", x),
    a = function(x) paste0("Ar = ", x),
    b = function(x) paste0("b = ", x))) +
  labs(x = "Years",y = "Number of Individuals scaled for sampled area",color = "Nemo Stage") +
  scale_x_continuous(limits = c(0,150))+
  scale_y_continuous(labels = function(x) format(x, big.mark = "", scientific = FALSE)) +
  scale_y_sqrt() # square root transformation to better visualization
pop_trend_a


## merge all panels
Fig_nemodemographic <- ggpubr::ggarrange(
  pop_trend_a, pop_trend_w, 
  nrow = 2,legend="bottom",align = "v", common.legend = T
) 

Fig_nemodemographic

path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/FigureS5_nemodemographic.png"), 
    width =7,height=7, res=600, units = "in")
plot(Fig_nemodemographic)
dev.off()

ggsave(paste0(path, "Manuscript/Figures/FigureS5_nemodemographic.pdf"),  plot = Fig_nemodemographic, 
       width =7,height=7, units = "in")

# -------- Figure S6 NEMO HI patterns Ar 50   ----

### composite figure allenwiller ar = 50
figure_hi_allenwiller_ar50 <- ggpubr::ggarrange(
  all_grid, 
  al_adults_50+ theme( axis.title.x = element_blank(), strip.text.x =element_blank(),legend.position ="none"), 
  al_adults_50_f+ theme( axis.title.x = element_blank(), axis.title.y = element_blank(),strip.text.x =element_blank(),legend.position ="none"), 
  legend_fill,
  al_juv_50+ theme( axis.title.x = element_blank(), strip.text.x =element_blank(),legend.position ="none"), 
  al_juv_50_f+ theme( axis.title.x = element_blank(), axis.title.y = element_blank(),strip.text.x =element_blank(),legend.position ="none"), 
  legend_alpha,
  al_seed_50+ theme(strip.text.x =element_text(size = 7),legend.position ="none"), 
  al_seed_50_f+ theme( strip.text.x =element_text(size = 7),axis.title.y = element_blank(),legend.position ="none"), 
  
  nrow = 3, ncol = 3,
  widths = c(1,2,2), 
  align = "h"
)  

top_labels <- plot_grid(
  ggdraw() + draw_label("1920", hjust = 0.5, size = 10, fontface = "bold"),
  ggdraw() + draw_label("2020", hjust = 0.2, size = 10, fontface = "bold"),
  ggdraw() + draw_label("2120", hjust =0.2, size = 10, fontface = "bold"),
  ncol = 3,
  rel_widths = c(1,2,2) 
)

row_labels <- plot_grid(
  ggdraw() + draw_label("Adults", angle = 270, hjust = 0.5, size = 10),
  ggdraw() + draw_label("Juveniles", angle = 270, hjust = 0.5, size = 10),
  ggdraw() + draw_label("Seedlings", angle = 270, hjust = 0.5, size = 10),
  ncol = 1,
  rel_heights = c(1, 1, 1)  
)

figure_hi_allenwiller_ar50_1 <- plot_grid(
  top_labels,figure_hi_allenwiller_ar50 ,
  nrow = 2,align = "v",
  rel_heights = c(0.02,1), rel_widths =  c(2, 1)
)

figure_hi_allenwiller_ar50_2 <- plot_grid(
  figure_hi_allenwiller_ar50_1,row_labels,
  ncol = 2,
  rel_widths = c( 1,0.02)  
)


path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/FigureS6a_allenwiller_ar50.png"), 
    width =7, height=7, units = "in", res=600)
plot(figure_hi_allenwiller_ar50_2)
dev.off()

ggsave(paste0(path, "Manuscript/Figures/FigureS6a_allenwiller_ar50.pdf"),  plot = figure_hi_allenwiller_ar50_2, 
       width =7,height=7,units = "in")


### composite figure waldi ar = 50
figure_hi_waldi_ar50 <- ggpubr::ggarrange(
  wal_grid, 
  wal_adults_50+ theme( axis.title.x = element_blank(), strip.text.x =element_blank(),legend.position ="none"), 
  wal_adults_50_f+ theme( axis.title.x = element_blank(), axis.title.y = element_blank(),strip.text.x =element_blank(),legend.position ="none"), 
  legend_fill,
  wal_juv_50+ theme( axis.title.x = element_blank(), strip.text.x =element_blank(),legend.position ="none"), 
  wal_juv_50_f+ theme( axis.title.x = element_blank(), axis.title.y = element_blank(),strip.text.x =element_blank(),legend.position ="none"), 
  legend_alpha,
  wal_seed_50+ theme(strip.text.x =element_text(size = 7),legend.position ="none"), 
  wal_seed_50_f+ theme( strip.text.x =element_text(size = 7),axis.title.y = element_blank(),legend.position ="none"), 
  
  nrow = 3, ncol = 3,
  widths = c(1,2,2),
  align = "h"
)  

figure_hi_waldi_ar50_1 <- plot_grid(
  top_labels,figure_hi_waldi_ar50 ,
  nrow = 2,align = "v",
  rel_heights = c(0.02,1), rel_widths =  c(2, 1)
)

figure_hi_waldi_ar50_2 <- plot_grid(
  figure_hi_waldi_ar50_1,row_labels,
  ncol = 2,
  rel_widths = c( 1,0.02)  
)


path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/FigureS6b_waldi_ar50.png"), 
    width =7, height=7, units = "in", res=600)
plot(figure_hi_waldi_ar50_2)
dev.off()

ggsave(paste0(path, "Manuscript/Figures/FigureS6b_waldi_ar50.pdf"),  plot = figure_hi_waldi_ar50_2, 
       width =7,height=7,units = "in")

# -------- Figure S7A_B. difference P_hi(SIM - OBS) in the circle plots -----------
sim_HI <- readRDS("NEMO/Simulated_HI_counts_circleplots_2020.RDS")
obs_HI <- fread("NEMO/Observed_HI_circleplots.csv")
sim_HI <- sim_HI %>%  mutate(
  generation = case_when(
    stage == 3~"Adult",
    stage == 2~"Juvenile",
    stage == 1~"Offspring"
  ))
# mean proportion across replicates for each circle plot and stage
mean_sim_HI <- sim_HI %>% group_by(scenario2, stand,generation, CirclePlotID) %>%
  summarize(prop_HI1 = mean(prop_HI1, na.rm = T), 
            prop_HI0.5 = mean(prop_HI0.5, na.rm = T),
            prop_HI0 = mean(prop_HI0, na.rm = T))

# merge and calculate difference
merged_HI <- obs_HI %>%  inner_join(mean_sim_HI, by = c("CirclePlotID", "stand", "generation"), suffix = c("_obs", "_sim"))
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

theme_p <- theme(
  axis.line = element_blank(),
  plot.background = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  plot.margin = margin(t = 8, r = 0, b = 1, l = 1),
  strip.background = element_blank(),
  strip.placement = "outside",
  strip.text.x = element_text(size = 7, margin = margin(b = 1)),
  strip.text.y.left = element_text(angle = 0, size = 7, margin = margin(r =1)),
  plot.title = element_text(size = 7, hjust = 0.5)
)


plot_HI_diff <- function(data_subset, title, diff_prop, HI_value, zoom = TRUE) {
  stand_name <- unique(data_subset$stand)
  
  data_subset <- data_subset %>%
    mutate(
      dp_label = as.character(dp),
      dp_label = ifelse(dp_label == "1000", "Flat\ndispersal", dp_label),
      dp_label = factor(dp_label, levels = c("Flat\ndispersal", "20", "50", "109"))
    )
  
  base_plot <- ggplot() +
   geom_sf( data = background_grid %>%filter(stand == stand_name, has_data == TRUE),fill = "grey90", color = NA) +
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
  
  # zoom in for Allenwiller only
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

# difference for F1
wal_seed_diff_30 <- plot_HI_diff(subset(HI_data_grid, stand == "Waldi" & Ar == 30 & generation == "Offspring"), "Ar = 30", diff_prop_HI1, 1, zoom = T)
wal_seed_diff_50 <- plot_HI_diff(subset(HI_data_grid, stand == "Waldi" & Ar == 50 & generation == "Offspring"), "Ar = 50", diff_prop_HI1, 1, zoom = T)
al_seed_diff_30 <- plot_HI_diff(subset(HI_data_grid, stand == "Allenwiller" & Ar == 30 & generation == "Offspring"), "Ar = 30", diff_prop_HI1, 1, zoom = T)
al_seed_diff_50 <- plot_HI_diff(subset(HI_data_grid, stand == "Allenwiller" & Ar == 50 & generation == "Offspring"), "Ar = 50", diff_prop_HI1, 1, zoom = T)

wal_juv_diff_30 <- plot_HI_diff(subset(HI_data_grid, stand == "Waldi" & Ar == 30 & generation == "Juvenile"), "Ar = 30", diff_prop_HI1, 1, zoom = T)
wal_juv_diff_50 <- plot_HI_diff(subset(HI_data_grid, stand == "Waldi" & Ar == 50 & generation == "Juvenile"), "Ar = 50", diff_prop_HI1, 1, zoom = T)
al_juv_diff_30 <- plot_HI_diff(subset(HI_data_grid, stand == "Allenwiller" & Ar == 30 & generation == "Juvenile"), "Ar = 30", diff_prop_HI1, 1, zoom = T)
al_juv_diff_50 <- plot_HI_diff(subset(HI_data_grid, stand == "Allenwiller" & Ar == 50 & generation == "Juvenile"), "Ar = 50", diff_prop_HI1, 1, zoom = T)

## difference for advanced generation hybrids
wal_seed_diff_30_2  <- plot_HI_diff(subset(HI_data_grid, stand == "Waldi" & Ar == 30& generation == "Offspring"), "Ar = 30", diff_prop_HI0.5, 0.5)
wal_seed_diff_50_2 <- plot_HI_diff(subset(HI_data_grid, stand == "Waldi" & Ar == 50& generation == "Offspring"), "Ar = 50", diff_prop_HI0.5, 0.5)
al_seed_diff_30_2 <- plot_HI_diff(subset(HI_data_grid, stand == "Allenwiller" & Ar == 30& generation == "Offspring"), "Ar = 30", diff_prop_HI0.5, 0.5)
al_seed_diff_50_2 <- plot_HI_diff(subset(HI_data_grid, stand == "Allenwiller" & Ar == 50& generation == "Offspring"), "Ar = 50", diff_prop_HI0.5, 0.5)

wal_juv_diff_30_2  <- plot_HI_diff(subset(HI_data_grid, stand == "Waldi" & Ar == 30& generation == "Juvenile"), "Ar = 30", diff_prop_HI0.5, 0.5)
wal_juv_diff_50_2 <- plot_HI_diff(subset(HI_data_grid, stand == "Waldi" & Ar == 50& generation == "Juvenile"), "Ar = 50", diff_prop_HI0.5, 0.5)
al_juv_diff_30_2 <- plot_HI_diff(subset(HI_data_grid, stand == "Allenwiller" & Ar == 30& generation == "Juvenile"), "Ar = 30", diff_prop_HI0.5, 0.5)
al_juv_diff_50_2 <- plot_HI_diff(subset(HI_data_grid, stand == "Allenwiller" & Ar == 50& generation == "Juvenile"), "Ar = 50", diff_prop_HI0.5, 0.5)


HI_diff_plot <- ggpubr::ggarrange(
  al_seed_diff_30+ theme(axis.title.x = element_blank(), strip.text.x =element_blank(),legend.position ="none"), 
  #al_juv_diff_30+ theme(axis.title = element_blank(),strip.text.x =element_blank(), strip.text.y.left = element_blank(),legend.position ="none"),
  wal_seed_diff_30 + theme(legend.position ="none"),
  wal_juv_diff_30+ theme(axis.title.y=element_blank(),  strip.text.y.left = element_blank(),legend.position ="none" ),
  
  
  nrow = 2, ncol = 2,
  common.legend = T, legend = "right",
  widths = c(1, 1), heights = c(1, 1),
  labels = c("Allenwiller", "", "Wäldi", ""),
  label.x = 0,  label.y = 1.02,
  font.label = list(size = 12, face = "bold"),
  align = "hv")
HI_diff_plot

path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/FigureS8A_HI1_diff.png"), 
    width =7.0866, height=7.0866, units = "in", res=600)
plot(HI_diff_plot)
dev.off()

ggsave(paste0(path, "Manuscript/Figures/FigureS7A_HI1_diff.pdf"),  plot = HI_diff_plot, 
       width =7.0866, height=7.0866, units = "in")

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
png(filename=paste0(path, "Manuscript/Figures/FigureS8B_HI0.5_diff.png"), 
    width =7.0866, height=7.0866, units = "in", res=600)
plot(HI_diff_plot)
dev.off()


ggsave(paste0(path, "Manuscript/Figures/FigureS7B_HI0.5_diff.pdf"),  plot = HI_diff_plot, 
       width =7.0866, height=7.0866, units = "in")

# -------- Figure S7 (C). plot SSE circle plots ----------
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
png(filename=paste0(path, "Manuscript/Figures/FigureS7C_SSE_grid.png"), 
    width =7.0866, height=7.0866, units = "in", res=600)
plot(SSE_plot)
dev.off()

ggsave(paste0(path, "Manuscript/Figures/FigureS7C_SSE_grid.pdf"),  plot = SSE_plot, 
       width =7.0866, height=7.0866, units = "in")



# -------- Figure X. Figure Nemo SSE  -----
# run script 10 first

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

# add the lable for axis
combined_df$scenario_label <- NA
combined_df$scenario_label[grepl("^flatdispersal", combined_df$scenario2)] <- paste0("Flat dispersal + Ar = ",gsub(".*_a([0-9]+)$", "\\1",combined_df$scenario2[grepl("^flatdispersal", combined_df$scenario2)]))
combined_df$scenario_label[ grepl("^spatial_", combined_df$scenario2)] <- paste0("dp = ", gsub(".*_dp([0-9]+)_.*", "\\1",combined_df$scenario2[ grepl("^spatial_", combined_df$scenario2)])," + Ar = ",  gsub(".*_a([0-9]+)$", "\\1", combined_df$scenario2[ grepl("^spatial_", combined_df$scenario2)]))
combined_df$scenario_label[grepl("^spatialselection", combined_df$scenario2)] <- paste0("dp = ", gsub(".*_dp([0-9]+)_.*", "\\1",combined_df$scenario2[grepl("^spatialselection", combined_df$scenario2)])," + s = ",  gsub(".*_s([0-9]\\.[0-9])_.*", "\\1", combined_df$scenario2[grepl("^spatialselection", combined_df$scenario2)])," + Ar = ",gsub(".*_a([0-9]+)$", "\\1",combined_df$scenario2[grepl("^spatialselection", combined_df$scenario2)]))

# sort the scenarios with increasing dp order - gradient color
levels(unique(factor(combined_df$scenario_label)))
combined_df$scenario_label <- factor(combined_df$scenario_label, levels=c("dp = 20 + Ar = 30", 
                                                                          "dp = 20 + Ar = 50", 
                                                                          "dp = 20 + s = 0.3 + Ar = 30", 
                                                                          "dp = 20 + s = 0.3 + Ar = 50" , 
                                                                          "dp = 20 + s = 0.5 + Ar = 30", 
                                                                          "dp = 20 + s = 0.5 + Ar = 50", 
                                                                          "dp = 20 + s = 0.8 + Ar = 30" , 
                                                                          "dp = 20 + s = 0.8 + Ar = 50", 
                                                                          "dp = 50 + Ar = 30", 
                                                                          "dp = 50 + Ar = 50", 
                                                                          "dp = 50 + s = 0.3 + Ar = 30", 
                                                                          "dp = 50 + s = 0.3 + Ar = 50" , 
                                                                          "dp = 50 + s = 0.5 + Ar = 30", 
                                                                          "dp = 50 + s = 0.5 + Ar = 50", 
                                                                          "dp = 50 + s = 0.8 + Ar = 30" , 
                                                                          "dp = 50 + s = 0.8 + Ar = 50",
                                                                          "dp = 109 + Ar = 30", 
                                                                          "dp = 109 + Ar = 50", 
                                                                          "dp = 109 + s = 0.3 + Ar = 30", 
                                                                          "dp = 109 + s = 0.3 + Ar = 50" , 
                                                                          "dp = 109 + s = 0.5 + Ar = 30", 
                                                                          "dp = 109 + s = 0.5 + Ar = 50", 
                                                                          "dp = 109 + s = 0.8 + Ar = 30" , 
                                                                          "dp = 109 + s = 0.8 + Ar = 50", 
                                                                          "Flat dispersal + Ar = 30", 
                                                                          "Flat dispersal + Ar = 50"))

combined_df <- combined_df %>%
  mutate(
    s_numeric = as.numeric(as.character(s)),
    dp_s_combo = paste0("dp", dp, "_s", s_numeric)
  )


color_map <- c(
  # Flat dispersal have same color(gray)
  "dpFlat\ndispersal_s0"   = "#636363",
  "dpFlat\ndispersal_s0.3" = "#636363",
  "dpFlat\ndispersal_s0.5" = "#636363",
  "dpFlat\ndispersal_s0.8" = "#636363",
  
  # dp = 109 (green)
  "dp109_s0"   = "#a1d99b",
  "dp109_s0.3" = "#74c476",
  "dp109_s0.5" = "#31a354",
  "dp109_s0.8" = "#006d2c",
  
  # dp = 50 (yellow-orange)
  "dp50_s0"   = "#F9D476",
  "dp50_s0.3" = "#F3A900",
  "dp50_s0.5" = "#E47400",
  "dp50_s0.8" = "#923B00",
  
  # dp = 20 (red)
  "dp20_s0"   = "#fcbba1",
  "dp20_s0.3" = "#fc9272",
  "dp20_s0.5" = "#fb6a4a",
  "dp20_s0.8" = "#cb181d"
)

# create my legend
## colors of legend
vector_dp <- c("Flat dispersal" = "#969696",
               "dp = 109" = "#31a354",
               "dp = 50" = "#fec44f",
               "dp = 20" = "#fb6a4a")

dummy_df <- data.frame(
  group = factor(names(vector_dp), levels = names(vector_dp)),
  x = 1:4,
  y = rep(1, 4),
  ymin = rep(0.7, 4),
  ymax = rep(1.3, 4)
)

# fake plot fir the legend
dummy_plot <- ggplot(dummy_df, aes(x = x, y = y)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax, color = group), width = 0.2) +
  geom_point(aes(color = group), size = 2) +
  scale_color_manual(values = vector_dp, name = "") +
  guides(color = guide_legend(override.aes = list(
    shape = 16, size = 2,
    linetype = 1, stroke = 0.5
  ), ncol = 1)) +
  theme_void(base_size = 7) +
  theme(
    legend.position = "right",
    legend.key = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.key.spacing.y = unit(0.3, 'cm'),
    legend.spacing.y = unit(0.05, "cm"),
    legend.margin = margin(2, 2, 2, 2)
  )


legend_dp <- cowplot::get_legend(dummy_plot)
grid.newpage()
grid.draw(legend_dp)

combined_plot <- ggplot(combined_df, aes(x = scenario_label, y = value)) +
  geom_point(aes(y = mean_value, color = dp_s_combo), size = 2) +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value, color = dp_s_combo), width = 0.2) +
  geom_hline(aes(yintercept = observed_value), linetype = "dashed", color = "black") +
  ggh4x::facet_grid2(metric ~ stand, labeller = labeller(stand = c("Allenwiller"="Allenwiller", "Waldi"="Wäldi")), scales = "free_y") +
  scale_color_manual(values = color_map, name = "Dispersal × Selection") +
  labs(x = "Simulation Scenario", y = "Value") +
  theme_bw(base_size = 7) +
  guides(fill = "none", color = "none")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20, unit = "pt"),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 9, margin = margin(b = 5), face = "bold"),
    strip.text.y.left = element_text(angle = 0, size = 9, margin = margin(r = 5)),
    axis.line = element_line(color = 'black'),
    plot.title = element_text(size = 10, hjust = 0.5)
  )


combined_plot <- ggarrange(combined_plot, legend_dp, ncol = 2, widths = c(6,1), align = "hv")

path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/FigureS9_NEMO_MI_SSE.png"), 
    width =7,height=5, res=600, units = "in")
plot(combined_plot)
dev.off()

ggsave(paste0(path, "Manuscript/Figures/FigureS9_NEMO_MI_SSE.pdf"),  plot = combined_plot, 
       width =7,height=5,units = "in")


# -------- Figure X. plot HI best scenarios current future ------
setwd("C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Data_clean/")

# function to process the hybrid index data - proportion per patch (ALL STAGES INCLUDING JUVENILES)
process_HI_data <- function(sim_HI_files, grid_files, stand_names) {
  # all the stages!
  sim_HI_data <- lapply(sim_HI_files, function(file) {
    readRDS(file) %>%
      mutate(scenario2 = paste0(scenario, "_a", s3_age)) %>%
      group_by(stand, stage, scenario2, replicate, pop) %>% # group by stage here
      summarize(
        prop_HI1 = sum(HI == 1) / n(),
        prop_HI0.5 = sum(HI > 0 & HI < 1) / n(),
        prop_HI0 = sum(HI == 0) / n(),
        .groups = "drop"
      ) %>%
      mutate(sum_check = prop_HI1 + prop_HI0.5 + prop_HI0) %>%
      ungroup()
  })
  
  prop_sim_HI_off <- bind_rows(sim_HI_data)
  
  # mean across replicates
  prop_sim_HI_off <- prop_sim_HI_off %>%
    group_by(stand,stage,  pop, scenario2) %>%
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
  
  # merge to grids for plotting
  grid_HI_data <- lapply(seq_along(grid_files), function(i) {
    grid <- read_sf(grid_files[i])
    stand_data <- prop_sim_HI_off %>% filter(stand == stand_names[i])
    grid_HI <- grid %>%
      left_join(stand_data, by = c("patch.ID" = "pop")) %>%
      filter(!is.na(scenario2))
  })
  
  prop_sim_HI_off_grid <- bind_rows(grid_HI_data)
  
  # add param, get the HI category with higher proporition
  prop_sim_HI_off_grid <- prop_sim_HI_off_grid %>%
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
  
  return(prop_sim_HI_off_grid)
}

# process both 2020 and 2120
grid_files <- c("Waldi_grid_4m_5396patches.shp", "Allenwiller_grid_r4m_9844patches.shp")
stand_names <- c("Waldi", "Allenwiller")

sim_HI_files_2020 <- c("NEMO/Quanti_data_waldi_current_HI.RDS", "NEMO/Quanti_data_allen_current_HI_cut.RDS")
grid_2020 <- process_HI_data(sim_HI_files_2020, grid_files, stand_names)

sim_HI_files_2120 <- c("NEMO/Quanti_data_waldi_future_HI.RDS", "NEMO/Quanti_data_allen_future_HI_cut.RDS")
grid_2120 <- process_HI_data(sim_HI_files_2120, grid_files, stand_names)

# get sampling data for the circle plots
samplings <- c("Waldi_sampling_scheme.shp", "Allenwiller_sampling_scheme.shp")
sampling_schemes <- list()
for (i in seq_along(samplings)) {
  sampling_schemes[[i]] <- read_sf(samplings[i]) %>%
    mutate(stand = stand_names[i])
  sampling_schemes[[i]] <- sampling_schemes[[i]][,-c(3:13)] 
}
sampling_all <- bind_rows(sampling_schemes)

# stage column 
stage_lab <- c("1" = "Offspring", "2" = "Juveniles", "3" = "Adults")
grid_2020$stage_lab <- stage_lab[grid_2020$stage]
grid_2120$stage_lab <- stage_lab[grid_2120$stage]

# plot function (one per stage)
plot_HI_split <- function(data_subset, sampling_data, stage_name) {
  stage_data <- data_subset %>% filter(stage_lab == stage_name)
  
  ggplot(stage_data) +
    geom_sf(aes(fill = dominant_HI, alpha = dominance_strength), color = NA) +
    geom_sf(data = sampling_data, shape = 21, color = "grey30", fill = NA,
            size = 0.5, stroke = 0.3, alpha = 0.6) +  
    scale_fill_manual(
      values = c("HI0" = "#FDE725FF", "HI0.5" = "#440154FF", "HI1" = "#25858EFF"),
      name = "HI",
      labels = c("HI = 0", "0 < HI < 1", "HI = 1")
    ) +
    scale_alpha_continuous(range = c(0.2, 1)) +
    theme_minimal(base_size = 10) +
    guides(alpha = "none") +labs(title = stage_name)+
    theme_p +
    theme(axis.title = element_blank(),
          strip.text = element_blank())
}

stages <- c("Adults", "Juveniles", "Offspring")

# Waldi
wal_2020 <- lapply(stages, function(stage)
  plot_HI_split(grid_2020 %>% filter(stand == "Waldi" & Ar == 30 & s == 0 & dp == 20),
                sampling_data = sampling_all %>% filter(stand == "Waldi"),
                stage_name = stage))

wal_2120 <- lapply(stages, function(stage)
  plot_HI_split(grid_2120 %>% filter(stand == "Waldi" & Ar == 30 & s == 0 & dp == 20),
                sampling_data = sampling_all %>% filter(stand == "Waldi"),
                stage_name = stage))
# Allenwiller
all_2020 <- lapply(stages, function(stage)
  plot_HI_split(grid_2020 %>% filter(stand == "Allenwiller" & Ar == 30 & s == 0.5 & dp == 109),
                sampling_data = sampling_all %>% filter(stand == "Allenwiller"),
                stage_name = stage))

all_2120 <- lapply(stages, function(stage)
  plot_HI_split(grid_2120 %>% filter(stand == "Allenwiller" & Ar == 30 & s == 0.5 & dp == 109),
                sampling_data = sampling_all %>% filter(stand == "Allenwiller"),
                stage_name = stage))

# get legend
legend_HI <- cowplot::get_legend(
  wal_2020[[1]] +
    theme(
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 6), 
      legend.key.size = unit(0.3, "cm")
    ) +
    guides(fill = guide_legend(direction = "horizontal"))
)

# remove legends 
all_2020 <- lapply(all_2020, function(p) p + theme(legend.position = "none"))
wal_2020 <- lapply(wal_2020, function(p) p + theme(legend.position = "none"))
all_2120 <- lapply(all_2120, function(p) p + theme(legend.position = "none"))
wal_2120 <- lapply(wal_2120, function(p) p + theme(legend.position = "none"))

# creating rows: add fake labels to the plots for proper alignment
row2a <- annotate_figure(ggarrange(plotlist = all_2020, ncol = 3, align = "h"),
                         left = text_grob("2020", rot = 0, size = 9))

row2w <- annotate_figure(ggarrange(plotlist = wal_2020, ncol = 3, align = "h"),
                         left = text_grob("   ", rot = 0, size = 9))
row2w <- ggarrange(row2w,labels = " ",heights = c(1),# fake label for alignment
                   font.label = list(size = 12, face = "bold"),label.x = 0,  label.y = 1,  hjust = 0, vjust = 1)

row3a <- annotate_figure(ggarrange(plotlist = all_2120, ncol = 3, align = "h"),
                         left = text_grob("2120", rot = 0, size = 9))

row3w <- annotate_figure(ggarrange(plotlist = wal_2120, ncol = 3, align = "h"),
                         left = text_grob("   ", rot = 0, size = 9))
row3w <- ggarrange(row3w,labels = " ",heights = c(1),# fake label for alignment
                   font.label = list(size = 12, face = "bold"),label.x = 0,  label.y = 1,  hjust = 0, vjust = 1)

# merge
row2 <- ggarrange(row2a, row2w, ncol = 2, align = "hv",   labels = c("Allenwiller", "Wäldi"),
                  label.x = 0.3,  label.y = 1,font.label = list(size = 12, face = "bold"))
row3 <- ggarrange(row3a, row3w, ncol = 2, align = "hv")

best_p <- ggarrange(ggarrange(row2,row3, ncol = 1,align = "hv"), legend_HI, ncol = 2,align = "hv", widths = c(5, 1)  
)
best_p

path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/Figure5_best_scenarios.png"), 
    width =7.0866, height=3, units = "in", res=600)
plot(best_p)
dev.off()

ggsave(paste0(path, "Manuscript/Figures/Figure5_best_scenarios.pdf"),  plot = best_p, 
       width =7.0866, height=3, units = "in")




