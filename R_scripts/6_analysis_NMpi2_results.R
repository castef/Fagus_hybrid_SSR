## NMpi2 analysis of results

library(plyr)
library(xlsx)
library(dplyr)
library(gplots)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(lubridate) ## working with data formats
library("tidyverse")
library("graph4lg")
library("raster")
library("tidyr")
library("terra")
library("data.table")

rm(list=ls())

#path = "C:/Users/Camilla/Dropbox (Old (1))/Dropbox/WSL_PhD/Projects/Hybridization/Parentage analysis/NMpi2/Camilla_data/"
path= "C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Data_clean"

setwd("C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Data_clean")
gc()

# custom colors for subspecies
mycols = c("#9E9AC8", "#F0E442")
names(mycols) = c("orientalis", "sylvatica")
colSet = scale_colour_manual(values = mycols)
fillSet = scale_fill_manual(values = mycols)

# ---- convert NMPi output files -----------------
output_files <- list.files(path = path, recursive = TRUE, full.names = TRUE, include.dirs = TRUE, 
                           pattern = "compcoeff.out")

output <- list()

for (i in seq_along(output_files)) {
  file <- output_files[i]
  
  tryCatch({
    data <- fread(file, fill = TRUE)
    # tidy the file and keep only relevant columns
    data <- data[c(3:5),]
    data <- transpose(data[, -1])
    colnames(data) <- c("Parameter", "Estimate", "SE")
    # stand
    stand <- ifelse(grepl("Allenwiller", basename(file)), "Allenwiller", "Waldi")
    data$stand <- stand
    output[[i]] <- data
  })
} 

output <- rbindlist(output, fill =TRUE)

# write results into csv
fwrite(output, "Mating_model/NMpi2_estimates.csv")

# ---- NMpi2 results dispersal curves #####

# Define the function K_jk 
K <- function(d, b, x) {
  # Calculate the scale parameter a
  a <- d * exp(lgamma(2 / b) - lgamma(3 / b))
  
  #  kernel function K_jk
  return(b * exp(-(x^b / a^b))  / (2* pi * a^2 * exp(lgamma(2/b))))
}

## Allenwiller
x <- 1:500  # Range for djk values (distance)
d <- c(119)  # Mean forward seed and pololen dispersal distance
b <- c(0.26282) # shape of seed and  pollen kernel
dispersal_labels <- "Seed"

# Create an empty data frame to store the results
results <- data.frame()

# Calculate K_jk for each combination of delta and b
for (i in 1:length(d)) {
  d_i <- d[i]
  b_i <- b[i]
  y_values <- K(d_i, b_i, x) # calculate the K value
  
  temp <- data.frame(x = x, K = y_values, 
                     d = as.factor(d_i), 
                     b = as.factor(b_i),
                     dispersal_type = dispersal_labels[i])
  results <- rbind(results, temp)
}


all_disp <- ggplot(results, aes(x = x, y = log10(K), color = paste(dispersal_type, "mean dispersal: d =", d, "m, b =", sprintf("%.2f", as.numeric(as.character(b))))))+
  scale_color_manual(values = "#05CDAC")+
  geom_line(size = 1) +
  labs(x = "Euclidean Distance",
       y = "Probability of dispersal (Log10)",
       color = "", 
       title = "Allenwiller") +
  theme(legend.position = c(0.98, 0.98),            # Fine-tune the legend position within the plot
        legend.justification = c(1, 1),   # Anchor legend to top-right
        legend.text = element_text(size = 12),
        panel.background = element_rect(fill = "white", colour = "black"), 
        legend.key = element_rect(fill = "transparent", color = NA)  # Removes border around legend keys
  )
all_disp

## Waldi
x <- 1:500  # Range for djk values (distance)
d <- c(28.95, 92.93680297)  # Mean forward dispersal distance (seed and pollen)
b <- c(0.53259, 0.8504)
dispersal_labels <- c("Seed", "Pollen")

# Create an empty data frame to store the results
results <- data.frame()

# Calculate K_jk for each combination of delta and b
for (i in 1:length(d)) {
  d_i <- d[i]
  b_i <- b[i]
  
  y_values <- K(d_i, b_i, x)
  
  temp <- data.frame(x = x, K = y_values, 
                     d = as.factor(d_i), 
                     b = as.factor(b_i),
                     dispersal_type = dispersal_labels[i])
  
  results <- rbind(results, temp)
}


wal_disp <- ggplot(results, aes(x = x, y = log(K), color = paste(dispersal_type, "mean dispersal: d =", sprintf("%.2f", as.numeric(as.character(d))), "m, b =", sprintf("%.2f", as.numeric(as.character(b))))))+
  scale_color_manual(values = c("#006666", "#05CDAC"))+
  geom_line(size = 1) +
  labs(x = "Euclidean Distance",
       y = "Probability of dispersal (Log10)",
       color = "", 
       title = "Waldi") +
  theme(legend.position = c(0.98, 0.98),            # Fine-tune the legend position within the plot
        legend.justification = c(1, 1),             # Anchor legend to top-right
        legend.text = element_text(size = 12),
        
        panel.background = element_rect(fill = "white", colour = "black"), 
        legend.key = element_rect(fill = "transparent", color = NA) ) # Removes border around legend keys

wal_disp

### combine plots
combined_plot <- plot_grid(all_disp, wal_disp, ncol = 2, rel_widths = c(1, 1))
combined_plot

path_fig = "C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Manuscript/Figures/"
png(filename=paste0(path_fig, "Dispersal_curves_log10.png"), width =23000,height=8000, res=2000)
plot(combined_plot)
dev.off()


### Allewniller and Waldi ML results - selection gradients ####
path = "C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Parentage analysis/NMpi2/Camilla_data/"

sel_grad = read.xlsx(paste0(path, "/NMpi2_overall_results.xlsx"), sheetName = "Selection gradients")
sel_grad$labels = factor(sel_grad$Sel_gradient, levels = c("DBH","Psyl", "Psyl biparental","Competition_coeff"),
                            labels = c("DBH","Species effect (syl)","Species effect assortative (syl)", "Competition"))

## effect of selection gradients on reproductinve success
sel_grad <- sel_grad %>% filter(!(Sel_gradient == "Psyl biparental" & Fem_male == "Female"))

sg_bar = ggplot(subset(sel_grad, Sel_gradient != "Psyl biparental"), 
       aes(x = as.numeric(Effect), y = factor(labels), color = factor(Fem_male))) + 
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = as.numeric(Effect) - SE, xmax = as.numeric(Effect) + SE), 
                 height = 0.1, size = 0.8, position = position_dodge(width = 0.1)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  facet_grid(Site ~ Fem_male) +
  scale_color_manual(values = c("#009CA8", "#006666"))+
  #geom_text(aes(label = paste0("Effect: ", round(as.numeric(Effect), 2), "\nZ: ", round(Z_score, 2)), x = 0), 
           # position = position_dodge(width = 0.5), vjust = 0.5, hjust = -0.1, size = 4) +
  labs(y = "Phenotypic character", x = "Effect on reproductive success (slope coefficients)") +
  xlim(-1, 1) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(size = 16),
    legend.position = "right",
    legend.direction = "vertical",
    strip.background = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    axis.line.x.bottom = element_line(color = "black")
  )

# bar of the Psyl biparental only
sg_bar2 = 
  ggplot(subset(sel_grad, Sel_gradient == "Psyl biparental" & Fem_male == "Male"), 
         aes(x = as.numeric(Effect), y = factor(labels), color = factor(Fem_male))) + 
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = as.numeric(Effect) - SE, xmax = as.numeric(Effect) + SE), 
                 height = 0.1, size = 0.8, position = position_dodge(width = 0.1)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  facet_grid(Site ~ Fem_male) +
  scale_color_manual(values = c("#009CA8", "#006666"))+
  #geom_text(aes(label = paste0("Effect: ", round(as.numeric(Effect), 2), "\nZ: ", round(Z_score, 2)), x = 0), 
  # position = position_dodge(width = 0.5), vjust = 0.5, hjust = -0.1, size = 4) +
  labs(y = "Phenotypic character", x = "Effect on reproductive success (slope coefficients)") +
  xlim(-6, 6) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(size = 16),
    legend.position = "right",
    legend.direction = "vertical",
    strip.background = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    axis.line.x.bottom = element_line(color = "black")
  )

sg_bar2

# combining plots
library(cowplot)
combined_plot = cowplot::plot_grid(
  sg_bar, sg_bar2,
  #plot_grid(lamb_a, lamb_w, ncol = 1), 
  ncol = 2, 
  rel_widths = c(2, 2)  # Adjust width ratio as needed
)
combined_plot

png(filename=paste0(path, "Figures/Sel_grad_all.png"), width =2000,height=1600, res=200)
plot(combined_plot)
dev.off()



########### to correct
# Define the functions for effect of phenotypic character on repr succ (female and male)
lambda_f <- function(z_j, gamma) {
  exp(sum(gamma * z_j))
}


lambda_k_m <- function(z_j, z_k, beta, c) {
  exp(sum(beta * (c * abs(z_j - z_k) + (1 - c) * z_k)))
}

# gamma, beta and c are vectors of length N = number of phenotypic characters
# z_j and z_k are matrices with each row representing a vector of length N = phenotypic characters

# getting tree phenotypic charaters full table
tree.pc = read.xlsx(paste0(path, "/Allenwiller/Allenwiller_NMpi2_data.xlsx"), sheetName = "all_data")

# getting the tree ID of the mother trees used in the candidate parents 
tree.id = read.xlsx(paste0(path, "NMpi2_overall_results.xlsx"), sheetName = "ML_genealogies")

# Filter for non-NA Mo1_ID and Fa1_ID
filtered_tree.id = tree.id[tree.id$Mo1_ID != "NA" & tree.id$Fa1_ID!= "NA", ]

# tree.id to include only rows where Test is "Allenwiller"
filtered_tree.id = filtered_tree.id[filtered_tree.id$Test == "Allenwiller", ]

# gettintg mother data
mothers.pc = merge(filtered_tree.id[, "Mo1_ID", drop = FALSE], 
                    tree.pc[, c("ID", "Z_DBH", "Z_Psyl", "Z_comp_coeff")], 
                    by.x = "Mo1_ID",by.y = "ID")

# getting father data
fathers.pc = merge(filtered_tree.id[, "Fa1_ID", drop = FALSE], 
                    tree.pc[, c("ID", "Z_DBH", "Z_Psyl", "Z_comp_coeff")], 
                    by.x = "Fa1_ID", by.y = "ID")
# convert to numeric
mothers.pc[, -1] = lapply(mothers.pc[, -1], as.numeric)
fathers.pc[, -1] = lapply(fathers.pc[, -1], as.numeric)

# convert to matrices with 3 rows (one for each Z column)
z_j = t(as.matrix(mothers.pc[, -1]))
z_k = t(as.matrix(fathers.pc[, -1]))
# other paraemters
gamma = c(0.56,0.23,-0.2)  # gamma values for the 3 phenotypic characters (DBH, Psyl, C)
beta = c(0.46, -1.36, 0.23)  # beta values for the 3 phenotypic characters (DBH, Psyl, C)
c = c(0, 1, 0)  # c values for the 3 phenotypic characters (DBH, Psyl, C)

# lambda_f and lambda_k_m for each column (tree) in z_j and z_k
lambda_f_values = apply(z_j, 2, lambda_f, gamma = gamma)
lambda_k_m_values = mapply(lambda_k_m, split(z_j, col(z_j)), split(z_k, col(z_k)), MoreArgs = list(beta = beta, c = c))
## function is splitting the matrix z_j and z_k into list of columns and then col generates column indices
## MoreArgs is adding additional arguments beta and c to the function

# dataframe with all together for male and female
lambda_values = data.frame(z_DBH_f = c(z_j[1,]),z_Psyl_f = c(z_j[2,]),z_C_f = c(z_j[3,]),lambda_f = c(lambda_f_values), f_ID = c(mothers.pc[,1]),
                           z_DBH_m = c(z_k[1,]),z_Psyl_m = c(z_k[2,]),z_C_m = c(z_k[3,]), lambda_m = c(lambda_k_m_values), m_ID = c(fathers.pc[,1]))

lambda_long = lambda_values %>%
  pivot_longer(cols = starts_with("z"), names_to = c(".value", "pc_gender"), names_pattern = "^(.*)_(.*)$") %>%
  pivot_longer(cols = c(lambda_f, lambda_m), names_to = "lambda_gender", values_to = "lambda", names_prefix = "lambda_") %>%
  pivot_longer(cols = c(z_DBH, z_Psyl, z_C), names_to = "phenotypic_char", values_to = "pc_value")

# range for the ticks
gender_colors = c("f" = "grey90", "m" = "grey60")
my_theme = theme(legend.title = element_blank(),
                 panel.grid = element_blank(), 
                 panel.border = element_blank(),  
                 axis.line = element_blank(),  
                 axis.ticks = element_blank(),  
                 axis.text.x = element_blank(),  
                 axis.text.y = element_blank()) 

# Create theoretical values for plotting
theoretical_values <- lambda_long %>%
  group_by(phenotypic_char, lambda_gender) %>%
  summarize(
    min_value = min(pc_value),
    max_value = max(pc_value)
  ) %>%
  ungroup() %>%
  mutate(
    pc_value = map2(min_value, max_value, ~seq(.x, .y, length.out = 100))
  ) %>%
  unnest(pc_value) %>%
  mutate(
    lambda_theoretical = case_when(
      lambda_gender == "f" ~ exp(gamma[match(phenotypic_char, c("z_DBH", "z_Psyl", "z_C"))] * pc_value),
      lambda_gender == "m" ~ exp(beta[match(phenotypic_char, c("z_DBH", "z_Psyl", "z_C"))] * pc_value)
    )
  )
lamb_a = ggplot(lambda_long, aes(x = pc_value, y = log(lambda), color = lambda_gender)) +
            geom_point(alpha = 0.5) +
            scale_color_manual(values = gender_colors) + 
            geom_line(data = theoretical_values, aes(x = pc_value, y = log(lambda_theoretical), color = lambda_gender), size = 1) + 
            theme_minimal() + 
            labs(x = "Standardized Phenotypic character",
                 y = "Log-fecundity (Lambda)") +
            my_theme +  ggtitle("Allenwiller")  +
            facet_wrap(~phenotypic_char)+
            geom_hline(yintercept = 0, color = "black") +  
            geom_vline(xintercept = 0, color = "black") +  
            geom_segment(data = data.frame(x = seq(floor(min(lambda_long$pc_value)), ceiling(max(lambda_long$pc_value)), by = 2)), 
               aes(x = x, xend = x, y = -0.1, yend = 0.1), color = "black") +
            geom_segment(data = data.frame(y = seq(floor(min(log(lambda_long$lambda))), ceiling(max(log(lambda_long$lambda))), by = 2)), 
               aes(x = -0.05, xend = 0.05, y = y, yend = y), color = "black") +
            geom_text(data = data.frame(x = seq(floor(min(lambda_long$pc_value)), ceiling(max(lambda_long$pc_value)), by = 2), y = -0.2), 
              aes(x = x, y = y, label = x), vjust = 1, color = "black") +
            geom_text(data = data.frame(x = -0.2, y = seq(floor(min(log(lambda_long$lambda))), ceiling(max(log(lambda_long$lambda))), by = 2)), 
              aes(x = x, y = y, label = y), hjust = 1, color = "black")

lamb_a          
rm(tree.id)
rm(tree.pc)
rm(lambda_long)

# Waldi
tree.pc = read.xlsx(paste0(path, "/Waldi/Waldi_NMpi2_data.xlsx"), sheetName = "wal_data")
# getting the tree ID of the mother trees used in the candidate parents 
tree.id = read.xlsx(paste0(path, "NMpi2_overall_results.xlsx"), sheetName = "ML_genealogies")

# Filter for non-NA Mo1_ID and Fa1_ID
filtered_tree.id = tree.id[tree.id$Mo1_ID != "NA" & tree.id$Fa1_ID!= "NA", ]

# tree.id to include only rows where Test is "Waldi"
filtered_tree.id = filtered_tree.id[filtered_tree.id$Test == "Waldi", ]

# gettintg mother data
mothers.pc = merge(filtered_tree.id[, "Mo1_ID", drop = FALSE], 
                   tree.pc[, c("ID", "Z_DBH", "Z_Psyl", "Z_comp_coeff")], 
                   by.x = "Mo1_ID",by.y = "ID")

# getting father data
fathers.pc = merge(filtered_tree.id[, "Fa1_ID", drop = FALSE], 
                   tree.pc[, c("ID", "Z_DBH", "Z_Psyl", "Z_comp_coeff")], 
                   by.x = "Fa1_ID", by.y = "ID")
# convert to numeric
mothers.pc[, -1] = lapply(mothers.pc[, -1], as.numeric)
fathers.pc[, -1] = lapply(fathers.pc[, -1], as.numeric)

# convert to matrices with 3 rows (one for each Z column)
z_j = t(as.matrix(mothers.pc[, -1]))
z_k = t(as.matrix(fathers.pc[, -1]))
# other paraemters
gamma = c(0.76164,-0.15678  ,-0.29216 )  # gamma values for the 3 phenotypic characters (DBH, Psyl, C)
beta = c(0.72197, -0.24152, -0.36405)  # beta values for the 3 phenotypic characters (DBH, Psyl, C)
c = c(0, 1, 0)  # c values for the 3 phenotypic characters (DBH, Psyl, C)

# lambda_f and lambda_k_m for each column (tree) in z_j and z_k
lambda_f_values = apply(z_j, 2, lambda_f, gamma = gamma)
lambda_k_m_values = mapply(lambda_k_m, split(z_j, col(z_j)), split(z_k, col(z_k)), MoreArgs = list(beta = beta, c = c))
## function is splitting the matrix z_j and z_k into list of columns and then col generates column indices
## MoreArgs is adding additional arguments beta and c to the function

## attach back the lambda values to the original df to find which trees have higher/lower lambda
# dataframe with all together for male and female
lambda_values = data.frame(z_DBH_f = c(z_j[1,]),z_Psyl_f = c(z_j[2,]),z_C_f = c(z_j[3,]),lambda_f = c(lambda_f_values), f_ID = c(mothers.pc[,1]),
                           z_DBH_m = c(z_k[1,]),z_Psyl_m = c(z_k[2,]),z_C_m = c(z_k[3,]), lambda_m = c(lambda_k_m_values), m_ID = c(fathers.pc[,1]))

lambda_long = lambda_values %>%
  pivot_longer(cols = starts_with("z"), names_to = c(".value", "pc_gender"), names_pattern = "^(.*)_(.*)$") %>%
  pivot_longer(cols = c(lambda_f, lambda_m), names_to = "lambda_gender", values_to = "lambda", names_prefix = "lambda_") %>%
  pivot_longer(cols = c(z_DBH, z_Psyl, z_C), names_to = "phenotypic_char", values_to = "pc_value")

# plots for each variable 
# range for the ticks
my_theme = theme(legend.title = element_blank(),
                 panel.grid = element_blank(), 
                 panel.border = element_blank(),  
                 axis.line = element_blank(),  
                 axis.ticks = element_blank(),  
                 axis.text.x = element_blank(),  
                 axis.text.y = element_blank()) 

# Create theoretical values for plotting
theoretical_values <- lambda_long %>%
  group_by(phenotypic_char, lambda_gender) %>%
  summarize(
    min_value = min(pc_value),
    max_value = max(pc_value)
  ) %>%
  ungroup() %>%
  mutate(
    pc_value = map2(min_value, max_value, ~seq(.x, .y, length.out = 100))
  ) %>%
  unnest(pc_value) %>%
  mutate(
    lambda_theoretical = case_when(
      lambda_gender == "f" ~ exp(gamma[match(phenotypic_char, c("z_DBH", "z_Psyl", "z_C"))] * pc_value),
      lambda_gender == "m" ~ exp(beta[match(phenotypic_char, c("z_DBH", "z_Psyl", "z_C"))] * pc_value)
    )
  )

lamb_w = ggplot(lambda_long, aes(x = pc_value, y = log(lambda), color = lambda_gender)) +
            geom_point(alpha = 0.5) +  
  geom_line(data = theoretical_values, aes(x = pc_value, y = log(lambda_theoretical), color = lambda_gender), size = 1) + 
            scale_color_manual(values =gender_colors) + 
            theme_minimal() +  # Clean theme
            labs(x = "Standardized phenotypic charater", y = "Log-fecundity (Lambda)") +
            my_theme + ggtitle("Waldi")  +
  facet_wrap(~phenotypic_char)+
  geom_hline(yintercept = 0, color = "black") +  
  geom_vline(xintercept = 0, color = "black") +  
  geom_segment(data = data.frame(x = seq(floor(min(lambda_long$pc_value)), ceiling(max(lambda_long$pc_value)), by = 2)), 
               aes(x = x, xend = x, y = -0.1, yend = 0.1), color = "black") +
  geom_segment(data = data.frame(y = seq(floor(min(log(lambda_long$lambda))), ceiling(max(log(lambda_long$lambda))), by = 2)), 
               aes(x = -0.05, xend = 0.05, y = y, yend = y), color = "black") +
  geom_text(data = data.frame(x = seq(floor(min(lambda_long$pc_value)), ceiling(max(lambda_long$pc_value)), by = 2), y = -0.2), 
            aes(x = x, y = y, label = x), vjust = 1, color = "black") +
  geom_text(data = data.frame(x = -0.2, y = seq(floor(min(log(lambda_long$lambda))), ceiling(max(log(lambda_long$lambda))), by = 2)), 
            aes(x = x, y = y, label = y), hjust = 1, color = "black")

lamb_w

## combining all the plots 
install.packages("cowplot")
library(cowplot)
combined_plot = cowplot::plot_grid(
  sg_bar, 
  plot_grid(lamb_a, lamb_w, ncol = 1), 
  ncol = 2, 
  rel_widths = c(1, 1)  # Adjust width ratio as needed
)

# Display the combined plot
print(combined_plot)


### Allenwiller and Waldi ML results - genealogies ####
path = "C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Parentage analysis/NMpi2/Camilla_data/"

# load final tests
all_par = read.table(paste0(path, "Allenwiller/Test_3sel_grad/Test49/Allenwiller_NMpi2_input_DBH_Psyl_Psylb_compcoeff.par"), header = T)
all_par$Test = "Allenwiller"
all_par = all_par[, -c(5:7)]

# get the ID from the input NMpi2 dataset 
all_ID  = read.xlsx(paste0(path, "Allenwiller/Allenwiller_NMpi2_data.xlsx"), sheetIndex = 1)
## get the DBH, coordinates and classification from NewHybrids
al_info <- fread("NewHybrids/Allenwiller_indiv_info_class.csv")

# associate NMpi2 code with ID and generation
id_map <- all_ID %>% dplyr::select(ID_number, generation, SampleID = ID)
all_par_ID <- all_par %>%
  mutate(Prog_gen = "Offspring",Mo1_gen = "Adult",Fa1_gen = "Adult")
# add progenies ID
all_par_ID <- all_par_ID %>%left_join(id_map, by = c("Prog" = "ID_number", "Prog_gen" = "generation")) %>%dplyr::rename(Prog_ID = SampleID)
# add mother ID
all_par_ID <- all_par_ID %>%left_join(id_map, by = c("Mo1" = "ID_number", "Mo1_gen" = "generation")) %>%dplyr::rename(Mo1_ID = SampleID)
# add father ID
all_par_ID <- all_par_ID %>%left_join(id_map, by = c("Fa1" = "ID_number", "Fa1_gen" = "generation")) %>%dplyr::rename(Fa1_ID = SampleID)
# add infos from info file
all_par_ID <- all_par_ID %>%
  left_join(al_info %>% dplyr::select(SampleID, Prog_DBH = DBH, Prog_X = x, Prog_Y = y, Prog_NH = NH_class_t80),by = c("Prog_ID" = "SampleID")) %>%
  left_join(al_info %>%  dplyr::select(SampleID, Mo1_DBH = DBH, Mo1_X = x, Mo1_Y = y, Mo1_NH = NH_class_t80),by = c("Mo1_ID" = "SampleID")) %>%
  left_join(al_info %>%  dplyr::select(SampleID, Fa1_DBH = DBH, Fa1_X = x, Fa1_Y = y, Fa1_NH = NH_class_t80),by = c("Fa1_ID" = "SampleID"))

## Waldi
# load final tests
wal_par = read.table(paste0(path, "Waldi/Test_3sel_grad/Test18/Waldi_NMpi2_input_DBH_Psyl_Psylb_compcoeff.par"), header = T)
wal_par$Test = "Waldi"
wal_par = wal_par[, -c(5:7)]

# get the ID from the input NMpi2 dataset 
wal_ID  = read.xlsx(paste0(path, "Waldi/Waldi_NMpi2_data.xlsx"), sheetIndex = 1)
## get the DBH, coordinates and classification from NewHybrids
wal_info <- fread("NewHybrids/Waldi_indiv_info_class.csv")

# associate NMpi2 code with ID and generation
id_map <- wal_ID %>% dplyr::select(ID_number, generation, SampleID = ID)
wal_par_ID <- wal_par %>%
  mutate(Prog_gen = "Offspring",Mo1_gen = "Adult",Fa1_gen = "Adult")
# add progenies ID
wal_par_ID <- wal_par_ID %>%left_join(id_map, by = c("Prog" = "ID_number", "Prog_gen" = "generation")) %>%dplyr::rename(Prog_ID = SampleID)
# add mother ID
wal_par_ID <- wal_par_ID %>%left_join(id_map, by = c("Mo1" = "ID_number", "Mo1_gen" = "generation")) %>%dplyr::rename(Mo1_ID = SampleID)
# add father ID
wal_par_ID <- wal_par_ID %>%left_join(id_map, by = c("Fa1" = "ID_number", "Fa1_gen" = "generation")) %>%dplyr::rename(Fa1_ID = SampleID)
# add infos from info file
wal_par_ID <- wal_par_ID %>%
  left_join(wal_info %>% dplyr::select(SampleID, Prog_DBH = DBH, Prog_X = x, Prog_Y = y, Prog_NH = NH_class_t80_corrected),by = c("Prog_ID" = "SampleID")) %>%
  left_join(wal_info %>%  dplyr::select(SampleID, Mo1_DBH = DBH, Mo1_X = x, Mo1_Y = y, Mo1_NH = NH_class_t80_corrected),by = c("Mo1_ID" = "SampleID")) %>%
  left_join(wal_info %>%  dplyr::select(SampleID, Fa1_DBH = DBH, Fa1_X = x, Fa1_Y = y, Fa1_NH = NH_class_t80_corrected),by = c("Fa1_ID" = "SampleID"))

par_results = rbind(all_par_ID, wal_par_ID)

#write.csv(par_results, "Mating_model/Parentage_results.csv")

# summarizing results
par_results = read.csv("Mating_model/Parentage_results.csv")
str(par_results)
par_results <- par_results %>%mutate(across(where(is.character), ~ gsub("orientalis", "hohenackeriana", .)))
par_results <- par_results %>%mutate(across(where(is.character), ~ gsub("Waldi", "Wäldi", .)))

## summary
par_sum = par_results %>% group_by(Test) %>% 
  summarise(not_passed_gen = sum(Pr1 < 0.8),
            # genealogies with probability lower than 0.8
            passed_gen = sum(Pr1 > 0.8), 
            # genealogies with probability higher than 0.8
            both_assigned = sum(Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1), 
            # offspring with both parents assigned
            single_parent = sum((Pr1 > 0.8 & Mo1 == -1 & Fa1 != -1) | (Pr1 > 0.8 & Mo1 != -1 & Fa1 == -1), na.rm = TRUE),
            # offspring with one parent assigned
            no_parents = sum(Pr1 > 0.8 & Mo1 == -1 & Fa1 == -1, na.rm = TRUE ), 
            # offspring with no parents assigned
            N_mo1_syl = length(unique(Mo1[Mo1_NH > 0.8 & !is.na(Mo1_NH)])),
            # number of sylvatica mothers
            N_mo1_ori = length(unique(Mo1[Mo1_NH < 0.1 & !is.na(Mo1_NH)])),
            # number of orientalis mothers
            N_fa1_syl = length(unique(Fa1[Fa1_NH > 0.8 & !is.na(Fa1_NH)])),
            # number of sylvatica fathers
            N_fa1_ori = length(unique(Fa1[Fa1_NH < 0.1 & !is.na(Fa1_NH)])),
            # number of orientalis fathers
            N_syl_syl = sum(Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1 & Mo1_NH >= 0.8 & Fa1_NH >= 0.8, na.rm = TRUE),
            # number of offspring with both sylvatica parents
            N_ori_ori = sum(Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1 & Mo1_NH <= 0.1 & Fa1_NH <= 0.1, na.rm = TRUE),
            # number of offspring with both orientalis parents
            N_syl_ori = sum(Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1 & Mo1_NH >= 0.8 & Fa1_NH <= 0.1, na.rm = TRUE),
            # number of offspring with sylvatica mother and orientalis father
            N_ori_syl = sum(Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1 & Mo1_NH <= 0.1 & Fa1_NH >= 0.8, na.rm = TRUE)
            # number of offspring with orientalis mother and sylvatica mother
  ) 

par_sum = base::t(par_sum)

# take only progenies with high probability
filtered <- par_results %>% filter(Pr1 > 0.8)
# remove unassigned
filtered <- par_results %>% filter(Prog_NH != "unassigned")

# classify offspring into no parents, 2 parents and both parents assigned
both_parents <- filtered %>%filter(!is.na(Mo1_ID) & !is.na(Fa1_ID))
one_parent <- filtered %>% filter(xor(is.na(Mo1_ID), is.na(Fa1_ID)))
no_parents <- filtered %>%filter(is.na(Mo1_ID) & is.na(Fa1_ID))
table(both_parents$Prog_NH, both_parents$Test)
table(one_parent$Prog_NH,one_parent$Test)
table(no_parents$Prog_NH, no_parents$Test)

# offspring genotype for each type of genealogy
ColourVector3 = c("Pure sylvatica" = "#FDE725FF", "Pure hohenackeriana" = "#482173FF", "F1" = "#25858EFF", "F2" = "#85D54AFF", "BC sylvatica" = "orange","BC hohenackeriana" ="#9C99C8", "unassigned" = "gray8")

combined <- bind_rows(
  both_parents %>% mutate(group = "Both parents"),
  one_parent %>% mutate(group = "One parent"),
  no_parents %>% mutate(group = "No parents")
)
combined <- combined[!(combined$Prog_NH == "unassigned"),]
combined$group <- factor(combined$group, levels = c("Both parents", "One parent", "No parents"))
combined$Prog_NH <- factor(combined$Prog_NH,levels = c("Pure sylvatica","Pure hohenackeriana","F1","F2","BC sylvatica", "BC hohenackeriana"))

overview <- ggplot(combined, aes(x = group, fill = Prog_NH)) +
  geom_bar(position = position_dodge2(width = 0.9, preserve = "single"), col = "black") +
  facet_wrap(~Test, scale = "free") +
  labs(title = "Offspring class per genealogy type",
       x = "", y = "Count", fill = "NH Class") +
  scale_fill_manual(values = ColourVector3, breaks = names(ColourVector3), drop = F)+
  theme_bw(base_size = 7) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.text = element_text(size = 10),
    legend.spacing.x = unit(0.5, "cm"),
    legend.key.size = unit(0.6, "cm"),
    legend.key.height = unit(0.6, "cm"),
    legend.box.just = "center"
  ) +
  guides(fill = guide_legend(nrow = 1))
overview

# proportion per type
combined$Prog_NH <- factor(combined$Prog_NH, levels = names(ColourVector3))

overview2 <- ggplot(combined, aes(x = Prog_NH, fill = group)) +
  geom_bar(position = "fill", col = "black", width = 0.5, linewidth = 0.2) +
  scale_fill_manual(values = c("white", "grey", "black")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "",
       x = "Offspring genotype", y = "Percentage", fill = "N. of assigned\nparents") +
  theme_bw(base_size = 7) + facet_wrap(~Test)+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid.minor = element_blank(),
        legend.position = "inside", 
        legend.position.inside = c(0.90, 0.80),
        legend.key = element_rect(colour = "white"),
        legend.text = element_text(size = 6),           # Text size in legend
        legend.title = element_text(size = 7),          # Legend title size
        legend.key.size = unit(0.3, "cm"),              # Size of legend keys (squares)
        legend.key.width = unit(0.3, "cm"),             # Width of legend keys
        legend.key.height = unit(0.3, "cm"),            # Height of legend keys
        legend.spacing.y = unit(0.1, "cm"),             # Vertical spacing between legend items
        legend.margin = margin(2, 2, 2, 2),             # Margin around entire legend (smaller)
        legend.box.spacing = unit(0.1, "cm"), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        plot.title = element_text(hjust = 0.5))
overview2

# genotype of the "single" parents (all offspring)
ggplot(one_parent, aes(x = Mo1_NH, fill = Mo1_NH)) +
  geom_bar(position = position_dodge2(width = 0.9, preserve = "single")) +facet_wrap(~Test, scale = "free") +
  scale_x_discrete(limits = names(ColourVector3)) +
  labs(title = "Single parent's class",
       x = "", y = "Count", fill = "Parent Type") +
  theme_bw(base_size = 7) +  scale_fill_manual(values = ColourVector3)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# genotpye of single parent of each offspring
one_parent$Prog_NH <- factor(one_parent$Prog_NH, levels = names(ColourVector3))

single_F1_parents <- ggplot(one_parent, aes(x = Prog_NH, fill = Mo1_NH)) +
  geom_bar(position = "fill", col = "black", width = 0.5, linewidth = 0.2) +
  scale_fill_manual(values = ColourVector3) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(#title = "Genotype of single parent of F1 offsprng",
       x = "Offspring genotype", y = "Percentage", fill = "Parent NH Class") +
  theme_bw(base_size = 7)+ facet_wrap(~Test)+
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.key = element_rect(colour = "white"),
        #legend.key.spacing.y = unit(0.4, 'cm'), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        plot.title = element_text(hjust = 0.5))

single_F1_parents

# spatial relationship of F1 with their parents
f1_both <- both_parents %>% filter(Prog_NH == "F1")
f1_one <- one_parent %>%filter(Prog_NH == "F1")

# euclidean distance between offspring and hohenackeriana and sylavtica parent
f1_both <- both_parents %>%
  filter(Prog_NH == "F1") %>%
  mutate(
    hohenackeriana_X = ifelse(Mo1_NH == "Pure hohenackeriana", Mo1_X,
                          ifelse(Fa1_NH == "Pure hohenackeriana", Fa1_X, NA)),
    hohenackeriana_Y = ifelse(Mo1_NH == "Pure hohenackeriana", Mo1_Y,
                          ifelse(Fa1_NH == "Pure hohenackeriana", Fa1_Y, NA)),
    
    sylvatica_X = ifelse(Mo1_NH == "Pure sylvatica", Mo1_X,
                         ifelse(Fa1_NH == "Pure sylvatica", Fa1_X, NA)),
    sylvatica_Y = ifelse(Mo1_NH == "Pure sylvatica", Mo1_Y,
                         ifelse(Fa1_NH == "Pure sylvatica", Fa1_Y, NA)),
    
    dist_to_hohenackeriana = sqrt((Prog_X - hohenackeriana_X)^2 + (Prog_Y - hohenackeriana_Y)^2),
    dist_to_sylvatica  = sqrt((Prog_X - sylvatica_X)^2 + (Prog_Y - sylvatica_Y)^2)
  )

dist_compare <- f1_both %>%
  dplyr::select(Test, dist_to_hohenackeriana, dist_to_sylvatica) %>%
  pivot_longer(cols = starts_with("dist_to_"), names_to = "Parent_NH", values_to = "Distance") %>%
  mutate(Parent_Type = recode(Parent_NH,
                              dist_to_hohenackeriana = "Pure hohenackeriana",
                              dist_to_sylvatica = "Pure sylvatica"))

dist_compare$Parent_Type <- factor(dist_compare$Parent_Type, levels = c("Pure sylvatica", "Pure hohenackeriana" ))

## check mean and max distnace of the two parents genotypes: 
dist_compare %>%
  group_by(Test, Parent_NH) %>%
  summarise(
    mean_distance = mean(Distance, na.rm = TRUE),
    max_distance = max(Distance, na.rm = TRUE),
    .groups = "drop"
  )

F1_parents_distance_both <- ggplot(dist_compare, aes(x = Parent_Type, y = Distance, fill = Parent_Type, colour = Parent_Type)) +
  geom_boxplot(alpha = 0.8, col = "black", linewidth = 0.3,outlier.size = 0.5) + 
  guides(fill = "none")+
  #geom_point()+
  scale_fill_manual(values = ColourVector3)+  
  scale_color_manual(values = ColourVector3)+  
  facet_wrap(~Test) +
  labs( #title = "Distance parent - F1 offspring for both assigned parents", 
        x = "Parent genotype", y = "Distance (m) from F1 offspring") +
  theme_bw(base_size = 7)+
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none", 
        legend.key = element_rect(colour = "white"),
        #legend.key.spacing.y = unit(0.4, 'cm'), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        plot.title = element_text(hjust = 0.5))

F1_parents_distance_both

# distance of the single parent from the offspring
f1_one <- f1_one %>%
  mutate(
    Parent_NH = ifelse(!is.na(Mo1_ID), Mo1_NH, Fa1_NH),
    Parent_Type = ifelse(!is.na(Mo1_ID), "Mother", "Father")
  )
f1_one <- f1_one %>%
  mutate( dist_to_parent = sqrt((Prog_X - Mo1_X)^2 + (Prog_Y - Mo1_Y)^2)  )

f1_one$Parent_NH <- factor(f1_one$Parent_NH, levels = names(ColourVector3))

## check mean and max distnace of the two parents genotypes:  
f1_one %>%
  group_by(Test, Parent_NH) %>%
  summarise(
    mean_distance = mean(dist_to_parent, na.rm = TRUE),
    max_distance = max(dist_to_parent, na.rm = TRUE),
    .groups = "drop"
  )


F1_parents_distance_one <- ggplot(f1_one, aes(x = Parent_NH, y = dist_to_parent, fill = Parent_NH, col = Parent_NH)) +
  geom_boxplot(alpha = 0.8, col = "black", linewidth = 0.3,outlier.size = 0.5) + 
  guides(fill = "none")+
  #geom_point()+
  scale_fill_manual(values = ColourVector3)+  
  scale_color_manual(values = ColourVector3)+  
  
  facet_wrap(~Test, scale = "free_y") +
  labs(#title = "Distance parent - F1 offspring for one assigned parent", 
       x = "Parent genotype", y = "Distance (m) from F1 offspring") +
  theme_bw(base_size = 7)+
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.key = element_rect(colour = "white"),
        #legend.key.spacing.y = unit(0.4, 'cm'), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        plot.title = element_text(hjust = 0.5))

F1_parents_distance_one

## combine plots
nmpi_genealogies2 <- ggpubr::ggarrange(overview2,single_F1_parents,
                                       F1_parents_distance_both,F1_parents_distance_one,
                                       heights = c(1, 1),common.legend = F,
                                       labels = c("A","B",  "C", "D"),font.label=list(size=8),
                                       nrow = 2, ncol = 2,align = "h") 

nmpi_genealogies2
# get the legend from the NH plot to get all the categories
nmpi_genealogies3 <- ggpubr::ggarrange(nmpi_genealogies2,legend_with_unassigned,
                                 nrow = 1, ncol = 2,widths = c(1, 0.02), align = "v") 
nmpi_genealogies3

path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/Figure_genealogies.png"), 
    width =6.2, height=4.8, units = "in", res=600)
plot(nmpi_genealogies3)
dev.off()

# mapping genealogies ##
filtered$Assignment <- case_when(
  !is.na(filtered$Mo1_ID) & !is.na(filtered$Fa1_ID) ~ "Both parents",
  xor(is.na(filtered$Mo1_ID), is.na(filtered$Fa1_ID)) ~ "One parent",
  is.na(filtered$Mo1_ID) & is.na(filtered$Fa1_ID) ~ "No parents"
)
f1_filtered <- filtered %>% filter(Prog_NH == "F1")

filtered$Assignment <- factor(filtered$Assignment, levels = c("Both parents", "One parent", "No parents"))
filtered$Prog_NH <- factor(filtered$Prog_NH,levels = c("Pure sylvatica","Pure hohenackeriana","F1","F2","BC sylvatica", "BC hohenackeriana"))
f1_filtered$Assignment <- factor(f1_filtered$Assignment, levels = c("Both parents", "One parent", "No parents"))

# split by stand to get constant coordinates
waldi_map <- ggplot() +
  geom_segment(data = f1_filtered %>% filter(Test == "Wäldi" & !is.na(Mo1_X)), 
               aes(x = Prog_X, y = Prog_Y, xend = Mo1_X, yend = Mo1_Y, color = Mo1_NH),
               arrow = arrow(length = unit(0.15, "cm"))) +
  geom_segment(data = f1_filtered %>% filter(Test == "Wäldi" &!is.na(Fa1_X)), 
               aes(x = Prog_X, y = Prog_Y, xend = Fa1_X, yend = Fa1_Y, color = Fa1_NH),
               arrow = arrow(length = unit(0.15, "cm"))) +
  
  geom_jitter(data = subset(filtered, Test == "Wäldi" ), 
              aes(x = Prog_X, y = Prog_Y, color = Prog_NH), 
              width = 0.8, height = 0.8, size = 2, alpha = 0.5) +
  
  facet_grid(Assignment ~ Test) +
  scale_color_manual(values = ColourVector3) +
 # labs(title = "Offspring spatial location and parent links by Test and assignment", x = "X", y = "Y", color = "NH Class") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.ticks = element_blank(),
        strip.text = element_text(size = 13),
        axis.text = element_blank(), axis.title = element_blank(),
        strip.background = element_rect(colour = "white", fill = "white"), 
        plot.title = element_text(hjust = 0.5))
waldi_map

# split by stand to get constant coordinates
allen_map <- ggplot() +
  geom_segment(data = f1_filtered %>% filter(Test == "Allenwiller" & !is.na(Mo1_X)), 
               aes(x = Prog_X, y = Prog_Y, xend = Mo1_X, yend = Mo1_Y, color = Mo1_NH),
               arrow = arrow(length = unit(0.15, "cm"))) +
  geom_segment(data = f1_filtered %>% filter(Test == "Allenwiller" &!is.na(Fa1_X)), 
               aes(x = Prog_X, y = Prog_Y, xend = Fa1_X, yend = Fa1_Y, color = Fa1_NH),
               arrow = arrow(length = unit(0.15, "cm"))) +
  
  geom_jitter(data = subset(filtered, Test == "Allenwiller" ), 
              aes(x = Prog_X, y = Prog_Y, color = Prog_NH), 
              width = 0.8, height = 0.8, size = 2, alpha = 0.5) +
  
  facet_grid(Assignment ~ Test) +
  scale_color_manual(values = ColourVector3) +
  labs( x = "", y = "", color = "") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        strip.text.y = element_blank(),
        strip.text = element_text(size=13),
        axis.text = element_blank(), 
        axis.title = element_blank(),
        strip.background = element_rect(colour = "white", fill = "white"), 
        plot.title = element_text(hjust = 0.5))
allen_map

map_genealogies <- ggpubr::ggarrange(allen_map, waldi_map, 
                                     nrow = 1, ncol = 2,align = "h") 
  
map_genealogies
map_genealogies2 <- ggpubr::ggarrange(map_genealogies,
                                      genotype_legend,
                                      nrow = 1, ncol = 2,widths = c(1, 0.2), align = "v") 
map_genealogies2

path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/"
png(filename=paste0(path, "Manuscript/Figures/Figure_map_genealogies.png"), width =12000,height=12000, res=1200)
plot(map_genealogies2)
dev.off()


#############################old########################

# N offspring for each parent divided by subspecies and adding as Sheet to result file
Noff_mothers = par_results %>% mutate(Mo1_sp = case_when(Mo1 != -1 & Mo1_Psyl > 0.8 ~ "Sylvatica",
                                                 Mo1 != -1 & Mo1_Psyl < 0.1 ~ "Orientalis", 
                                                 Mo1 == -1 ~ NA))%>% 
  group_by(Test,Mo1_ID, Mo1_X, Mo1_Y, Mo1_DBH, Mo1_sp ) %>% 
  summarise(N_off = n()) %>% 
  ungroup()

Noff_fathers = par_results %>% mutate(Fa1_sp = case_when(Fa1 != -1 & Fa1_Psyl > 0.8 ~ "Sylvatica",
                                                             Fa1 != -1 & Fa1_Psyl < 0.1 ~ "Orientalis", 
                                                             Fa1 == -1 ~ NA))%>% 
  group_by(Test,Fa1_ID, Fa1_X, Fa1_Y, Fa1_DBH, Fa1_sp ) %>% 
  summarise(N_off = n()) %>% 
  ungroup()

write.xlsx(Noff_mothers, "NMpi2_overall_results.xlsx", sheetName="Mothers_Noffspring", append=TRUE)
write.xlsx(Noff_fathers, "NMpi2_overall_results.xlsx", sheetName="Fathers_Noffspring", append=TRUE)

# summary statistics mothers 
Noff_mothers = read.xlsx("NMpi2_overall_results.xlsx", sheetName = "Mothers_Noffspring")

Noff_mothers$N_off = as.numeric(as.character(Noff_mothers$N_off))
hist(Noff_mothers[which(Noff_mothers$Test == "Allenwiller" & Noff_mothers$Mo1_sp == "Sylvatica"),]$N_off)
hist(Noff_mothers[which(Noff_mothers$Test == "Allenwiller" & Noff_mothers$Mo1_sp == "Orientalis"),]$N_off)
hist(Noff_mothers[which(Noff_mothers$Test == "Waldi" & Noff_mothers$Mo1_sp == "Sylvatica"),]$N_off)
hist(Noff_mothers[which(Noff_mothers$Test == "Waldi" & Noff_mothers$Mo1_sp == "Orientalis"),]$N_off)

Noff_mothers_sum = Noff_mothers %>% group_by(Test, Mo1_sp) %>%
                      summarise(N = n(), 
                                Mean_Noff = mean(N_off),
                                SD_Noff = sd(N_off))

# test of the means between N offspring of mother sylvatica and mother orientalis
wilcox.test(Noff_mothers[which(Noff_mothers$Mo1_sp == "Sylvatica"),]$N_off, 
       Noff_mothers[which(Noff_mothers$Mo1_sp == "Orientalis"),]$N_off)

wilcox.test(Noff_mothers[which(Noff_mothers$Test == "Waldi" & Noff_mothers$Mo1_sp == "Sylvatica"),]$N_off,
       Noff_mothers[which(Noff_mothers$Test == "Waldi" & Noff_mothers$Mo1_sp == "Orientalis"),]$N_off)

wilcox.test(Noff_mothers[which(Noff_mothers$Test == "Allenwiller" & Noff_mothers$Mo1_sp == "Sylvatica"),]$N_off, 
       Noff_mothers[which(Noff_mothers$Test == "Allenwiller" & Noff_mothers$Mo1_sp == "Orientalis"),]$N_off)

# summary statistics fathers 
Noff_fathers = read.xlsx("NMpi2_overall_results.xlsx", sheetName = "Fathers_Noffspring")

Noff_fathers$N_off = as.numeric(as.character(Noff_fathers$N_off))
hist(Noff_fathers[which(Noff_fathers$Test == "Allenwiller" & Noff_fathers$Fa1_sp == "Sylvatica"),]$N_off)
hist(Noff_fathers[which(Noff_fathers$Test == "Allenwiller" & Noff_fathers$Fa1_sp == "Orientalis"),]$N_off)
hist(Noff_fathers[which(Noff_fathers$Test == "Waldi" & Noff_fathers$Fa1_sp == "Sylvatica"),]$N_off)
hist(Noff_fathers[which(Noff_fathers$Test == "Waldi" & Noff_fathers$Fa1_sp == "Orientalis"),]$N_off)

Noff_fathers_sum = Noff_fathers %>% group_by(Test, Fa1_sp) %>%
  summarise(N = n(), 
            Mean_Noff = mean(N_off),
            SD_Noff = sd(N_off))

# test of the means between N offspring of mother sylvatica and mother orientalis
wilcox.test(Noff_fathers[which(Noff_fathers$Fa1_sp == "Sylvatica"),]$N_off, 
            Noff_fathers[which(Noff_fathers$Fa1_sp == "Orientalis"),]$N_off)

wilcox.test(Noff_fathers[which(Noff_fathers$Test == "Waldi" & Noff_fathers$Fa1_sp == "Sylvatica"),]$N_off,
            Noff_fathers[which(Noff_fathers$Test == "Waldi" & Noff_fathers$Fa1_sp == "Orientalis"),]$N_off)

wilcox.test(Noff_fathers[which(Noff_fathers$Test == "Allenwiller" & Noff_fathers$Fa1_sp == "Sylvatica"),]$N_off, 
            Noff_fathers[which(Noff_fathers$Test == "Allenwiller" & Noff_fathers$Fa1_sp == "Orientalis"),]$N_off)

# plotting Mothers Allenwiller and Waldi
Noff_mothers[4:7] = sapply(Noff_mothers[4:7], as.numeric)
Noff_fathers[4:7] = sapply(Noff_fathers[4:7], as.numeric)

hist(Noff_mothers$N_off)
size_scale = scale_size_continuous(limits=c(1,80),breaks=c(5,20,40,60,80))

ggplot() + theme_bw() +
      geom_point(subset(Noff_mothers, Test != "Allenwiller_ori_rem"), mapping = aes(x = Mo1_X, y = Mo1_Y, size = N_off, col = Mo1_sp)) + 
     #scale_color_manual(values = c("Sylvatica" = "#F0E442", "Orientalis" = "#9E9AC8", "NA"= "grey"))+size_scale+
    #  geom_text(data = subset(Noff_mothers, Test != "Waldi"), aes(x = X, y = Y,label = N_off), vjust = 1.5, color = "black", size = 3) +
      facet_wrap(~Test, scales = "free") + labs(title = "NMpi2 genealogies - Mother trees. Size = N offspring") +  
      theme(legend.position="none")

# plotting Fathers Allenwiller and Waldi
which.max(Noff_fathers$N_off)
size_scale = scale_size_continuous(limits=c(1,30),breaks=c(1,5,10,15,20))

ggplot() + theme_bw() + 
      geom_point(subset(Noff_fathers, Test != "Allenwiller_ori_rem"), mapping = aes(x = Fa1_X, y = Fa1_Y, size = N_off, col = Fa1_sp)) + 
      #scale_color_manual(values = c("Sylvatica" = "#F0E442", "Orientalis" = "#9E9AC8")) +size_scale+
   #   geom_text(data = subset(Noff_fathers, Test != "Waldi"), aes(x = X, y = Y,label = N_off), vjust = 1.5, color = "black", size = 3) +
  facet_wrap(~Test, scales = "free") + labs(title = "NMpi2 genealogies - Father trees. Size = N offspring") 
      theme(legend.position="none")

 ### comparing results of NH with genealogies by NMpi2 ####
NH = read.csv("C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/NewHybrids/Results/Allenwiller/Allenwiller_coord_circumf_class80ppt.csv")
NMpi = read.xlsx("C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Parentage analysis/NMpi2/Camilla_data/NMpi2_overall_results.xlsx", sheetIndex = 4)

# attaching coordinates and NH genotype of the offspring to NMpi2 file
colnames(NMpi)[9] = "IndivName"
NMpi = left_join(NMpi, NH[, c(2,7,8, 15, 16)], by = "IndivName")

#### checking the genealogies of each class

## TOT Pure1 in NH and  NMpi (Pr > 0.8 and Pr < 0.8)
nrow(subset(NH, max_class == "Pure1"))
nrow(subset(NMpi,Test == "Allenwiller"&  max_class == "Pure1"& Pr1 > 0.8))
nrow(subset(NMpi,Test == "Allenwiller"&  max_class == "Pure1"))
# classified Pure1 with no parents
Pure1_0par = inner_join(subset(NH,max_class == "Pure1"), subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 & Mo1 == -1 & Fa1 == -1 ), by = "IndivName")
nrow(Pure1_0par)
# classified Pure1 with 1 parents
Pure1_1par = inner_join(subset(NH,max_class == "Pure1"), subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 & Mo1 == -1 & Fa1 != -1 | Test == "Allenwiller" & Pr1 > 0.8 & Mo1 != -1 & Fa1 == -1), by = "IndivName")
nrow(Pure1_1par)
# classified Pure1 with 2 parents
Pure1_2par = inner_join(subset(NH,max_class == "Pure1"), subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1), by = "IndivName")
nrow(Pure1_2par)
# plotting 
plot(vect(subset(NMpi, Test == "Allenwiller"), geom=c("X", "Y"), crs="+init=epsg:3035"))
plot(vect(subset(NMpi, Test == "Allenwiller"& max_class == "Pure1"),geom=c("X", "Y"), crs="+init=epsg:3035"), add = T, col = "blue") # all Pure1
plot(vect(Pure1_0par, geom=c("X.x", "Y.x"), crs="+init=epsg:3035"), add = T, col = "red") # Pure1 with 0 parents assigned
plot(vect(Pure1_1par, geom=c("X.x", "Y.x"), crs="+init=epsg:3035"), add = T, col = "orange") # Pure1 with 1 parents assigned
plot(vect(Pure1_2par, geom=c("X.x", "Y.x"), crs="+init=epsg:3035"), add = T, col = "green") # Pure1 with 2 parents assigned
plot(vect(subset(NMpi,Test == "Allenwiller"& max_class == "Pure1"& Pr1 < 0.8),geom=c("X", "Y"), crs="+init=epsg:3035"), add = T, col = "grey") # pure1 with parentage with low probability"
legend("right", col=c("blue", "red", "orange", "green", "gray"), 
       legend=c("all Pure1 by NH", "Pure1 with 0 parents assigned", "Pure1 with 1 parents assigned", "Pure1 with 2 parents assigned", "Pure1 with parentage with low probability"), 
       lty = 1, lwd = 2)

## TOT Pure2  in NH and  NMpi (Pr > 0.8 and Pr < 0.8)
nrow(subset(NH, max_class == "Pure2"))
nrow(subset(NMpi,Test == "Allenwiller"&  max_class == "Pure2"& Pr1 > 0.8))
nrow(subset(NMpi,Test == "Allenwiller"&  max_class == "Pure2"& Pr1 < 0.8))
# classified Pure2 with no parents
Pure2_0par = inner_join(subset(NH,max_class == "Pure2"), subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 & Mo1 == -1 & Fa1 == -1 ), by = "IndivName")
nrow(Pure2_0par)
# classified Pure2 with 1 parents
Pure1_1par = inner_join(subset(NH,max_class == "Pure2"), subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 & Mo1 == -1 & Fa1 != -1 | Test == "Allenwiller" & Pr1 > 0.8 & Mo1 != -1 & Fa1 == -1), by = "IndivName")
nrow(Pure1_1par)
# classified Pure2 with 2 parents
Pure1_2par = inner_join(subset(NH,max_class == "Pure2"), subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1), by = "IndivName")
nrow(Pure1_2par)
# plotting 
plot(vect(subset(NMpi, Test == "Allenwiller"), geom=c("X", "Y"), crs="+init=epsg:3035"))
plot(vect(subset(NMpi, Test == "Allenwiller"& max_class == "Pure1"),geom=c("X", "Y"), crs="+init=epsg:3035"), add = T, col = "blue") # all Pure2
plot(vect(Pure1_0par, geom=c("X.x", "Y.x"), crs="+init=epsg:3035"), add = T, col = "red") # Pure2 with 0 parents assigned
plot(vect(Pure1_1par, geom=c("X.x", "Y.x"), crs="+init=epsg:3035"), add = T, col = "orange") # Pure2 with 1 parents assigned
plot(vect(Pure1_2par, geom=c("X.x", "Y.x"), crs="+init=epsg:3035"), add = T, col = "green") # Pure2 with 2 parents assigned
plot(vect(subset(NMpi,Test == "Allenwiller"& max_class == "Pure1"& Pr1 < 0.8),geom=c("X", "Y"), crs="+init=epsg:3035"), add = T, col = "grey") # Pure2 with parentage with low probability"
legend("right", col=c("blue", "red", "orange", "green", "gray"), 
       legend=c("all Pure2 by NH", "Pure1 with 0 parents assigned", "Pure2 with 1 parents assigned", "Pure2 with 2 parents assigned", "Pure2 with parentage with low probability"), 
       lty = 1, lwd = 2)

## TOT F1 in NH and  NMpi (Pr > 0.8 and Pr < 0.8)
nrow(subset(NH, max_class == "F1"))
nrow(subset(NMpi,Test == "Allenwiller"&  max_class == "F1"& Pr1 > 0.8))
nrow(subset(NMpi,Test == "Allenwiller"&  max_class == "F1"& Pr1 < 0.8))
# classified F1 with no parents
F1_0par = inner_join(subset(NH,max_class == "F1"), subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 & Mo1 == -1 & Fa1 == -1 ), by = "IndivName")
nrow(F1_0par)
# classified F1 with 1 parents
F1_1par = inner_join(subset(NH,max_class == "F1"), subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 & Mo1 == -1 & Fa1 != -1 | Test == "Allenwiller" & Pr1 > 0.8 & Mo1 != -1 & Fa1 == -1), by = "IndivName")
nrow(F1_1par)
# classified F1 with 2 parents
F1_2par = inner_join(subset(NH,max_class == "F1"), subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1), by = "IndivName")
nrow(F1_2par)
# plotting 
plot(vect(subset(NMpi, Test == "Allenwiller"), geom=c("X", "Y"), crs="+init=epsg:3035"))
plot(vect(subset(NMpi, Test == "Allenwiller"& max_class == "F1"),geom=c("X", "Y"), crs="+init=epsg:3035"), add = T, col = "blue") # all F1
plot(vect(F1_0par, geom=c("X.x", "Y.x"), crs="+init=epsg:3035"), add = T, col = "red") # F1 with 0 parents assigned
plot(vect(F1_1par, geom=c("X.x", "Y.x"), crs="+init=epsg:3035"), add = T, col = "orange") # F1 with 1 parents assigned
plot(vect(F1_2par, geom=c("X.x", "Y.x"), crs="+init=epsg:3035"), add = T, col = "green") # F1 with 2 parents assigned
plot(vect(subset(NMpi,Test == "Allenwiller"& max_class == "F1"& Pr1 < 0.8),geom=c("X", "Y"), crs="+init=epsg:3035"), add = T, col = "grey") #  F1 with low probability of genealogy
legend("right", col=c("blue", "red", "orange", "green", "gray"), 
       legend=c("all F1 by NH", "F1 with 0 parents assigned", "F1 with 1 parents assigned", "F1 with 2 parents assigned", "F1 with parentage with low probability"), 
       lty = 1, lwd = 2)


## TOT F2 in NH and  NMpi (Pr > 0.8 and Pr < 0.8)
nrow(subset(NH, max_class == "F2"))
nrow(subset(NMpi,Test == "Allenwiller"&  max_class == "F2"& Pr1 > 0.8))
nrow(subset(NMpi,Test == "Allenwiller"&  max_class == "F2"& Pr1 < 0.8))
# classified F2 with no parents
F2_0par = inner_join(subset(NH,max_class == "F2"), subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 & Mo1 == -1 & Fa1 == -1 ), by = "IndivName")
nrow(F2_0par)
# classified F2 with 1 parents
F2_1par = inner_join(subset(NH,max_class == "F2"), subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 & Mo1 == -1 & Fa1 != -1 | Test == "Allenwiller" & Pr1 > 0.8 & Mo1 != -1 & Fa1 == -1), by = "IndivName")
nrow(F2_1par)
# classified F2 with 2 parents
F2_2par = inner_join(subset(NH,max_class == "F2"), subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1), by = "IndivName")
nrow(F2_2par)
# plotting 
plot(vect(subset(NMpi, Test == "Allenwiller"), geom=c("X", "Y"), crs="+init=epsg:3035"))

plot(vect(subset(NMpi,Test == "Allenwiller"& max_class == "F2"),geom=c("X", "Y"), crs="+init=epsg:3035"), add = T, col = "blue") # all F2
plot(vect(F2_0par, geom=c("X.x", "Y.x"), crs="+init=epsg:3035"), add = T, col = "red") # F2 with 0 parents assigned
plot(vect(F2_1par, geom=c("X.x", "Y.x"), crs="+init=epsg:3035"), add = T, col = "orange") # F2 with 1 parents assigned
plot(vect(F2_2par, geom=c("X.x", "Y.x"), crs="+init=epsg:3035"), add = T, col = "green") # F2 with 2 parents assigned
plot(vect(subset(NMpi,Test == "Allenwiller"& max_class == "F2"& Pr1 < 0.8),geom=c("X", "Y"), crs="+init=epsg:3035"), add = T, col = "grey") #"F2 with parentage with low probability"
legend("right", col=c("blue", "red", "orange", "green", "gray"), 
       legend=c("all F2 by NH", "F2 with 0 parents assigned", "F2 with 1 parents assigned", "F2 with 2 parents assigned", "F2 with parentage with low probability"), 
       lty = 1, lwd = 2)

## TOT BC1 in NH and  NMpi (Pr > 0.8 and Pr < 0.8)
nrow(subset(NH, max_class == "BC1"))
nrow(subset(NMpi,Test == "Allenwiller"&  max_class == "BC1"& Pr1 > 0.8))
nrow(subset(NMpi,Test == "Allenwiller"&  max_class == "BC1"& Pr1 < 0.8))
# classified BC1 with no parents
BC1_0par = inner_join(subset(NH,max_class == "BC1"), subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 & Mo1 == -1 & Fa1 == -1 ), by = "IndivName")
nrow(BC1_0par)
# classified BC1 with 1 parents
BC1_1par = inner_join(subset(NH,max_class == "BC1"), subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 & Mo1 == -1 & Fa1 != -1 | Test == "Allenwiller" & Pr1 > 0.8 & Mo1 != -1 & Fa1 == -1), by = "IndivName")
nrow(BC1_1par)
# classified BC1 with 2 parents
BC1_2par = inner_join(subset(NH,max_class == "BC1"), subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1), by = "IndivName")
nrow(BC1_2par)
# plotting 
plot(vect(subset(NMpi, Test == "Allenwiller"), geom=c("X", "Y"), crs="+init=epsg:3035"))

plot(vect(subset(NMpi,Test == "Allenwiller"& max_class == "BC1"),geom=c("X", "Y"), crs="+init=epsg:3035"), add = T, col = "blue") # all BC1
plot(vect(BC1_0par, geom=c("X.x", "Y.x"), crs="+init=epsg:3035"), add = T, col = "red") # BC1 with 0 parents assigned
plot(vect(BC1_1par, geom=c("X.x", "Y.x"), crs="+init=epsg:3035"), add = T, col = "orange") # BC1 with 1 parents assigned
plot(vect(BC1_2par, geom=c("X.x", "Y.x"), crs="+init=epsg:3035"), add = T, col = "green") # BC1 with 2 parents assigned
plot(vect(subset(NMpi,Test == "Allenwiller"& max_class == "BC1"& Pr1 < 0.8),geom=c("X", "Y"), crs="+init=epsg:3035"), add = T, col = "grey") #"BC1 with parentage with low probability"
legend("right", col=c("blue", "red", "orange", "green", "gray"), 
       legend=c("all BC1 by NH", "BC1 with 0 parents assigned", "BC1 with 1 parents assigned", "BC1 with 2 parents assigned", "BC1 with parentage with low probability"), 
       lty = 1, lwd = 2)


### checking the NH classification of all the genealogies offspring types
# NH genotype of the offpsring with no parents assigned by NMpi2
table(subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 & Mo1 == -1 & Fa1 == -1 )$max_class)
# NH genotype of the offpsring with both sylvatica parents assigned by NMpi2
table(subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 & Mo1_Psyl >= 0.8 & Fa1 != -1 & Fa1_Psyl >= 0.8)$max_class)
# NH genotype of the offpsring with both orientalis parents assigned by NMpi2
table(subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 & Mo1_Psyl <= 0.1 & Fa1 != -1 & Fa1_Psyl <= 0.1)$max_class)
# NH genotype of the offpsring with sylvatica mother and orientalis father
table(subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 &  Mo1_Psyl >= 0.8  & Fa1 != -1 & Fa1_Psyl <= 0.1)$max_class)
# NH genotype of the offpsring with orientalis mother and sylvatica father
table(subset(NMpi, Test == "Allenwiller" & Pr1 > 0.8 &  Mo1_Psyl <= 0.1  & Fa1 != -1 & Fa1_Psyl  >= 0.8)$max_class)



##### OLD STUFF
### Allenwiller and Waldi correlation between selection gradients ####
install.packages("corrplot")
install.packages("Hmisc")
library("corrplot")
library(Hmisc)

all_sg = read.xlsx(paste0(path,"Allenwiller/Allenwiller NMpi2 data.xlsx"), sheetIndex = 1)
all_sg = all_sg[which(all_sg$generation == "Adult"), c(40:50)]
plot(all_sg[, !grepl("Z_", names(all_sg))])
matrix = rcorr(as.matrix( all_sg[, !grepl("Z_", names(all_sg))]), type="spearman")
corrplot(matrix$r, method = "circle", type = "upper", order = "hclust", insig = "p-value",
         p.mat = matrix$P, sig.level = 0.005, tl.cex = 1.25, tl.col = "black",
         cl.cex=1.25, title = "")

wal_sg = read.xlsx(paste0(path,"Waldi/Wal_NMpi2_maindata.xlsx"), sheetIndex = 1)
wal_sg = wal_sg[which(wal_sg$generation == "Adult"), c(41:46)]
plot(wal_sg[, !grepl("Z_", names(wal_sg))])
matrix = rcorr(as.matrix( wal_sg[, !grepl("Z_", names(wal_sg))]), type="spearman")
corrplot(matrix$r, method = "circle", type = "upper", order = "hclust", insig = "p-value",
         p.mat = matrix$P, sig.level = 0.005, tl.cex = 1.25, tl.col = "black",
         cl.cex=1.25, title = "")


### overview genealogies tests Allenwiller####

### tests with 0 phenotypic traits
# list of all the files with genealogies results
files_gen = list.files("Allenwiller/Test_0sel_grad/", recursive = T, full.names = T, include.dirs = T, pattern = "*.par")
# reading files into lists with name of the folder and combining into a dataframe
files_gen_ls = lapply(files_gen, function(x) {read.table(file = x, header = T, sep = "", quote = "")})
files_gen_ls = setNames(files_gen_ls, substr(files_gen, 28,33))
gen_df = bind_rows(files_gen_ls, .id = "Test")
# summarizing data for each test
test_sum = gen_df %>% group_by(Test) %>% 
  summarize(passed_gen = sum(Pr1 > 0.8), 
            both_assigned = sum(Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1), 
            single_parent = sum(Pr1 > 0.8 & Mo1 == -1 & Fa1 != -1 | Pr1 > 0.8 & Mo1 != -1 & Fa1 == -1), 
            no_parents = sum(Pr1 > 0.8 & Mo1 == -1 & Fa1 == -1 )) %>% 
  mutate(N_phen_traits = 0)

### tests with 1 phenotypic trait
files_gen = list.files("Allenwiller/Test_1sel_grad/", recursive = T, full.names = T, include.dirs = T, pattern = "*.par")
# reading files into lists with name of the folder and combining into a dataframe
files_gen_ls = lapply(files_gen, function(x) {read.table(file = x, header = T, sep = "", quote = "")})
files_gen_ls = setNames(files_gen_ls, substr(files_gen, 28,33))
gen_df_1 = bind_rows(files_gen_ls, .id = "Test")

# summarizing data for each test
str(gen_df_1)
test_sum_1 = gen_df_1 %>% group_by(Test) %>% 
  summarize(passed_gen = sum(Pr1 > 0.8), 
            both_assigned = sum(Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1), 
            single_parent = sum(Pr1 > 0.8 & Mo1 == -1 & Fa1 != -1 | Pr1 > 0.8 & Mo1 != -1 & Fa1 == -1), 
            no_parents = sum(Pr1 > 0.8 & Mo1 == -1 & Fa1 == -1 )) %>% 
  mutate(N_phen_traits = 1)

### tests with 2 phenotypic traits
files_gen = list.files("Allenwiller/Test_2sel_grad/", recursive = T, full.names = T, include.dirs = T, pattern = "*.par")
# reading files into lists with name of the folder and combining into a dataframe
files_gen_ls = lapply(files_gen, function(x) {read.table(file = x, header = T, sep = "", quote = "")})
files_gen_ls = setNames(files_gen_ls, substr(files_gen, 28,33))
gen_df_2 = bind_rows(files_gen_ls, .id = "Test")

# summarizing data for each test
str(gen_df_2)
test_sum_2 = gen_df_2 %>% group_by(Test) %>% 
  summarize(passed_gen = sum(Pr1 > 0.8), 
            both_assigned = sum(Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1), 
            single_parent = sum(Pr1 > 0.8 & Mo1 == -1 & Fa1 != -1 | Pr1 > 0.8 & Mo1 != -1 & Fa1 == -1), 
            no_parents = sum(Pr1 > 0.8 & Mo1 == -1 & Fa1 == -1 )) %>% 
  mutate(N_phen_traits = 2)

### tests with 3 phenotypic trait
files_gen = list.files("Allenwiller/Test_3sel_grad/", recursive = T, full.names = T, include.dirs = T, pattern = "*.par")
# reading files into lists with name of the folder and combining into a dataframe
files_gen_ls = lapply(files_gen, function(x) {read.table(file = x, header = T, sep = "", quote = "")})
files_gen_ls = setNames(files_gen_ls, substr(files_gen, 28,33))
gen_df_3= bind_rows(files_gen_ls, .id = "Test")

# summarizing data for each test
str(gen_df_3)
test_sum_3 = gen_df_3 %>% group_by(Test) %>% 
  summarize(passed_gen = sum(Pr1 > 0.8), 
            both_assigned = sum(Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1), 
            single_parent = sum(Pr1 > 0.8 & Mo1 == -1 & Fa1 != -1 | Pr1 > 0.8 & Mo1 != -1 & Fa1 == -1), 
            no_parents = sum(Pr1 > 0.8 & Mo1 == -1 & Fa1 == -1 )) %>% 
  mutate(N_phen_traits = 3)

# combining data
geneal_results = rbind(test_sum, test_sum_1, test_sum_2, test_sum_3)
write.csv(geneal_results, "Allenwiller/Allenwiller_NMpi2_genealogies.csv")   

### overview genealogies tests Waldi ####

### tests with 0 phenotypic traits
# list of all the files with genealogies results
files_gen = list.files("Waldi/Test_0sel_grad/", recursive = T, full.names = T, include.dirs = T, pattern = "*.par")
# reading files into lists with name of the folder and combining into a dataframe
files_gen_ls = lapply(files_gen, function(x) {read.table(file = x, header = T, sep = "", quote = "")})
files_gen_ls = setNames(files_gen_ls, substr(files_gen, 22,26))
gen_df = bind_rows(files_gen_ls, .id = "Test")
# summarizing data for each test
test_sum = gen_df %>% group_by(Test) %>% 
  summarize(passed_gen = sum(Pr1 > 0.8), 
            both_assigned = sum(Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1), 
            single_parent = sum(Pr1 > 0.8 & Mo1 == -1 & Fa1 != -1 | Pr1 > 0.8 & Mo1 != -1 & Fa1 == -1), 
            no_parents = sum(Pr1 > 0.8 & Mo1 == -1 & Fa1 == -1 )) %>% 
  mutate(N_phen_traits = 0)

### tests with 1 phenotypic trait
files_gen = list.files("Waldi/Test_1sel_grad/", recursive = T, full.names = T, include.dirs = T, pattern = "*.par")
# reading files into lists with name of the folder and combining into a dataframe
files_gen_ls = lapply(files_gen, function(x) {read.table(file = x, header = T, sep = "", quote = "")})
files_gen_ls = setNames(files_gen_ls, substr(files_gen, 22,27))
gen_df_1 = bind_rows(files_gen_ls, .id = "Test")

# summarizing data for each test
str(gen_df_1)
test_sum_1 = gen_df_1 %>% group_by(Test) %>% 
  summarize(passed_gen = sum(Pr1 > 0.8), 
            both_assigned = sum(Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1), 
            single_parent = sum(Pr1 > 0.8 & Mo1 == -1 & Fa1 != -1 | Pr1 > 0.8 & Mo1 != -1 & Fa1 == -1), 
            no_parents = sum(Pr1 > 0.8 & Mo1 == -1 & Fa1 == -1 )) %>% 
  mutate(N_phen_traits = 1)

### tests with 2 phenotypic traits
files_gen = list.files("Waldi/Test_2sel_grad/", recursive = T, full.names = T, include.dirs = T, pattern = "*.par")
# reading files into lists with name of the folder and combining into a dataframe
files_gen_ls = lapply(files_gen, function(x) {read.table(file = x, header = T, sep = "", quote = "")})
files_gen_ls = setNames(files_gen_ls, substr(files_gen, 22,27))
gen_df_2 = bind_rows(files_gen_ls, .id = "Test")

# summarizing data for each test
str(gen_df_2)
test_sum_2 = gen_df_2 %>% group_by(Test) %>% 
  summarize(passed_gen = sum(Pr1 > 0.8), 
            both_assigned = sum(Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1), 
            single_parent = sum(Pr1 > 0.8 & Mo1 == -1 & Fa1 != -1 | Pr1 > 0.8 & Mo1 != -1 & Fa1 == -1), 
            no_parents = sum(Pr1 > 0.8 & Mo1 == -1 & Fa1 == -1 )) %>% 
  mutate(N_phen_traits = 2)

### tests with 3 phenotypic trait
files_gen = list.files("Waldi/Test_3sel_grad/", recursive = T, full.names = T, include.dirs = T, pattern = "*.par")
# reading files into lists with name of the folder and combining into a dataframe
files_gen_ls = lapply(files_gen, function(x) {read.table(file = x, header = T, sep = "", quote = "")})
files_gen_ls = setNames(files_gen_ls, substr(files_gen, 22,27))
gen_df_3= bind_rows(files_gen_ls, .id = "Test")

# summarizing data for each test
str(gen_df_3)
test_sum_3 = gen_df_3 %>% group_by(Test) %>% 
  summarize(passed_gen = sum(Pr1 > 0.8), 
            both_assigned = sum(Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1), 
            single_parent = sum(Pr1 > 0.8 & Mo1 == -1 & Fa1 != -1 | Pr1 > 0.8 & Mo1 != -1 & Fa1 == -1), 
            no_parents = sum(Pr1 > 0.8 & Mo1 == -1 & Fa1 == -1 )) %>% 
  mutate(N_phen_traits = 3)

# combining data
geneal_results = rbind(test_sum, test_sum_1, test_sum_2, test_sum_3)
write.csv(geneal_results, "Waldi/Waldi_NMpi2_genealogies.csv")

### overview genealogies test Merged ####
merged_gen = read.table(paste0(path, "merged test/Sites_merged_NMpi2.par"), header = T, sep = "", quote = "")
merged_ID  = read.xlsx(paste0(path, "merged test/Sites_merged_NMpi2.xlsx"), sheetIndex = 1)

# adding ID and Site to offspring
merged_gen_ID = merged_gen
merged_gen_ID$Prog_ID = ifelse(merged_gen_ID$Prog == merged_ID[which(merged_ID$X442 == 1), "X1062"], merged_ID[which(merged_ID$X442 == 1), "ID"], NA) 
merged_gen_ID$Prog_Site = ifelse(merged_gen_ID$Prog_ID == merged_ID[which(merged_ID$X442 == 1), "ID"], merged_ID[which(merged_ID$X442 == 1), "Site"], NA) 
# add ID and Site to mother and father
merged_gen_ID$X1062 = merged_gen_ID$Mo1
merged_gen_ID = merged_gen_ID %>% left_join(merged_ID[which(merged_ID$X442 == 0),c(1:6)], by = "X1062")
colnames(merged_gen_ID)[11:15] = c("Mo1_Site", "Mo1_ID", "Mo1_gen", "Mo1_X", "Mo1_Y")
merged_gen_ID$X1062 = merged_gen_ID$Fa1
merged_gen_ID = merged_gen_ID %>% left_join(merged_ID[which(merged_ID$X442 == 0),c(1:6)], by = "X1062")
colnames(merged_gen_ID)[16:20] = c("Fa1_Site", "Fa1_ID", "Fa1_gen", "Fa1_X", "Fa1_Y")

# ANALYSIS # 
# genealogies with P > 0.8
nrow(subset(merged_gen_ID, Pr1 > 0.8))   
# offspring Allenwiller and Waldi
table(merged_gen_ID$Prog_Site)

# Allenwiller offspring with mother and father from Allenwiller or Waldi
nrow(subset(merged_gen_ID, Prog_Site == "Allenwiller" & Mo1_Site == "Allenwiller" & Fa1_Site == "Allenwiller"))
nrow(subset(merged_gen_ID, Prog_Site == "Allenwiller" & Mo1_Site == "Allenwiller" & Fa1_Site == "Waldi"))
nrow(subset(merged_gen_ID, Prog_Site == "Allenwiller" & Mo1_Site == "Waldi" & Fa1_Site == "Allenwiller"))
nrow(subset(merged_gen_ID, Prog_Site == "Allenwiller" & Mo1_Site == "Waldi" & Fa1_Site == "Waldi"))

# Waldi offspring with mother and father from Allenwiller or Waldi
nrow(subset(merged_gen_ID, Prog_Site == "Waldi" & Mo1_Site == "Allenwiller" & Fa1_Site == "Allenwiller"))
nrow(subset(merged_gen_ID, Prog_Site == "Waldi" & Mo1_Site == "Allenwiller" & Fa1_Site == "Waldi"))
nrow(subset(merged_gen_ID, Prog_Site == "Waldi" & Mo1_Site == "Waldi" & Fa1_Site == "Allenwiller"))
nrow(subset(merged_gen_ID, Prog_Site == "Waldi" & Mo1_Site == "Waldi" & Fa1_Site == "Waldi"))

# Allenwiller offspring with both, one (or none) or no parents assigned
nrow(subset(merged_gen_ID, Prog_Site == "Allenwiller" & Mo1 != -1 & Fa1 != -1))
nrow(subset(merged_gen_ID, Prog_Site == "Allenwiller" & Mo1 == -1 | Prog_Site == "Allenwiller" & Fa1 == -1))
nrow(subset(merged_gen_ID, Prog_Site == "Allenwiller" & Mo1 == -1 & Fa1 == -1))

# Waldi offspring with both, one (or none) or no parents assigned
nrow(subset(merged_gen_ID, Prog_Site == "Waldi" & Mo1 != -1 & Fa1 != -1))
nrow(subset(merged_gen_ID, Prog_Site == "Waldi" & Mo1 == -1 | Prog_Site == "Waldi" & Fa1 == -1))
nrow(subset(merged_gen_ID, Prog_Site == "Waldi" & Mo1 == -1 & Fa1 == -1))


### removing orientalis triangle from Allenwiller ####
#importing location trees
trees = read.csv( "C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/NewHybrids/Results/Allenwiller/Allenwiller_coord_circumf_class80ppt.csv")
young = subset(trees,generation!="Adult")
trees = subset(trees,generation=="Adult")
tree_points = terra::vect(trees, geom=c("X", "Y"), crs="+init=epsg:3035")

# selecting orientalis individuals to remove
ori_removed = subset(tree_points,tree_points$Pure2 > 0.8 &
                       grepl("VA0", tree_points$IndivName) &
                       tree_points$IndivName != c("VA001", "VA002", "VA003", "VA004", "VA033", "VA034",
                                                  "VA035", "VA036", "VA037", "VA038", "VA050", "VA053",
                                                  "VA054", "VA056", "VA057", "VA058", "VA059", "VA060", "VA061" ))
# visual check (need to load raster library)
plot(tree_points)
plot(ori_removed, add = T, col = "red")

# read NMpi2 input file to remove individuals
file = read.xlsx(paste0(path,"Allenwiller/Allenwiller NMpi2 data.xlsx"), sheetIndex = 1)
file_removed = vect(subset(file, !(ID %in% ori_removed$IndivName)), geom=c("X", "Y"), crs="+init=epsg:3035")
plot(file_removed[which(file_removed$generation == "Adult")])

write.xlsx(as.data.frame(file_removed,geom="XY"), paste0(path,"Allenwiller/Test_ori_removed/All_NMPi2_input_ori_removed2.xlsx"))

#### analysing the results ##

## final tests results check ##


### Allenwiller and Waldi Bayesian results - checking convergence ####
library(reshape2)
data = read.table(paste0(path, "Allenwiller/Bayesian_test/Allenwiller_NMpi2_input_DBH_Psyl_compcoeff.MCMC"), header = FALSE)
# extract iterations and parameters
iter = data$V1
param = data[, 5:ncol(data)]
# converting to long format for plotting 
long_data = melt(param)
long_data$iter = rep(iter, times = ncol(param))

# trace plots for each parameter
ggplot(long_data, aes(x = iter, y = value, color = variable)) +
  geom_line() + guides(color = "none")+
  facet_wrap(~variable, scales = "free_y") +
  theme_minimal() +
  labs(title = "Trace Plots of MCMC Parameters", x = "Iteration", y = "Parameter Value")


### Allenwiller and Waldi Bayesian results - genealogies ####
all_par = read.table(paste0(path, "Allenwiller/Bayesian_test/Allenwiller_NMpi2_input_DBH_Psyl_compcoeff.par"), header = T)
wal_par = read.table(paste0(path, "Waldi/Bayesian_test/Waldi_NMpi2_input_DBH_Psyl_compcoeff.par"), header = T)
all_par$Test = "Allenwiller"
wal_par$Test = "Waldi"

# adding the ID and coordinates of offspring and DBH, genotype classification and coordinates of the parents 
all_ID  = read.xlsx(paste0(path, "Allenwiller/Allenwiller_NMpi2_data.xlsx"), sheetIndex = 1)
wal_ID  = read.xlsx(paste0(path, "Waldi/Waldi_NMpi2_data.xlsx"), sheetIndex = 1)

all_par_ID = all_par
all_par_ID$Prog_ID = ifelse(all_par_ID$Prog == all_ID[which(all_ID$generation == "Offspring"), "ID_number"], all_ID[which(all_ID$generation == "Offspring"), "ID"], NA) 
all_par_ID$Prog_X = ifelse(all_par_ID$Prog == all_ID[which(all_ID$generation == "Offspring"), "ID_number"], all_ID[which(all_ID$generation == "Offspring"), "X"], NA) 
all_par_ID$Prog_Y = ifelse(all_par_ID$Prog == all_ID[which(all_ID$generation == "Offspring"), "ID_number"], all_ID[which(all_ID$generation == "Offspring"), "Y"], NA) 
all_par_ID$ID_number = all_par_ID$Mo1
all_par_ID = all_par_ID %>% left_join(all_ID[which(all_ID$generation == "Adult"),c(3:6,41,45)], by = "ID_number")
colnames(all_par_ID)[13:17] = c("Mo1_ID","Mo1_X", "Mo1_Y", "Mo1_DBH", "Mo1_Psyl")
all_par_ID$ID_number = all_par_ID$Fa1
all_par_ID = all_par_ID %>% left_join(all_ID[which(all_ID$generation == "Adult"),c(3:6,41,45)], by = "ID_number")
colnames(all_par_ID)[18:22] = c("Fa1_ID", "Fa1_X", "Fa1_Y","Fa1_DBH","Fa1_Psyl")

wal_par_ID = wal_par
wal_par_ID$Prog_ID = ifelse(wal_par_ID$Prog == wal_ID[which(wal_ID$generation == "Offspring"), "ID_number"], wal_ID[which(wal_ID$generation == "Offspring"), "ID"], NA) 
wal_par_ID$Prog_X = ifelse(wal_par_ID$Prog == wal_ID[which(wal_ID$generation == "Offspring"), "ID_number"], wal_ID[which(wal_ID$generation == "Offspring"), "X"], NA) 
wal_par_ID$Prog_Y = ifelse(wal_par_ID$Prog == wal_ID[which(wal_ID$generation == "Offspring"), "ID_number"], wal_ID[which(wal_ID$generation == "Offspring"), "Y"], NA) 
wal_par_ID$ID_number = wal_par_ID$Mo1
wal_par_ID = wal_par_ID %>% left_join(wal_ID[which(wal_ID$generation == "Adult"),c(1,5:7,41,43)], by = "ID_number")
colnames(wal_par_ID)[13:17] = c("Mo1_ID","Mo1_X", "Mo1_Y", "Mo1_DBH", "Mo1_Psyl")
wal_par_ID$ID_number = wal_par_ID$Fa1
wal_par_ID = wal_par_ID %>% left_join(wal_ID[which(wal_ID$generation == "Adult"),c(1,5:7,41,43)], by = "ID_number")
colnames(wal_par_ID)[18:22] = c("Fa1_ID", "Fa1_X", "Fa1_Y","Fa1_DBH","Fa1_Psyl")

par_results = rbind(all_par_ID, wal_par_ID)
write.csv(par_results, "Parentage_results_Bayesian.csv")

# write.xlsx(par_results, "NMpi2_overall_results.xlsx", sheetName="Bayesian_genealogies", append=TRUE)

# summarizing results
par_results = read.xlsx("NMpi2_overall_results.xlsx", sheetName = "Bayesian_genealogies")

par_sum = par_results %>% group_by(Test) %>% 
  summarise(not_passed_gen = sum(Pr1 < 0.8),
            # genealogies with probability lower than 0.8
            passed_gen = sum(Pr1 > 0.8), 
            # genealogies with probability higher than 0.8
            both_assigned = sum(Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1), 
            # offspring with both parents assigned
            single_parent = sum((Pr1 > 0.8 & Mo1 == -1 & Fa1 != -1) | (Pr1 > 0.8 & Mo1 != -1 & Fa1 == -1), na.rm = TRUE),
            # offspring with one parent assigned
            no_parents = sum(Pr1 > 0.8 & Mo1 == -1 & Fa1 == -1, na.rm = TRUE ), 
            # offspring with no parents assigned
            N_mo1_syl = length(unique(Mo1[Mo1_Psyl > 0.8 & !is.na(Mo1_Psyl)])),
            # number of sylvatica mothers
            N_mo1_ori = length(unique(Mo1[Mo1_Psyl < 0.1 & !is.na(Mo1_Psyl)])),
            # number of orientalis mothers
            N_fa1_syl = length(unique(Fa1[Fa1_Psyl > 0.8 & !is.na(Fa1_Psyl)])),
            # number of sylvatica fathers
            N_fa1_ori = length(unique(Fa1[Fa1_Psyl < 0.1 & !is.na(Fa1_Psyl)])),
            # number of orientalis fathers
            N_syl_syl = sum(Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1 & Mo1_Psyl >= 0.8 & Fa1_Psyl >= 0.8, na.rm = TRUE),
            # number of offspring with both sylvatica parents
            N_ori_ori = sum(Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1 & Mo1_Psyl <= 0.1 & Fa1_Psyl <= 0.1, na.rm = TRUE),
            # number of offspring with both orientalis parents
            N_syl_ori = sum(Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1 & Mo1_Psyl >= 0.8 & Fa1_Psyl <= 0.1, na.rm = TRUE),
            # number of offspring with sylvatica mother and orientalis father
            N_ori_syl = sum(Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1 & Mo1_Psyl <= 0.1 & Fa1_Psyl >= 0.8, na.rm = TRUE)
            # number of offspring with orientalis mother and sylvatica mother
  ) 

# plotting the offspring with both assigned parents
par(mfrow = c(1,2))
all_off_bothassigned = subset(par_results, Test == "Allenwiller"& Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1)$Prog_ID
plot(vect(all_ID,geom=c("X", "Y"), crs="+init=epsg:3035"), main = "Allenwiller - offspring with both assigned parents")
plot(vect(all_ID[which(all_ID$ID %in% all_off_bothassigned),],geom=c("X", "Y"), crs="+init=epsg:3035"), add = T, col = "red")

wal_off_bothassigned = subset(par_results, Test == "Waldi"& Pr1 > 0.8 & Mo1 != -1 & Fa1 != -1)$Prog_ID
plot(vect(wal_ID,geom=c("X", "Y"), crs="+init=epsg:3035"), main = "Waldi - offspring with both assigned parents")
plot(vect(wal_ID[which(wal_ID$ID %in% wal_off_bothassigned),],geom=c("X", "Y"), crs="+init=epsg:3035"), add = T, col = "red")

# N offspring for each parent divided by subspecies and adding as Sheet to result file
Noff_mothers = par_results %>% mutate(Mo1_sp = case_when(Mo1 != -1 & Mo1_Psyl > 0.8 ~ "Sylvatica",
                                                         Mo1 != -1 & Mo1_Psyl < 0.1 ~ "Orientalis", 
                                                         Mo1 == -1 ~ NA))%>% 
  group_by(Test,Mo1_ID, Mo1_X, Mo1_Y, Mo1_DBH, Mo1_sp ) %>% 
  summarise(N_off = n()) %>% 
  ungroup()
write.xlsx(Noff_mothers, "NMpi2_overall_results.xlsx", sheetName="Bayesian_Mot_Noff", append=TRUE)

Noff_fathers = par_results %>% mutate(Fa1_sp = case_when(Fa1 != -1 & Fa1_Psyl > 0.8 ~ "Sylvatica",
                                                         Fa1 != -1 & Fa1_Psyl < 0.1 ~ "Orientalis", 
                                                         Fa1 == -1 ~ NA))%>% 
  group_by(Test,Fa1_ID, Fa1_X, Fa1_Y, Fa1_DBH, Fa1_sp ) %>% 
  summarise(N_off = n()) %>% 
  ungroup()

write.xlsx(Noff_fathers, "NMpi2_overall_results.xlsx", sheetName="Bayesian_Fat_Noff", append=TRUE)

# summary statistics mothers 
Noff_mothers = read.xlsx("NMpi2_overall_results.xlsx", sheetName = "Bayesian_Mot_Noff")
Noff_fathers = read.xlsx("NMpi2_overall_results.xlsx", sheetName = "Bayesian_Fat_Noff")

par(mfrow = c(2,2))
# number of offspring for sylvatica and orientalis mothers
Noff_mothers$N_off = as.numeric(as.character(Noff_mothers$N_off))
hist(Noff_mothers[which(Noff_mothers$Test == "Allenwiller" & Noff_mothers$Mo1_sp == "Sylvatica"),]$N_off)
hist(Noff_mothers[which(Noff_mothers$Test == "Allenwiller" & Noff_mothers$Mo1_sp == "Orientalis"),]$N_off)
hist(Noff_mothers[which(Noff_mothers$Test == "Waldi" & Noff_mothers$Mo1_sp == "Sylvatica"),]$N_off)
hist(Noff_mothers[which(Noff_mothers$Test == "Waldi" & Noff_mothers$Mo1_sp == "Orientalis"),]$N_off)

Noff_mothers_sum = Noff_mothers %>% group_by(Test, Mo1_sp) %>%  
  summarise(N = n(), Mean_Noff = mean(N_off), SD_Noff = sd(N_off))

# test of the means between N offspring of mother sylvatica and mother orientalis
wilcox.test(Noff_mothers[which(Noff_mothers$Mo1_sp == "Sylvatica"),]$N_off, 
            Noff_mothers[which(Noff_mothers$Mo1_sp == "Orientalis"),]$N_off)

wilcox.test(Noff_mothers[which(Noff_mothers$Test == "Waldi" & Noff_mothers$Mo1_sp == "Sylvatica"),]$N_off,
            Noff_mothers[which(Noff_mothers$Test == "Waldi" & Noff_mothers$Mo1_sp == "Orientalis"),]$N_off)

wilcox.test(Noff_mothers[which(Noff_mothers$Test == "Allenwiller" & Noff_mothers$Mo1_sp == "Sylvatica"),]$N_off, 
            Noff_mothers[which(Noff_mothers$Test == "Allenwiller" & Noff_mothers$Mo1_sp == "Orientalis"),]$N_off)

# summary statistics fathers 
Noff_fathers$N_off = as.numeric(as.character(Noff_fathers$N_off))
hist(Noff_fathers[which(Noff_fathers$Test == "Allenwiller" & Noff_fathers$Fa1_sp == "Sylvatica"),]$N_off)
hist(Noff_fathers[which(Noff_fathers$Test == "Allenwiller" & Noff_fathers$Fa1_sp == "Orientalis"),]$N_off)
hist(Noff_fathers[which(Noff_fathers$Test == "Waldi" & Noff_fathers$Fa1_sp == "Sylvatica"),]$N_off)
hist(Noff_fathers[which(Noff_fathers$Test == "Waldi" & Noff_fathers$Fa1_sp == "Orientalis"),]$N_off)

Noff_fathers_sum = Noff_fathers %>% group_by(Test, Fa1_sp) %>%
  summarise(N = n(), 
            Mean_Noff = mean(N_off),
            SD_Noff = sd(N_off))

# test of the means between N offspring of mother sylvatica and mother orientalis
wilcox.test(Noff_fathers[which(Noff_fathers$Fa1_sp == "Sylvatica"),]$N_off, 
            Noff_fathers[which(Noff_fathers$Fa1_sp == "Orientalis"),]$N_off)

wilcox.test(Noff_fathers[which(Noff_fathers$Test == "Waldi" & Noff_fathers$Fa1_sp == "Sylvatica"),]$N_off,
            Noff_fathers[which(Noff_fathers$Test == "Waldi" & Noff_fathers$Fa1_sp == "Orientalis"),]$N_off)

wilcox.test(Noff_fathers[which(Noff_fathers$Test == "Allenwiller" & Noff_fathers$Fa1_sp == "Sylvatica"),]$N_off, 
            Noff_fathers[which(Noff_fathers$Test == "Allenwiller" & Noff_fathers$Fa1_sp == "Orientalis"),]$N_off)

# plotting Mothers Allenwiller and Waldi
Noff_mothers[4:6] = sapply(Noff_mothers[4:6], as.numeric)
Noff_fathers[4:6] = sapply(Noff_fathers[4:6], as.numeric)

hist(Noff_mothers$N_off)
size_scale = scale_size_continuous(limits=c(1,80),breaks=c(5,20,40,60,80))

ggplot() + theme_bw() +
  geom_point(subset(Noff_mothers, Test != "Allenwiller_ori_rem"), mapping = aes(x = Mo1_X, y = Mo1_Y, size = N_off, col = Mo1_sp)) + 
  scale_color_manual(values = c("Sylvatica" = "#F0E442", "Orientalis" = "#9E9AC8"))+size_scale+
  #  geom_text(data = subset(Noff_mothers, Test != "Waldi"), aes(x = X, y = Y,label = N_off), vjust = 1.5, color = "black", size = 3) +
  facet_wrap(~Test, scales = "free") + labs(title = "NMpi2 genealogies - Mother trees. Size = N offspring") 
theme(legend.position="none")

# plotting Fathers Allenwiller and Waldi
which.max(Noff_fathers$N_off)
size_scale = scale_size_continuous(limits=c(1,30),breaks=c(1,5,10,15,20))

ggplot() + theme_bw() + 
  geom_point(subset(Noff_fathers, Test != "Allenwiller_ori_rem"), mapping = aes(x = Fa1_X, y = Fa1_Y, size = N_off, col = Fa1_sp)) + 
  scale_color_manual(values = c("Sylvatica" = "#F0E442", "Orientalis" = "#9E9AC8")) +size_scale+
  #   geom_text(data = subset(Noff_fathers, Test != "Waldi"), aes(x = X, y = Y,label = N_off), vjust = 1.5, color = "black", size = 3) +
  facet_wrap(~Test, scales = "free") + labs(title = "NMpi2 genealogies - Father trees. Size = N offspring") 


###  correlation between Number of offspring and other traits (DBH and fitness traits) #####
install.packages("smplot2")
library(smplot2)

### correlation between number of offspring and DBH - age of the trees
Noff_mothers = read.xlsx(paste0(path, "/NMpi2_overall_results.xlsx"), sheetName = "ML_Mot_Noff")
# removing extra test and unknown mothers
Noff_mothers = Noff_mothers[which(Noff_mothers$Test != "Allenwiller_ori_rem" & Noff_mothers$Mo1_ID != "NA"),]

# calculating the tree age (converting cm to inches and multiplying for growth factor for beech)
Noff_mothers$Mo1_age = as.numeric(Noff_mothers$Mo1_DBH) * 0.393701 * 4

ggplot(data = Noff_mothers, aes(x=as.numeric(Mo1_age),  y=N_off, col = Mo1_sp))+
  scale_color_manual(values = c("Sylvatica" = "#F0E442", "Orientalis" = "#9E9AC8"))+
  geom_point(stat='identity')+
  ylab("Number of offspring")+
  xlab("Tree age ")+
  facet_wrap(~Test, scales="free_y", ncol=1)+sm_statCorr()

# creating different age classes and checking if the fecundity is different
install.packages("multcomp")
library(multcomp)

# creating age bins 
# calculating back the original number of seedlings: the assigned offspring is only about 60% of the total, and the total is only the genotyped seedlings, which are 20% of the seedlings in the area unit (circle plot)
summary(Noff_mothers$Mo1_age)
# proportion of assigned Allenwiller and Waldi on the total
prop_assigned = ( (340 + 300) / (639 + 423) ) * 100

Noff_mothers = Noff_mothers %>% 
  mutate(AgeClass = cut(Mo1_age,breaks = c(min(Noff_mothers$Mo1_age),50, 70, 90, 110, max(Noff_mothers$Mo1_age)),include.lowest = T), 
         Noff_original = (N_off * 100 / prop_assigned ) * 100  / 20 )

# check potential differences between subspecies
ggplot(subset(Noff_mothers, Test !="Allenwiller_ori_rem"), aes(x = AgeClass, y = Noff_original, group = Mo1_sp, fill = Mo1_sp))+
  scale_fill_manual(values = c("Sylvatica" = "#F0E442", "Orientalis" = "#9E9AC8"))+
  stat_summary(fun.y = mean, geom = "bar",width=0.6, position=position_dodge(width = 0.7)) + theme_bw() +
  stat_summary(fun.data = mean_se, geom = "errorbar", width=.2,  position=position_dodge(width = 0.7))+
  stat_compare_means(method = "t.test", label.x = 1.4, label.y = 25, aes(label = ..p.signif..))+
  ylab("Number of Offspring of Mother trees")


# ANOVA to check if the fecundity is different
anova = aov(Noff_original~AgeClass, data = Noff_mothers )
summary(anova)
# Tukey's HSD post-hoc test 
tukey = TukeyHSD(anova)
print(tukey)
par(mar = c(2,10,3,4), las = 1)  # las = 1 makes labels horizontal
plot(tukey)

# mean number of offspring for each age class
ggplot(subset(Noff_mothers, Test !="Allenwiller_ori_rem"), aes(x = AgeClass, y = Noff_original))+
  stat_summary(fun.y = mean, geom = "bar",width=0.6, position=position_dodge(width = 0.7)) + theme_bw() +
  stat_summary(fun.data = mean_se, geom = "errorbar", width=.2,  position=position_dodge(width = 0.7))+
  stat_compare_means(method = "anova", label.x = 1.4, label.y = 100, aes(label = ..p.signif..))+
  ylab("Number of Offspring of Mother trees")

#### from these results, seems reasonable to split the population into age classes 30-70, 70-110, 110-max
Noff_mothers = Noff_mothers %>% 
  mutate(AgeClass_final = cut(Mo1_age,breaks = c(min(Noff_mothers$Mo1_age, na.rm = T),70, 110, max(Noff_mothers$Mo1_age, na.rm = T)), include.lowest = T))

# mean number of offspring for each age class
ggplot(subset(Noff_mothers, Test !="Allenwiller_ori_rem"), aes(x = AgeClass_final, y = Noff_original))+
  stat_summary(fun.y = mean, geom = "bar",width=0.6, position=position_dodge(width = 0.7)) + theme_bw() +
  stat_summary(fun.data = mean_se, geom = "errorbar", width=.2,  position=position_dodge(width = 0.7))+
  stat_compare_means(method = "anova", label.x = 2, label.y = 100, aes(label = ..p.signif..))+
  ylab("Number of Offspring of Mother trees")

# ANOVA to check if the fecundity is different
anova = aov(Noff_original~AgeClass_final, data = Noff_mothers )
summary(anova)
# Tukey's HSD post-hoc test 
tukey = TukeyHSD(anova)
print(tukey)
par(mar = c(2,10,3,4), las = 1)  # las = 1 makes labels horizontal
plot(tukey)

#### mean N offspring for NEMO 

Noff_mothers_sum = Noff_mothers %>% group_by(AgeClass_final) %>% 
  summarise(mean_Noff = mean(Noff_original, na.rm = T), n = n())

# correlation number of offspring and fitness traits ###
traits = read.xlsx("C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Data/Dodorico_Data_for_LME_model.xlsx", sheetIndex = 1)
# changing column names for dataset combination
colnames(Noff_mothers)[2] = "tree_ID"
colnames(Noff_fathers)[2] = "tree_ID"

# taking the mean of the measurements of july and august
traits = traits %>% group_by(site, year, tree_ID) %>%  summarise(across(5:14, \(x) mean(x, na.rm = TRUE)))
traits[traits =="Allen"] = "Allenwiller"
traits = left_join(traits, Noff_mothers[,2:4], by = "tree_ID", relationship = "many-to-many")
traits = left_join(traits, Noff_fathers[,2:4], by = "tree_ID", relationship = "many-to-many")
colnames(traits)[14:17] = c("Noff_mother", "Mother_sp", "Noff_father", "Father_sp")

write.csv(traits, "C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Parentage analysis/NMpi2/Camilla_data/Traits_fecundity.csv")

# checking correlation between all the traits
plot(traits[, c(4:14, 16)])

# correlogram with white cells = not significant correlation, number = coefficient
### allenwiller
matrix_all = rcorr(as.matrix(traits[traits$site == "Allenwiller", c(4:14, 16)], type="spearman"))
diag(matrix_all$P) = diag(matrix_all$r)
setwd("C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Parentage analysis/NMpi2/Camilla_data/Figures/")
png(filename=paste0(path, "Correl_All_Noff_traits.png"), width =2000,height=1400, res=200)
corrplot(matrix_all$r, method = "number", type = "upper", order = "hclust", insig = "blank",
         p.mat = matrix_all$P, sig.level = 0.005, tl.cex = 1.25, tl.col = "black", cl.cex=1.25, title = "")
dev.off()
### waldi
matrix_wal = rcorr(as.matrix(traits[traits$site == "Waldi", c(4:14, 16)], type="spearman"))
diag(matrix_wal$P) = diag(matrix_wal$r)
png(filename=paste0(path, "Correl_Wal_Noff_traits.png"), width =2000,height=1400, res=200)
corrplot(matrix_wal$r, method = "number", type = "upper", order = "hclust", insig = "blank",
         p.mat = matrix_wal$P, sig.level = 0.005, tl.cex = 1.25, tl.col = "black", cl.cex=1.25, title = "")
dev.off()
# sylvatica (mother trees)
matrix_syl = rcorr(as.matrix(traits[traits$Mother_sp == "Sylvatica", c(4:14, 16)], type="spearman"))
diag(matrix_syl$P) = diag(matrix_syl$r)
png(filename=paste0(path, "Correl_Syl_Noff_traits.png"), width =2000,height=1400, res=200)
corrplot(matrix_syl$r, method = "number", type = "upper", order = "hclust", insig = "blank",
         p.mat = matrix_syl$P, sig.level = 0.005, tl.cex = 1.25, tl.col = "black", cl.cex=1.25, title = "")
dev.off()
# orientalis (mother trees)
matrix_ori = rcorr(as.matrix(traits[traits$Mother_sp == "Orientalis", c(4:14, 16)], type="spearman"))
diag(matrix_ori$P) = diag(matrix_ori$r)
png(filename=paste0(path, "Correl_Ori_Noff_traits.png"), width =2000,height=1400, res=200)
corrplot(matrix_ori$r, method = "number", type = "upper", order = "hclust", insig = "blank",
         p.mat = matrix_ori$P, sig.level = 0.005, tl.cex = 1.25, tl.col = "black", cl.cex=1.25, title = "")
dev.off()



### removing surrounding syl in Allenwiller, keeping only the core ####
dat <- read.xlsx(paste0(path, "Allenwiller/Allenwiller_NMPi2_data.xlsx"), sheetName = "all_data")
dat_v <- vect(dat, geom=c("X", "Y"), crs="+init=epsg:3035")
plot(dat_v)
# removing trees in the West and in the South part of the plot
hist(dat$X)
subset <- subset(dat, X>4126650 & Y >2840050)
subset_v <- vect(subset, geom=c("X", "Y"), crs="+init=epsg:3035")

# visual check
plot(dat_v)
plot(subset_v, add = T, col="red")

# write the reduced dataframe
write.xlsx(subset, paste0(path,"Allenwiller/Allenwiller_NMPi2_data.xlsx"), sheetName = "reduced_data", append = T)

### analysis ##

