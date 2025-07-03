# analysis of SSR panel with hybriddetective #

# installing packages in a directory where you have permissions
.libPaths("C:/Users/stefanin/R/win-library/4.4")

devtools::install_github("bwringe/parallelnewhybrid") #required
devtools::install_github("kkeenan02/diveRsity", dependencies = TRUE)
devtools::install_github("rystanley/genepopedit", dependencies = TRUE) #required to convert to NH format
devtools::install_github("bwringe/hybriddetective", dependencies = TRUE) #The original package has a bug, install the one below
devtools::install_github("stevemussmann/hybriddetective", force = TRUE) ## fixed version (found here https://github.com/bwringe/hybriddetective/issues/12 )

library("parallel")
library("adegenet")
library("diveRsity")
library("genepopedit") ## needed for converting files
library("parallelnewhybrid")
library("hybriddetective")

path="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Data_clean"
setwd(dir="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Data_clean")
path.nh <- "C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/NewHybrids/packages/newhybrids-master/bin/PC/newhybrids/"

####  identifying the most informative loci for panel analysis ####

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
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip coordinates for easier reading
  labs(
    title = "Informativeness of Loci",
    x = "Locus",
    y = "Informativeness Score"
  ) +
  theme_minimal()
ssr_panel

png(filename="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Manuscript/Figures/SSR_panel_informativescore.png", 
    width =8000,height=8000, res=1500)
plot(ssr_panel)
dev.off()


#### preparing datasets for validation ####
## Waldi
dir.create("Hybriddetective")

# load genotype of all the individuals
df_val <- read.xlsx("NewHybrids/Waldi_all.xlsx", sheetIndex = 1)
#load circumference data
circumference <- read.csv("Waldi_indiv_info.csv")
df_val <- merge(df_val, circumference[, c("SampleID", "circumference")], by = "SampleID")
# select the adult sylvatica with biggest circumference
syl <- df_val[grep("^A", df_val$SampleID), ]
biggest_syl <- syl[order(syl$circumference, decreasing = TRUE), ][1:30, ]
# select orietnalis
ori <- df_val[grep("^O", df_val$SampleID), ]
# create reference dataframe
reference <- rbind(biggest_syl, subset(ori, gen == "Adult"))

# creating genepop and nh format
reference <- reference[,-c(2:4, 37)] #removing extra columns
for (j in 2:33) {
  for (i in 1:nrow(reference)){
    value <- as.numeric(reference[i,j])
    if (value < 100){
      reference[i,j]= as.character(paste0(0,value))
    } else {
      next
    }
  }
}

#merging loci columns
reference$casolfagus_29 <- paste0(reference$casolfagus_29,reference$NA..1)
reference$concat14 <- paste0(reference$concat14,reference$NA..2)
reference$csolfagus_05 <- paste0(reference$csolfagus_05,reference$NA..3)
reference$csolfagus_06 <- paste0(reference$csolfagus_06,reference$NA..4)
reference$csolfagus_19 <- paste0(reference$csolfagus_19,reference$NA..5)
reference$csolfagus_31 <- paste0(reference$csolfagus_31, reference$NA..6)
reference$DE576 <- paste0(reference$DE576,reference$NA..7)
reference$DUKCT <- paste0(reference$DUKCT,reference$NA..8) 
reference$DZ447 <- paste0(reference$DZ447,reference$NA..9)
reference$EEU75 <- paste0(reference$EEU75,reference$NA..10)
reference$EJV8T <- paste0(reference$EJV8T,reference$NA..11)
reference$EMILY <- paste0(reference$EMILY,reference$NA..12)
reference$ERHBI <- paste0(reference$ERHBI, reference$NA..13)
reference$FS1_15 <- paste0(reference$FS1_15,reference$NA..14)
reference$sfc_0036 <- paste0(reference$sfc_0036, reference$NA..15)
reference$sfc_1143 <- paste0(reference$sfc_1143, reference$NA..16)
reference <- reference[, -grep("NA.", colnames(reference))]
reference[reference=="0-10-1"] <- as.character("000000") #missing data

## adding a line to set orientalis as population 2
add_pop_line <- function(filename) {
  # Read the file
  lines <- readLines(filename)
  
  # Find the first line that starts with "O" (this will be the divide between A and O samples)
  o_start <- which(grepl("^O", lines))[1]
  
  # Create new lines with "Pop" inserted before O samples
  new_lines <- c(
    lines[1:(o_start-1)],  # A samples
    "Pop",                 # Add Pop line
    lines[o_start:length(lines)]  # O samples
  )
  
  # Write the modified file back
  writeLines(new_lines, filename)
}

# creating files only with most and less informative loci
reference_most <- data.frame(SampleID = reference$SampleID)
reference_most <- cbind(reference_most, reference[, c(colnames(reference) %in% most_inf)])
reference_less <- data.frame(SampleID = reference$SampleID)
reference_less <- cbind(reference_less, reference[, c(colnames(reference) %in% less_inf)])

## creating files in the respective directories
dir.create("Hybriddetective/Waldi_panel_all")
dir.create("Hybriddetective/Waldi_panel_most")
dir.create("Hybriddetective/Waldi_panel_less")

# Genepop file of all the reference files for validation
genepopedit::genepop_unflatten(reference, "Hybriddetective/Waldi_panel_all/Waldi_genepop_validation.txt")
add_pop_line("Hybriddetective/Waldi_panel_all/Waldi_genepop_validation.txt")
genepopedit::genepop_unflatten(reference_most, "Hybriddetective/Waldi_panel_most/Waldi_genepop_validation_locimost.txt")
add_pop_line("Hybriddetective/Waldi_panel_most/Waldi_genepop_validation_locimost.txt")
genepopedit::genepop_unflatten(reference_less, "Hybriddetective/Waldi_panel_less/Waldi_genepop_validation_lociless.txt")
add_pop_line("Hybriddetective/Waldi_panel_less/Waldi_genepop_validation_lociless.txt")


### Allenwiller 
# load Kurz results 
kurz_class <- read.xlsx("20201030_VosgesAll(V9_C2)_ClumppResults_Species.xlsx", sheetIndex = 1)
# load my dataset and merge
df_val <- read.csv("Allenwiller_indiv_info.csv")
df_val <- merge(df_val, kurz_class[, c(1,6)], by.x = "SampleID", by.y = "ProbeID", all.x = T)
# select biggest orientalis based on Kurz initial classification
syl <- subset(df_val, generation == "Adult" & Species == "Sylvatica")
biggest_syl <- syl[order(syl$circumference, decreasing = TRUE), ][1:30, ]
# select biggest orientalis from the triangle (we are more sure they are orientalis) and based on Kurz classification
ori <-df_val[grep("VA0", df_val$SampleID), ]
ori <- subset(ori, generation == "Adult" & Species == "Orientalis")
biggest_ori <- ori[order(ori$circumference, decreasing = TRUE), ][1:30, ]
# create reference dataframe
referenceID <- c(biggest_syl$SampleID, biggest_ori$SampleID)

# load genotypes of all the individuals
df_all <- read.xlsx("NewHybrids/Allenwiller_all.xlsx", sheetIndex = 1)
# subset only the reference individuals
reference <- df_all[(df_all$SampleID %in% referenceID),]
# reorder with sylvatica first 
reference <- reference[order(reference$subspecies, decreasing = TRUE),]

reference <- reference[,-c(1,3:6)] #removing extra columns
# converting to numeric and filling with 0s
for (j in 2:33) {
  for (i in 1:nrow(reference)){
    value <- as.numeric(reference[i,j])
    if (value < 100){
      reference[i,j]= as.character(paste0(0,value))
    } else {
      next
    }
  }
}

#merging loci columns
reference$casolfagus_29 <- paste0(reference$casolfagus_29,reference$NA..1)
reference$concat14 <- paste0(reference$concat14,reference$NA..2)
reference$csolfagus_05 <- paste0(reference$csolfagus_05,reference$NA..3)
reference$csolfagus_06 <- paste0(reference$csolfagus_06,reference$NA..4)
reference$csolfagus_19 <- paste0(reference$csolfagus_19,reference$NA..5)
reference$csolfagus_31 <- paste0(reference$csolfagus_31, reference$NA..6)
reference$DE576 <- paste0(reference$DE576,reference$NA..7)
reference$DUKCT <- paste0(reference$DUKCT,reference$NA..8) 
reference$DZ447 <- paste0(reference$DZ447,reference$NA..9)
reference$EEU75 <- paste0(reference$EEU75,reference$NA..10)
reference$EJV8T <- paste0(reference$EJV8T,reference$NA..11)
reference$EMILY <- paste0(reference$EMILY,reference$NA..12)
reference$ERHBI <- paste0(reference$ERHBI, reference$NA..13)
reference$FS1_15 <- paste0(reference$FS1_15,reference$NA..14)
reference$sfc_0036 <- paste0(reference$sfc_0036, reference$NA..15)
reference$sfc_1143 <- paste0(reference$sfc_1143, reference$NA..16)
reference <- reference[, -grep("NA.", colnames(reference))]
reference[reference=="0-10-1"] <- as.character("000000") #missing data

# function to wirte the Pop line to separate sylvatica from orientalis
add_pop_line <- function(filename) {
  # Read the file
  lines <- readLines(filename)
  
  # Find the first line that starts with "O" (this will be the divide between A and O samples)
  o_start <- which(grepl("VA0", lines))[1]
  
  # Create new lines with "Pop" inserted before O samples
  new_lines <- c(
    lines[1:(o_start-1)],  # A samples
    "Pop",                 # Add Pop line
    lines[o_start:length(lines)]  # O samples
  )
  
  # Write the modified file back
  writeLines(new_lines, filename)
}

# creating files only with most and less informative loci
reference_most <- data.frame(SampleID = reference$SampleID)
reference_most <- cbind(reference_most, reference[, c(colnames(reference) %in% most_inf)])
reference_less <- data.frame(SampleID = reference$SampleID)
reference_less <- cbind(reference_less, reference[, c(colnames(reference) %in% less_inf)])

## creating files in the respective directories
dir.create("Hybriddetective/Allenwiller_panel_all")
dir.create("Hybriddetective/Allenwiller_panel_most")
dir.create("Hybriddetective/Allenwiller_panel_less")

# Genepop file of all the reference files for validation
genepopedit::genepop_unflatten(reference, "Hybriddetective/Allenwiller_panel_all/Allenwiller_genepop_validation.txt")
add_pop_line("Hybriddetective/Allenwiller_panel_all/Allenwiller_genepop_validation.txt")
genepopedit::genepop_unflatten(reference_most, "Hybriddetective/Allenwiller_panel_most/Allenwiller_genepop_validation_locimost.txt")
add_pop_line("Hybriddetective/Allenwiller_panel_most/Allenwiller_genepop_validation_locimost.txt")
genepopedit::genepop_unflatten(reference_less, "Hybriddetective/Allenwiller_panel_less/Allenwiller_genepop_validation_lociless.txt")
add_pop_line("Hybriddetective/Allenwiller_panel_less/Allenwiller_genepop_validation_lociless.txt")


#### Hybriddetective workflow ####

### for each dataset (full panel, reduced panels):
# 1. simulation of multigenerational hybrids using the training dataset (validation dataset)
# 2. analysis of multigenerational hybrids training dataset with parallelnh
# 3. check convergence
# 4. power quantification and analysis

## Waldi

# simulation of multigenerational hybrids
freqbasedsim_AlleleSample(GPD = "Hybriddetective/Waldi_panel_most/Waldi_genepop_validation_locimost.txt", NumSims = 3, NumReps = 3)
freqbasedsim_AlleleSample(GPD = "Hybriddetective/Waldi_panel_less/Waldi_genepop_validation_lociless.txt", NumSims = 3, NumReps = 3)
freqbasedsim_AlleleSample(GPD = "Hybriddetective/Waldi_panel_all/Waldi_genepop_validation.txt", NumSims = 3, NumReps = 3)

# run newhybrids on each simulation
### change the number of iterations
parallelnh_WIN(folder.data = "Hybriddetective/Waldi_panel_most/",where.NH = path.nh,burnin = 10000, sweeps = 500000) 
parallelnh_WIN(folder.data = "Hybriddetective/Waldi_panel_less/",where.NH = path.nh,burnin = 10000, sweeps = 500000) 
parallelnh_WIN(folder.data = "Hybriddetective/Waldi_panel_all/",where.NH = path.nh,burnin = 10000, sweeps = 500000) 

# delete the first unnecessary folder otherwise it doesn't work
unlink("Hybriddetective/Waldi_panel_most/NH.Results/Waldi_genepop_validation_locimost.txt_Results", recursive = T)
unlink("Hybriddetective/Waldi_panel_less/NH.Results/Waldi_genepop_validation_lociless.txt_Results", recursive = T)
unlink("Hybriddetective/Waldi_panel_all/NH.Results/Waldi_genepop_validation.txt_Results", recursive = T)

# check for convergence of MCMC 
nh_preCheckR(PreDir = "Hybriddetective/Waldi_panel_most/NH.Results/")
nh_preCheckR(PreDir = "Hybriddetective/Waldi_panel_less/NH.Results/") # not converging, trying to increase iterations
nh_preCheckR(PreDir = "Hybriddetective/Waldi_panel_all/NH.Results/")

parallelnh_WIN(folder.data = "Hybriddetective/Waldi_panel_less/",where.NH = path.nh,burnin = 10000, sweeps = 500000) 

#creating plots to check the results
nh_multiplotR(NHResults = "Hybriddetective/Waldi_panel_most/NH.Results/") 
nh_multiplotR(NHResults = "Hybriddetective/Waldi_panel_less/NH.Results/") 
nh_multiplotR(NHResults = "Hybriddetective/Waldi_panel_all/NH.Results/") 

# check efficacy and accuracy of marker panel 
# run the forked function first
my_hybridPowerComp(dir = "Hybriddetective/Waldi_panel_most/NH.Results/") 
my_hybridPowerComp(dir = "Hybriddetective/Waldi_panel_less/NH.Results/") 
my_hybridPowerComp(dir = "Hybriddetective/Waldi_panel_all/NH.Results/") 


## Allenwiller
# simulation of multigenerational hybrids
freqbasedsim_AlleleSample(GPD = "Hybriddetective/Allenwiller_panel_most/Allenwiller_genepop_validation_locimost.txt", NumSims = 3, NumReps = 3)
freqbasedsim_AlleleSample(GPD = "Hybriddetective/Allenwiller_panel_less/Allenwiller_genepop_validation_lociless.txt", NumSims = 3, NumReps = 3)
freqbasedsim_AlleleSample(GPD = "Hybriddetective/Allenwiller_panel_all/Allenwiller_genepop_validation.txt", NumSims = 3, NumReps = 3)

# run newhybrids on each simulation
### change the number of iterations
parallelnh_WIN(folder.data = "Hybriddetective/Allenwiller_panel_most/",where.NH = path.nh,burnin = 10000, sweeps = 500000) 
parallelnh_WIN(folder.data = "Hybriddetective/Allenwiller_panel_less/",where.NH = path.nh,burnin = 10000, sweeps = 500000) # not converging
parallelnh_WIN(folder.data = "Hybriddetective/Allenwiller_panel_all/",where.NH = path.nh,burnin = 10000, sweeps = 500000) 

# delete the first unnecessary folder otherwise it doesn't work
unlink("Hybriddetective/Allenwiller_panel_most/NH.Results/Allenwiller_genepop_validation_locimost.txt_Results", recursive = T)
unlink("Hybriddetective/Allenwiller_panel_less/NH.Results/Allenwiller_genepop_validation_lociless.txt_Results", recursive = T)
unlink("Hybriddetective/Allenwiller_panel_all/NH.Results/Allenwiller_genepop_validation.txt_Results", recursive = T)

# check for convergence of MCMC 
nh_preCheckR(PreDir = "Hybriddetective/Allenwiller_panel_most/NH.Results/")
nh_preCheckR(PreDir = "Hybriddetective/Allenwiller_panel_less/NH.Results/")
nh_preCheckR(PreDir = "Hybriddetective/Allenwiller_panel_all/NH.Results/")

#creating plots to check the results
nh_multiplotR(NHResults = "Hybriddetective/Allenwiller_panel_most/NH.Results/") 
nh_multiplotR(NHResults = "Hybriddetective/Allenwiller_panel_less/NH.Results/") 
nh_multiplotR(NHResults = "Hybriddetective/Allenwiller_panel_all/NH.Results/") 

# run the forked function first
# check efficacy and accuracy of marker panel 
my_hybridPowerComp(dir = "Hybriddetective/Allenwiller_panel_most/NH.Results/") 
my_hybridPowerComp(dir = "Hybriddetective/Allenwiller_panel_less/NH.Results/") # working only after removing the non converged simulations
my_hybridPowerComp(dir = "Hybriddetective/Allenwiller_panel_all/NH.Results/") 


##### plotting efficiency, accuracy and power 

# load all the data files for a specific plot and merge the results
process_files <- function(base_directory, plot_number) {
  # read a plot
  plot_pattern <- paste0("Plot_", plot_number, ".csv$")
  
  # list files for the plot
  plot_files <- list.files(
    path = base_directory,
    pattern = plot_pattern,
    recursive = TRUE,
    full.names = TRUE
  )

  all_data_list <- list()

  for(file_path in plot_files) {
    # find folder containing "panel"
    path_parts <- strsplit(file_path, "/")[[1]]
    folder_name <- path_parts[grep("panel", path_parts)]
    
    # stand from folder
    stand <- ifelse(grepl("Waldi", folder_name), "Waldi", "Allenwiller")
    
    # define panel type
    if(grepl("_all", folder_name)) {
      panel_label <- "all loci"
    } else if(grepl("_less", folder_name)) {
      panel_label <- "8 less informative loci"
    } else if(grepl("_most", folder_name)) {
      panel_label <- "8 most informative loci"
    }

    data <- read.csv(file_path)
    data$stand <- stand
    data$panel <- panel_label

    all_data_list[[length(all_data_list) + 1]] <- data
  }
  
  # combine all
  all_data <- do.call(rbind, all_data_list)
  
  return(all_data)
}

# define the directoy path
base_directory <- paste0(path, "/Hybriddetective")
# process results for the 3 measures specifying plot number
accuracy_res <- process_files(base_directory, 7)
efficiency_res <- process_files(base_directory, 8)
power_res <- process_files(base_directory, 9)

# rename columns 
colnames(accuracy_res)[c(1,2)] <- c("PostProb", "Accuracy")
colnames(efficiency_res)[c(1,2)] <- c("PostProb", "Efficiency")
colnames(power_res)[c(1,2)] <- c("PostProb", "Power")

# save file
write.csv(accuracy_res, "Hybriddetective/Accuracy_results.csv")
write.csv(efficiency_res, "Hybriddetective/Efficiency_results.csv")
write.csv(power_res, "Hybriddetective/Power_results.csv")

# Accuracy Lineplot 
# renaming for the plot
accuracy_res$Genotype_class <- factor(accuracy_res$colour, 
                           levels = c("P1", "P2", "F1", "F2", "BC1", "BC2"),
                           labels = c("Pure sylvatica", "Pure orientalis", "F1", "F2","BC sylvatica", "BC orientalis"))

# reorder panels
accuracy_res$panel <- factor(accuracy_res$panel, levels = c("all loci", "8 most informative loci", "8 less informative loci"))

ColourVector = c( "#FDE725FF","#482173FF", "#25858EFF", "#85D54AFF", "orange","#9C99C8",  "gray8")

# vertical line corresponding to the Posterior Probability threshold chosen for the classification
accuracy_plot <- ggplot(accuracy_res) +
  geom_ribbon(aes(x = PostProb, 
                  ymin = sd.lower, 
                  ymax = sd.upper, 
                  fill = Genotype_class,
                  group = Genotype_class),
              alpha = 0.2) +  # Controls transparency of the shaded area
  geom_line(aes(x = PostProb, y = Accuracy, colour = Genotype_class), lwd = 1.25) +
  geom_vline(xintercept = 0.8, linetype = 2, linewidth = 1) +
  facet_grid(stand ~ panel)   +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        plot.background = element_rect(colour = "white"), 
        panel.grid.major = element_line(colour = "grey90"),
        legend.position="bottom", 
        legend.key = element_blank(),
        strip.background = element_rect(colour = "black", fill = "white")) +
  scale_colour_manual(values = ColourVector) +
  scale_fill_manual(values = ColourVector) +  # Match the fill colors to the line colors
  labs(x = "Posterior Probability Threshold", 
       y = expression("Accuracy "%+-%"sd"), 
       col = "Genotype class",
       fill = "Genotype class") +  # Add label for the fill legend
  coord_cartesian(ylim = c(0, 1))

accuracy_plot

path.fig = "C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Manuscript/Figures/"
png(filename= paste0(path.fig, "accuracy_plot.png"), width =10000,height=5000, res=800)
plot(accuracy_plot)
dev.off()

# Efficiency line plot
# renaming for the plot
efficiency_res$Genotype_class <- factor(efficiency_res$colour, 
                                      levels = c("Pure1", "Pure2", "F1", "F2", "BC1", "BC2"),
                                      labels = c("Pure sylvatica", "Pure orientalis", "F1", "F2","BC sylvatica", "BC orientalis"))

# reorder panels
efficiency_res$panel <- factor(efficiency_res$panel, levels = c("all loci", "8 most informative loci", "8 less informative loci"))

ColourVector = c( "#FDE725FF","#482173FF", "#25858EFF", "#85D54AFF", "orange","#9C99C8",  "gray8")

# vertical line corresponding to the Posterior Probability threshold chosen for the classification

efficiency_plot <- ggplot(efficiency_res) +
  geom_ribbon(aes(x = PostProb, 
                  ymin = sd.lower, 
                  ymax = sd.upper, 
                  fill = Genotype_class,
                  group = Genotype_class),
              alpha = 0.2) +  # Controls transparency of the shaded area
  geom_line(aes(x = PostProb, y = Efficiency, colour = Genotype_class), lwd = 1.25) +
  geom_vline(xintercept = 0.8, linetype = 2, linewidth = 1) +
  facet_grid(stand ~ panel)   +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        plot.background = element_rect(colour = "white"), 
        panel.grid.major = element_line(colour = "grey90"),
        legend.position="bottom", 
        legend.key = element_blank(),
        strip.background = element_rect(colour = "black", fill = "white")) +
  scale_colour_manual(values = ColourVector) +
  scale_fill_manual(values = ColourVector) +  # Match the fill colors to the line colors
  labs(x = "Posterior Probability Threshold", 
       y = expression("Accuracy "%+-%"sd"), 
       col = "Genotype class",
       fill = "Genotype class") +  # Add label for the fill legend
  coord_cartesian(ylim = c(0, 1))

efficiency_plot


path.fig = "C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Manuscript/Figures/"
png(filename= paste0(path.fig, "efficiency_plot.png"), width =10000,height=5000, res=800)
plot(efficiency_plot)
dev.off()

# power line plot
# renaming for the plot
power_res$Genotype_class <- factor(power_res$colour, 
                                        levels = c("P1", "P2", "F1", "F2", "BC1", "BC2"),
                                        labels = c("Pure sylvatica", "Pure orientalis", "F1", "F2","BC sylvatica", "BC orientalis"))

# reorder panels
power_res$panel <- factor(power_res$panel, levels = c("all loci", "8 most informative loci", "8 less informative loci"))

ColourVector = c( "#FDE725FF","#482173FF", "#25858EFF", "#85D54AFF", "orange","#9C99C8",  "gray8")

# vertical line corresponding to the Posterior Probability threshold chosen for the classification

power_plot <- ggplot(power_res) +
  geom_ribbon(aes(x = PostProb, 
                  ymin = sd.lower, 
                  ymax = sd.upper, 
                  fill = Genotype_class,
                  group = Genotype_class),
              alpha = 0.2) +  # Controls transparency of the shaded area
  geom_line(aes(x = PostProb, y = Power, colour = Genotype_class), lwd = 1.25) +
  geom_vline(xintercept = 0.8, linetype = 2, linewidth = 1) +
  facet_grid(stand ~ panel)   +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        plot.background = element_rect(colour = "white"), 
        panel.grid.major = element_line(colour = "grey90"),
        legend.position="bottom", 
        legend.key = element_blank(),
        strip.background = element_rect(colour = "black", fill = "white")) +
  scale_colour_manual(values = ColourVector) +
  scale_fill_manual(values = ColourVector) +  # Match the fill colors to the line colors
  labs(x = "Posterior Probability Threshold", 
       y = expression("Power "%+-%"sd"), 
       col = "Genotype class",
       fill = "Genotype class") +  # Add label for the fill legend
  coord_cartesian(ylim = c(0, 1))

power_plot

path.fig = "C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Manuscript/Figures/"
png(filename= paste0(path.fig, "power_plot.png"), width =10000,height=5000, res=800)
plot(power_plot)
dev.off()
