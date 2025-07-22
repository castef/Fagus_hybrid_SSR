#### creating input files for NMpi2 ####

library("xlsx")
library("dplyr")
library("tidyverse")
library("graph4lg")
library("sf")
library("sp")
library("raster")
library("ggplot2")
library("tidyverse")
library("terra")

setwd(dir="~/Data_clean")

## each file must contain:
# generation code (0 = adult, 1, offspring)
# ID number (seq 1: nindividuals, 2 sequences per generation)
# X coordinate
# Y coordinate
# cyDNA (-1 if missing)
# alleles column (one column per allele)
# (N standardized covariate for adults OR femaleness), mother ID for offpsring (-1 if unknown)
# (N standardized covariate for adults OR femaleness), father ID for offspring (-1 if unkonw)
# same for adults, dispersal index for offspring (1 if dispersed)

# Allenwiller

# load individuals data
data <- read.csv("Allenwiller_genotypes_withcoord.csv")

# work on adults data
adults <- subset(data, generation == "Adult")

# remove adults with Na coordinates
adults <- adults[!is.na(adults$x),]
adults <- adults[!is.na(adults$y),]

#create basic datafrme
input <- data.frame(gen = 0, 
                       seq = 1:nrow(adults), 
                       x = adults$x, 
                       y = adults$y, 
                       cyDNA = -1)
# add the genotype
input <- cbind(input, adults[, c(5:36)])

# add the femaleness
input$fem = 0.5

# OPTIONAL: add the covariates
input2 <- input
# DBH
input2$z_DBH = scale(as.numeric(adults$circumference)/pi)

# posterior probability of sylvatica from NH
nh_res <- read.csv("Newhybrids/Allenwiller_parallelnh/Allenwiller_parallelnh_results.csv")
nh_res <- merge(adults, nh_res[, c(2,4)], by.x = "SampleID", by.y = "IndivName", all.x = TRUE)
input2$z_psyl = scale(nh_res$X0.000.0.000.0.000.1.000)

# add another column if you want to test the assortative mating (psyl = biparental factor)
input2$z_psyl_bip = scale(nh_res$X0.000.0.000.0.000.1.000)

# add the competition coefficient

# calculate competition coefficient first (formula from De Sauvage 2023 - Hagyi 1974)

# for each tree, compute the euclidean distance with all the other trees
eucl_trees = graph4lg::mat_geo_dist(adults, ID = "SampleID", x = "x", y = "y")
hist(eucl_trees)

# radius of competition 
radius = 20
# array to store the competition coefficient
adults$DBH = as.numeric(adults$circumference)/pi
comp_coeff = numeric(length(adults$DBH))

# iterate over each tree the competition coefficient
for (i in 1:length(adults$DBH)) {
  competition_sum = 0
  # for each tree i, iterate over all trees j
  for (j in 1:length(adults$DBH)) {
    # if tree j is within 20 meters and is not the same as tree i, compute the competition component
    if (i != j && eucl_trees[i, j] <= radius) {
      # sum all the components to get the competition coefficient for the tree i
      competition_sum = competition_sum + (as.numeric(adults$DBH)[j] / as.numeric(adults$DBH)[i]) / eucl_trees[i, j]
    }
  }
  comp_coeff[i] = competition_sum
}

hist(comp_coeff)

# add as covairate
input2$z_comp_coeff = scale(comp_coeff)

# check if there is any collinerity between the covariates!
cor(input2[ c(39:42)])
plot(input2[ c(39:42)])

# reloacte all the covariates before femaleness
input2 = input2 %>%   relocate(starts_with("z_"), .before = fem)

### offsprings AND JUVENILES!
offspring <- subset(data, generation == "Offspring" | generation == "Juvenile")
# remove offpsring with Na coordinates
offspring <- offspring[!is.na(offspring$x),]
offspring <- offspring[!is.na(offspring$y),]

#create basic datafrme
input3 <- data.frame(gen = 1, 
                    seq = 1:nrow(offspring), 
                    x = offspring$x, 
                    y = offspring$y, 
                    cyDNA = -1)
# add the genotype
input3 <- cbind(input3, offspring[, c(5:36)])
# add mother, father and dispersed type
input3$mother = -1
input3$father = -1
input3$disp = 1

# chec n columns
ncol(input2)
ncol(input3)

# save the data for merging later with Waldi
adults_a = input2
offspring_a = input3

# add header (number of adults, number of offspring, number of loci, number of covariates)
header <- c(nrow(input2), nrow(input3), 16, sum(grepl("^z_", colnames(input2))))

# write header and after adults and offspring dataframes 
output_file <- "Mating_model/Input/Allenwiller_NMpi2_input_DBH_Psyl_PsylB_compcoeff.txt"

writeLines(paste(header, collapse = " "), con = output_file)
write.table(input2, output_file, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

# since offspring have 2 column less, need to add empty column first
input3$empty <- ""
input3$empty2 <- ""
write.table(input3, output_file, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)




# Waldi

# load individuals data
data <- read.csv("Waldi_genotypes_withcoord.csv")

# work on adults data
adults <- subset(data, generation == "Adult")

# remove adults with Na coordinates
adults <- adults[!is.na(adults$x),]
adults <- adults[!is.na(adults$y),]

#create basic datafrme
input <- data.frame(gen = 0, 
                    seq = 1:nrow(adults), 
                    x = adults$x, 
                    y = adults$y, 
                    cyDNA = -1)
# add the genotype
input <- cbind(input, adults[, c(5:36)])

# add the femaleness
input$fem = 0.5

# OPTIONAL: add the covariates
input2 <- input
# DBH
input2$z_DBH = scale(as.numeric(adults$circumference)/pi)

# posterior probability of sylvatica from NH
nh_res <- read.csv("Newhybrids/Waldi_parallelnh/Waldi_parallelnh_results.csv")
nh_res <- merge(adults, nh_res[, c(2,4)], by.x = "SampleID", by.y = "IndivName", all.x = TRUE)
input2$z_psyl = scale(nh_res$X0.000.0.000.0.000.1.000)

# add another column if you want to test the assortative mating (psyl = biparental factor)
input2$z_psyl_bip = scale(nh_res$X0.000.0.000.0.000.1.000)

# add the competition coefficient

# calculate competition coefficient first (formula from De Sauvage 2023 - Hagyi 1974)

# for each tree, compute the euclidean distance with all the other trees
eucl_trees = graph4lg::mat_geo_dist(adults, ID = "SampleID", x = "x", y = "y")
hist(eucl_trees)

# radius of competition 
radius = 20
# array to store the competition coefficient
adults$DBH = as.numeric(adults$circumference)/pi
comp_coeff = numeric(length(adults$DBH))

# iterate over each tree the competition coefficient
for (i in 1:length(adults$DBH)) {
  competition_sum = 0
  # for each tree i, iterate over all trees j
  for (j in 1:length(adults$DBH)) {
    # if tree j is within 20 meters and is not the same as tree i, compute the competition component
    if (i != j && eucl_trees[i, j] <= radius) {
      # sum all the components to get the competition coefficient for the tree i
      competition_sum = competition_sum + (as.numeric(adults$DBH)[j] / as.numeric(adults$DBH)[i]) / eucl_trees[i, j]
    }
  }
  comp_coeff[i] = competition_sum
}

hist(comp_coeff)

# add as covairate
input2$z_comp_coeff = scale(comp_coeff)

# check if there is any collinerity between the covariates!
cor(input2[ c(39:42)])
plot(input2[ c(39:42)])

# reloacte all the covariates before femaleness
input2 = input2 %>%   relocate(starts_with("z_"), .before = fem)

### offsprings
offspring <- subset(data, generation == "Offspring" | generation == "Juvenile")
# remove offpsring with Na coordinates
offspring <- offspring[!is.na(offspring$x),]
offspring <- offspring[!is.na(offspring$y),]

#create basic datafrme
input3 <- data.frame(gen = 1, 
                     seq = 1:nrow(offspring), 
                     x = offspring$x, 
                     y = offspring$y, 
                     cyDNA = -1)
# add the genotype
input3 <- cbind(input3, offspring[, c(5:36)])
# add mother, father and dispersed type
input3$mother = -1
input3$father = -1
input3$disp = 1


ncol(input2)
ncol(input3)
# again saving waldi to merge with allenwiller
adults_w = input2
offspring_w = input3

# add header (number of adults, number of offspring, number of loci, number of covariates)
header <- c(nrow(input2), nrow(input3), 16, sum(grepl("^z_", colnames(input2))))

# write header and after adults and offspring dataframes 
output_file <- "Mating_model/Input/Waldi_NMpi2_input_DBH_Psyl_PsylB_compcoeff.txt"

writeLines(paste(header, collapse = " "), con = output_file)
write.table(input2, output_file, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

# since offspring have 2 column less, need to add empty column first
input3$empty <- ""
input3$empty2 <- ""
write.table(input3, output_file, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)



### writing the combined input file
colnames(adults_a) == colnames(adults_w)
colnames(adults_a) = colnames(adults_w)
adults_all = rbind(adults_a, adults_w)

# adjust seq column
adults_all$seq = 1:nrow(adults_all)

# remove phenotypic characters
adults_all = adults_all[, -c(38:41)]

# offspring
colnames(offspring_a) == colnames(offspring_w)
colnames(offspring_a) = colnames(offspring_w)
offspring_all = rbind(offspring_a, offspring_w)
offspring_all$seq = 1:nrow(offspring_all)

# adjust header
header <- c(nrow(adults_all), nrow(offspring_all), 16, sum(grepl("^z_", colnames(adults_all))))

ncol(adults_all)
ncol(offspring_all)

# since adults have 3 column less, need to add empty column first
adults_all$empty <- ""
adults_all$empty2 <- ""

# write header and after adults and offspring dataframes 
output_file <- "Mating_model/Input/Sites_merged_NMpi2.txt"

writeLines(paste(header, collapse = " "), con = output_file)
write.table(adults_all, output_file, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(offspring_all, output_file, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
