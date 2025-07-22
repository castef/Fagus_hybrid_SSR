# genotype classes with NewHybrids # 

# requires rtools to install old packages

install.packages("remotes")
library(remotes)
install.packages("devtools")
library("devtools")

# if the commands don't work try copying the package folder in the Program Files
devtools::install_github("bwringe/parallelnewhybrid") #required
devtools::install_github("kkeenan02/diveRsity", dependencies = TRUE)
devtools::install_github("rystanley/genepopedit", dependencies = TRUE) #required to convert to NH format
devtools::install_github("bwringe/hybriddetective", dependencies = TRUE) #The original package has a bug, install the one below
devtools::install_github("stevemussmann/hybriddetective", force = TRUE) ## fixed version (found here https://github.com/bwringe/hybriddetective/issues/12 )

install.packages("qgraph")
install.packages("xlsx")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("grid")
install.packages("parallel")
install.packages("plyr")
install.packages("stringr")
install.packages("scales")
install.packages("reshape2")
install.packages("tidyr")
install.packages("adegenet")

library("xlsx")
library("qgraph")
library("plyr")
library("dplyr")
library("ggplot2")
library("grid")
library("stringr")
library("scales")
library("reshape2")
library("tidyr")
library("parallel")
library("adegenet")
library("diveRsity")
library("genepopedit") ## needed for converting files
library("parallelnewhybrid")
library("hybriddetective")

setwd(dir="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Data_clean")
path = "C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Data_clean"

#----  Waldi dataset preparation ----

## writing files for NewHybrid analysis on all the individuals
df_val <- read.xlsx("NewHybrids/Waldi_all.xlsx", sheetIndex = 1)
colnames(df_val)

# changing all alleles under 100 to standard 3-digits format
for (j in 5:36) {
  for (i in 1:nrow(df_val)){
    value <- as.numeric(df_val[i,j])
    if (value < 100){
      df_val[i,j]= as.character(paste0(0,value))
    } else {
      next
    }
  }
}
#merging loci columns
df_val$casolfagus_29 <- paste0(df_val$casolfagus_29,df_val$NA..1)
df_val$concat14 <- paste0(df_val$concat14,df_val$NA..2)
df_val$csolfagus_05 <- paste0(df_val$csolfagus_05,df_val$NA..3)
df_val$csolfagus_06 <- paste0(df_val$csolfagus_06,df_val$NA..4)
df_val$csolfagus_19 <- paste0(df_val$csolfagus_19,df_val$NA..5)
df_val$csolfagus_31 <- paste0(df_val$csolfagus_31, df_val$NA..6)
df_val$DE576 <- paste0(df_val$DE576,df_val$NA..7)
df_val$DUKCT <- paste0(df_val$DUKCT,df_val$NA..8) 
df_val$DZ447 <- paste0(df_val$DZ447,df_val$NA..9)
df_val$EEU75 <- paste0(df_val$EEU75,df_val$NA..10)
df_val$EJV8T <- paste0(df_val$EJV8T,df_val$NA..11)
df_val$EMILY <- paste0(df_val$EMILY,df_val$NA..12)
df_val$ERHBI <- paste0(df_val$ERHBI, df_val$NA..13)
df_val$FS1_15 <- paste0(df_val$FS1_15,df_val$NA..14)
df_val$sfc_0036 <- paste0(df_val$sfc_0036, df_val$NA..15)
df_val$sfc_1143 <- paste0(df_val$sfc_1143, df_val$NA..16)
df_val <- df_val[, -grep("NA.", colnames(df_val))]
df_val[df_val=="0-10-1"] <- as.character("000000") #missing data

# Genepop file of all the individuals
genepopedit::genepop_unflatten(df_val[, -c(2,3)], "NewHybrids/Waldi_genepop_all.txt")

# write inviduals files 
write.table(df_val$SampleID, "NewHybrids/Waldi_all_individuals.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# NH format 
genepopedit::genepop_newhybrids(genepop="NewHybrids/Waldi_genepop_all.txt",path = "Newhybrids/Waldi_NH_all.txt")



#---- Allenwiller dataset preparation  ----  
# all individuals for the final run parallelnewhybrid
df_val <- read.csv("NewHybrids/Allenwiller_all.csv")

# write individuals file
write.table(df_val$SampleID, "NewHybrids/Allenwiller_all_individuals.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

df_val <- df_val[,-c(1,3:6)]
for (j in 2:33) {
  for (i in 1:nrow(df_val)){
    value <- as.numeric(df_val[i,j])
    if (value < 100){
      df_val[i,j]= as.character(paste0(0,value))
    } else {
      next
    }
  }
}
#merging loci columns
df_val$casolfagus_29 <- paste0(df_val$casolfagus_29,df_val$X.1)
df_val$concat14 <- paste0(df_val$concat14,df_val$X.2)
df_val$csolfagus_05 <- paste0(df_val$csolfagus_05,df_val$X.3)
df_val$csolfagus_06 <- paste0(df_val$csolfagus_06,df_val$X.4)
df_val$csolfagus_19 <- paste0(df_val$csolfagus_19,df_val$X.5)
df_val$csolfagus_31 <- paste0(df_val$csolfagus_31, df_val$X.6)
df_val$DE576 <- paste0(df_val$DE576,df_val$X.7)
df_val$DUKCT <- paste0(df_val$DUKCT,df_val$X.8)
df_val$DZ447 <- paste0(df_val$DZ447,df_val$X.9)
df_val$EEU75 <- paste0(df_val$EEU75,df_val$X.10)
df_val$EJV8T <- paste0(df_val$EJV8T,df_val$X.11)
df_val$EMILY <- paste0(df_val$EMILY,df_val$X.12)
df_val$ERHBI <- paste0(df_val$ERHBI, df_val$X.13)
df_val$FS1_15 <- paste0(df_val$FS1_15,df_val$X.14)
df_val$sfc_0036 <- paste0(df_val$sfc_0036, df_val$X.15)
df_val$sfc_1143 <- paste0(df_val$sfc_1143, df_val$X.16)
df_val <- df_val[, -grep("X", colnames(df_val))]
df_val[df_val=="0-10-1"] <- as.character("000000") #missing data

# writing genepop format with all the individuals (removing generation column)
genepopedit::genepop_unflatten(df_val,"NewHybrids/Allenwiller_genepop_all.txt" )

# Newhybrids file all individuals 
genepopedit::genepop_newhybrids(genepop="NewHybrids/Allenwiller_genepop_all.txt",path = "NewHybrids/Allenwiller_NH_all.txt")  

#---- running Newhybrids (parallelnewhybrids) on the whole dataset ----
path.hold <- getwd()
# set the path of newhybrids
path.nh <- "C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/NewHybrids/packages/newhybrids-master/bin/PC/newhybrids/"


# Waldi # 
dir.create("NewHybrids/Waldi_parallelnh")

#copy file NH format and individuals file to the directory for the analysis
file.copy(from = c("NewHybrids/Waldi_NH_all.txt", "NewHybrids/Waldi_all_individuals.txt"),
          to = c("NewHybrids/Waldi_parallelnh/Waldi_NH_all.txt", "NewHybrids/Waldi_parallelnh/Waldi_all_individuals.txt"))

#running NewHybrids on all data
parallelnh_WIN(folder.data ="NewHybrids/Waldi_parallelnh/",where.NH = path.nh, burnin = 10000, sweeps = 50000) # sweep = iteration of MCMC algorithm

# check plot
nh_plotR(NHResults = "NewHybrids/Waldi_parallelnh/NH.Results/Waldi_NH_all.txt_Results/Waldi_NH_all.txt_PofZ.txt", ReversePure = 2)

# replacing names in results file and saving results (removing first useless row)
nh_res <-  read.table("NewHybrids/Waldi_parallelnh/NH.Results/Waldi_NH_all.txt_Results/Waldi_NH_all.txt_PofZ.txt")
nh_res$V2[-1] <- as.data.frame(read.table("NewHybrids/Waldi_parallelnh/Waldi_all_individuals.txt"))$V1
write.table(nh_res, "NewHybrids/Waldi_parallelnh/Waldi_parallelnh_results.csv", sep=",",row.names = F, col.names=F)

### Allenwiller 
dir.create("NewHybrids/Allenwiller_parallelnh")

#copy file NH format and individuals file to the directory for the analysis
file.copy(from = c("NewHybrids/Allenwiller_NH_all.txt", "NewHybrids/Allenwiller_all_individuals.txt"),
          to = c("NewHybrids/Allenwiller_parallelnh/Allenwiller_NH_all.txt", "NewHybrids/Allenwiller_parallelnh/Allenwiller_all_individuals.txt"))

#running NewHybrids
parallelnh_WIN(folder.data ="NewHybrids/Allenwiller_parallelnh/",where.NH = path.nh, burnin = 10000, sweeps = 500000) # sweep = iteration of MCMC algorithm

# check plot (first individuals are orientalis, therefore colors are reversed to keep consistency with Waldi)
# Pure 1 = red = sylvatica, Pure 2 = blue = orientalis
nh_plotR(NHResults = "NewHybrids/Allenwiller_parallelnh/NH.Results/Allenwiller_NH_all.txt_Results/Allenwiller_NH_all.txt_PofZ.txt", 
         ReversePure = 2)

# replacing names in results file and saving results (removing first useless row)
nh_res <- read.table("NewHybrids/Allenwiller_parallelnh/NH.Results/Allenwiller_NH_all.txt_Results/Allenwiller_NH_all.txt_PofZ.txt")
nh_res$V2[-1] <- as.data.frame(read.table("NewHybrids/Allenwiller_parallelnh/Allenwiller_all_individuals.txt"))$V1
write.csv(nh_res, "NewHybrids/Allenwiller_parallelnh/Allenwiller_parallelnh_results.csv", sep=",",row.names = F, col.names=F)

#---- merge with tree info ----

### Allenwiller
# importing output newhybrid and individuals ID and merging
nh_output = read.table("NewHybrids/Allenwiller_parallelnh/NH.Results/Allenwiller_NH_all.txt_Results/Allenwiller_NH_all.txt_PofZ.txt", header = TRUE)[,-1]
Indiv_name = read.table("NewHybrids/Allenwiller_parallelnh/NH.Results/Allenwiller_NH_all.txt_Results/Allenwiller_all_individuals.txt")
nh_output$IndivName = as.factor(Indiv_name$V1)
colnames(nh_output) = c("IndivName", "Pure2", "Pure1","F1", "F2", "BC2", "BC1") # assigning names with Pure1 = sylvatica, Pure2 = orientalis (reversed!!!)

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

# add information of generation 
generation = read.csv("Allenwiller_indiv_info.csv")
NH_melt = merge(NH_melt, generation[, 1:2], by.x = "IndivName", by.y = "SampleID")

# renaming for the plot
NH_melt$PopProb2 <- factor(NH_melt$PopProb, 
  levels = c("Pure1", "Pure2", "F1", "F2", "BC1", "BC2"),
  labels = c("Pure sylvatica", "Pure orientalis", "F1", "F2","BC sylvatica", "BC orientalis"))

# plot all the individuals
ColourVector = c( "#FDE725FF","#482173FF", "#25858EFF", "#85D54AFF", "orange","#9C99C8",  "gray8")

al_nh_plot = ggplot(NH_melt, aes(x = IndivName, y=as.numeric(CumProb), fill=PopProb2))+ # CumProb changed to numeric
  geom_bar(stat="identity", position = "stack", width = 1) + 
  scale_fill_manual(values=ColourVector)+
  scale_colour_manual(values=ColourVector)+ylab("Cumulative Probability")+xlab("Individual") +
  scale_y_continuous(limits = c(0, 1.05), expand=c(0, 0)) + 
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black") + # Add dotted horizontal line
  theme(axis.title.x = element_text(size = 10), axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 10),
        panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"),
        legend.position="top", strip.background = element_blank()) + guides(fill=guide_legend(nrow=1, title = "")) +
  facet_wrap(~generation, drop = T, scales = "free", ncol = 1)
al_nh_plot

### Waldi
# importing output newhybrid and individuals ID and merging
nh_output = read.table("NewHybrids/Waldi_parallelnh/NH.Results/Waldi_NH_all.txt_Results/Waldi_NH_all.txt_PofZ.txt", header = TRUE)[,-1]
Indiv_name = read.table("NewHybrids/Waldi_parallelnh/NH.Results/Waldi_NH_all.txt_Results/Waldi_all_individuals.txt")
nh_output$IndivName = as.factor(Indiv_name$V1)
colnames(nh_output) = c("IndivName", "Pure2", "Pure1","F1", "F2", "BC2", "BC1") # assigning names with Pure1 = sylvatica, Pure2 = orientalis (reversed!!!)

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

# add information of generation 
generation = read.csv("Waldi_indiv_info.csv")
NH_melt = merge(NH_melt, generation[, 1:3], by.x = "IndivName", by.y = "SampleID")

# renaming for the plot
NH_melt$PopProb2 <- factor(NH_melt$PopProb, 
                           levels = c("Pure1", "Pure2", "F1", "F2", "BC1", "BC2"),
                           labels = c("Pure sylvatica", "Pure orientalis", "F1", "F2","BC sylvatica", "BC orientalis"))

# plot all the individuals
ColourVector = c( "#FDE725FF","#482173FF", "#25858EFF", "#85D54AFF", "orange","#9C99C8",  "gray8")

wal_nh_plot = ggplot(NH_melt, aes(x = IndivName, y=as.numeric(CumProb), fill=PopProb2))+ # CumProb changed to numeric
  geom_bar(stat="identity", position = "stack", width = 1) + 
  scale_fill_manual(values=ColourVector)+
  scale_colour_manual(values=ColourVector)+ylab("Cumulative Probability")+xlab("Individual") +
  scale_y_continuous(limits = c(0, 1.05), expand=c(0, 0)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black") + # Add dotted horizontal line
  theme(axis.title.x = element_text(size = 10), axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 10),
        panel.background = element_rect(fill = "white", colour = "black"), 
        plot.background = element_rect(colour = "white"),
        legend.position="top", strip.background = element_blank()) + 
  guides(fill=guide_legend(nrow=1, title = "")) +
  facet_wrap(~generation, drop = T, scales = "free", ncol = 1)
wal_nh_plot

# plotting all together and saving (need to align or save separately and do with inkscape)
fig_nh_plot=ggpubr::ggarrange(al_nh_plot, wal_nh_plot,
                              common.legend = T,
                              #align = "h", 
                              labels = c("Allenwiller", "Waldi"))
fig_nh_plot

png(filename="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Manuscript/Figures/Nh_plot_last.png", width =19000,height=8000, res=1500)
plot(fig_nh_plot)
dev.off()

# save singular plots
png(filename="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Manuscript/Figures/All_Nh_plot_last.png", width =10000,height=9000, res=1500)
plot(al_nh_plot)
dev.off()
png(filename="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Manuscript/Figures/Wal_Nh_plot_last.png", width =10000,height=7000, res=1500)
plot(wal_nh_plot)
dev.off()

#### merging results with dataset ###

## CLASSIFICATION: 
# Pure1 = Pure sylvatica
# Pure2 = Pure orientalis
# F1 = Pure sylvatica x Pure orientalis
# F2 = F1 x F1
# BC1 = Pure sylvatica x F1
# BC2 = Pure orientalis x F1

### Waldi ###
# loading nh results
nh_res <- read.csv("NewHybrids/Waldi_parallelnh/Waldi_parallelnh_results.csv")

# NB: original classification for Allenwiller: pure sylvatica = Pure2, pure orientalis = Pure 1 --> renaming columns so that pure sylvatica=Pure1 
colnames(nh_res) <- c("500000_sweeps", "IndivName", "Pure2", "Pure1", "F1", "F2", "BC2", "BC1")

#finding the class with the higher posterior probability
nh_res$NH_max_class <- colnames(nh_res[, 3:8])[max.col(nh_res[, 3:8], ties.method = "first")]
# getting the posterior prob of the max class
nh_res <- nh_res %>%
  rowwise() %>%  mutate(max_class_pp = max(Pure2, Pure1, F1, F2, BC2, BC1))

# classifying as assigned or unassigned to hybrid class depending to the critical pp threshold
nh_res$t75 <- ifelse(nh_res$max_class_pp > 0.75, "assigned", "unassigned")
nh_res$t80 <- ifelse(nh_res$max_class_pp > 0.8, "assigned", "unassigned")
nh_res$t85 <- ifelse(nh_res$max_class_pp > 0.85, "assigned", "unassigned")
nh_res$t90 <- ifelse(nh_res$max_class_pp > 0.9, "assigned", "unassigned")

table(nh_res$t75)
table(nh_res$t80)
table(nh_res$t85)
table(nh_res$t90)

# classification according to critical posterior probability threshold = 0.8
nh_res$NH_class_t80 <- ifelse(nh_res$t80 == "assigned", nh_res$NH_max_class, nh_res$t80)

# merging classification results with tree info file
indiv_info <- read.csv("Waldi_indiv_info_withcoord.csv")
indiv_class <- merge(indiv_info, nh_res[,c("IndivName", "NH_max_class", "NH_class_t80")], by.x = "SampleID", by.y = "IndivName", all.x = TRUE, all.y = TRUE)

# check dimensions
dim(indiv_info)
dim(indiv_class)

# renaming classes
indiv_class$NH_class_t80_corrected <- factor(indiv_class$NH_class_t80, 
                                   levels = c("Pure1", "Pure2", "F1", "F2", "BC1", "BC2", "unassigned"),
                                   labels = c("Pure sylvatica", "Pure orientalis", "F1", "F2","BC sylvatica", "BC orientalis", "unassigned"))

# check how many individuals per class
table(indiv_class$Ageclass_name, indiv_class$NH_class_t80_corrected)

# check the adult BC orientalis
indiv_class[which(indiv_class$NH_class_t80_corrected == "BC orientalis"),]
nh_res[which(nh_res$NH_class_t80 == "BC2"),]

# correct and classify as unassigned
indiv_class[which(indiv_class$NH_class_t80_corrected == "BC orientalis"),"NH_class_t80_corrected"] = "unassigned"
table(indiv_class$Ageclass_name, indiv_class$NH_class_t80_corrected)

# tidy and write 
indiv_class = indiv_class[, -11]
write.csv(indiv_class, "NewHybrids/Waldi_indiv_info_class.csv")

### Allenwiller ###
# loading nh results
nh_res <- read.csv("NewHybrids/Allenwiller_parallelnh/Allenwiller_parallelnh_results.csv")

# NB: original classification for Allenwiller: pure sylvatica = Pure2, pure orientalis = Pure 1 --> renaming columns so that pure sylvatica=Pure1 
colnames(nh_res) <- c("500000_sweeps", "IndivName", "Pure2", "Pure1", "F1", "F2", "BC2", "BC1")

#finding the class with the higher posterior probability
nh_res$NH_max_class <- colnames(nh_res[, 3:8])[max.col(nh_res[, 3:8], ties.method = "first")]
# getting the posterior prob of the max class
nh_res <- nh_res %>%
  rowwise() %>%  mutate(max_class_pp = max(Pure2, Pure1, F1, F2, BC2, BC1))

# classifying as assigned or unassigned to hybrid class depending to the critical pp threshold
nh_res$t75 <- ifelse(nh_res$max_class_pp > 0.75, "assigned", "unassigned")
nh_res$t80 <- ifelse(nh_res$max_class_pp > 0.8, "assigned", "unassigned")
nh_res$t85 <- ifelse(nh_res$max_class_pp > 0.85, "assigned", "unassigned")
nh_res$t90 <- ifelse(nh_res$max_class_pp > 0.9, "assigned", "unassigned")

table(nh_res$t75)
table(nh_res$t80)
table(nh_res$t85)
table(nh_res$t90)

# classification according to critical posterior probability threshold = 0.8
nh_res$NH_class_t80 <- ifelse(nh_res$t80 == "assigned", nh_res$NH_max_class, nh_res$t80)

# merging classification results with tree info file
indiv_info <- read.csv("Allenwiller_indiv_info_withcoord.csv")
indiv_class <- merge(indiv_info, nh_res[,c("IndivName", "NH_max_class", "NH_class_t80")], by.x = "SampleID", by.y = "IndivName", all.x = TRUE, all.y = TRUE)

# check dimensions to make sure don't miss individuals
dim(indiv_info)
dim(indiv_class)

# renaming classes
indiv_class$NH_class_t80_2 <- factor(indiv_class$NH_class_t80, 
                                     levels = c("Pure1", "Pure2", "F1", "F2", "BC1", "BC2", "unassigned"),
                                     labels = c("Pure sylvatica", "Pure orientalis", "F1", "F2","BC sylvatica", "BC orientalis", "unassigned"))

# check how many individuals per class
table(indiv_class$generation, indiv_class$NH_class_t80_2)

# check BC
indiv_class[which(indiv_class$NH_class_t80_2 == "BC sylvatica"),]
nh_res[which(nh_res$NH_class_t80 == "BC1"),]

# tidy and write 
indiv_class = indiv_class[, -11]
colnames(indiv_class)[11] <- "NH_class_t80"
#write.csv(indiv_class, "NewHybrids/Allenwiller_indiv_info_class.csv")

# subset
subset <- fread("Allenwiller_indiv_info_withcoord_subset.csv")
indiv_class_subset <- indiv_class[indiv_class$SampleID %in% subset$SampleID, ]
#write.csv(indiv_class_subset, "NewHybrids/Allenwiller_indiv_info_class_subset.csv")

#---- calculating proportions of classes ----
## full dataset
obs1 <- fread("NewHybrids/Allenwiller_indiv_info_class.csv")
obs1$stand = "Allenwiller"
obs2 <- fread("NewHybrids/Waldi_indiv_info_class.csv")
obs2$stand <- "Waldi"

# keep the correct classification column
obs2$NH_class_t80 <- obs2$NH_class_t80_corrected
obs2$NH_class_t80_corrected <- NULL
# reorder columns and merge
obs1$circumference <- as.integer(obs1$circumference)
obs <- bind_rows(obs1, obs2)

# chekc how many unassigned and assigned per stand and age class
obs %>%
  group_by(stand, generation) %>%
  summarize(
    total = n(),
    unassigned = sum(NH_class_t80 == "unassigned", na.rm = TRUE),
    prop_unassigned = unassigned / total
  )
  
# calculating proportions
obs_proportion <- obs %>%
  group_by(stand, Ageclass_name, NH_class_t80) %>%
  summarize(genot_n = n(), .groups = "drop") %>% 
  group_by(stand, Ageclass_name) %>% # total number of individuals per stand and stage_group
  mutate(stage_tot = sum(genot_n),
         genot_prop = genot_n / stage_tot) %>% # proportion of each genotype within the stage group
  ungroup()
obs_proportion

# fwrite(obs_proportion, "Observed_proportions.csv")

## considering subset Allenwiller
obs1 <- fread("NewHybrids/Allenwiller_indiv_info_class_subset.csv")
obs1$stand = "Allenwiller"
obs2 <- fread("NewHybrids/Waldi_indiv_info_class.csv")
obs2$stand <- "Waldi"
# keep the correct classification column
obs2$NH_class_t80 <- obs2$NH_class_t80_corrected
obs2$NH_class_t80_corrected <- NULL
# reorder columns and merge
obs1$circumference <- as.integer(obs1$circumference)
obs <- bind_rows(obs1, obs2)

# calculating proportions
obs_proportion <- obs %>%
  group_by(stand, Ageclass_name, NH_class_t80) %>%
  summarize(genot_n = n(), .groups = "drop") %>% 
  group_by(stand, Ageclass_name) %>% # total number of individuals per stand and age class
  mutate(stage_tot = sum(genot_n),
         genot_prop = genot_n / stage_tot) %>% # proportion of each genotype within the stage group
  ungroup()
obs_proportion

# fwrite(obs_proportion, "Observed_proportions_subset.csv")
