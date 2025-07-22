### Dataset preparation ###

library(dplyr)
#library(gsl)
library(raster)
library(terra)
library(ggplot2)
library(readr)
library(xlsx)
library(data.table)

setwd(dir="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Data_clean")

# ---- Waldi ----
# load coordinates metric file
wal_m <- fread("Waldi_circlecplots_Coordinates_metric.csv")
names(wal_m)[2] <- "CirclePlotID"
wal_m=wal_m[,-c("V1")]

#check for duplicates
dup_ID=wal_m$CirclePlotID[duplicated(wal_m$CirclePlotID)]
wal_m[wal_m$CirclePlotID %in% dup_ID,]

# load circle plots info
wal_cp <- fread("Waldi_circleplot_info.csv")
# two columns are identical
which(wal_cp$SampleID!=wal_cp$CirclePlotTreeID)
#remove one of them
wal_cp=wal_cp[,-c("SampleID")]

# merge Circular plot info with coordinates of both circle plots and individual adults
wal_cp_full=merge(wal_cp,wal_m,by="CirclePlotID",all=TRUE)

#coordinates of adults only
ad=fread("Waldi_adults_Coordinates_metric.csv")
ad=ad[,-1]
names(ad)=c("CirclePlotTreeID","x_mothertree","y_mothertree")

# add the coordinates of the mother tree (CirclePlotTreeID)
wal_cp_full2=merge(wal_cp_full,ad,by="CirclePlotTreeID",all.x=TRUE,all.y=FALSE)

#fwrite(subset(wal_cp_full2, !is.na("CirclePlotTreeID")),"Waldi_circleplot_info_withcoord.csv" )

# import individuals file
indiv = fread("Waldi_indiv_info.csv")
# create column adult and circle plot id to merge the coordinates
indiv$CirclePlotID = case_when(
  substr(indiv$SampleID, 1, 1) == "O" ~ substr(indiv$SampleID, 1, 6),
  substr(indiv$SampleID, 1, 1) %in% c("B", "R") ~ substr(indiv$SampleID, 1, 2),
  TRUE ~ substr(indiv$SampleID, 1, 5)
)

# add coordinates of circle plots
indiv_merge = merge(indiv, wal_m, by = "CirclePlotID", all.x = T)

# add coordinates of adults
names(ad)=c("CirclePlotID","x_mothertree","y_mothertree")
indiv_merge = merge(indiv_merge, ad, by = "CirclePlotID", all.x = T)

# combine coordinates columns into one
indiv_merge$x <- coalesce(indiv_merge$x, indiv_merge$x_mothertree)
indiv_merge$y <- coalesce(indiv_merge$y, indiv_merge$y_mothertree)
indiv_merge <- indiv_merge[, -c( "x_mothertree", "y_mothertree")]

# add Ageclass classification
levels(unique(factor(indiv_merge$Ageclass)))
indiv_merge <- indiv_merge %>%
  mutate(Ageclass_name = case_when(
    Ageclass == 1 | Ageclass == 2  ~ "Seedling",
    Ageclass == 3  ~ "Sapling",
    Ageclass == 4 | Ageclass == 5 ~ "Juvenile",
    Ageclass == 10 ~ "Adult",
    TRUE           ~ NA_character_  
  ))

# correct the generation (split Juveniles from Offspring)
levels(unique(factor(indiv_merge$generation)))
indiv_merge[which(indiv_merge$Ageclass_name == "Juvenile"), "generation"] = "Juvenile"
levels(unique(factor(indiv_merge$generation)))
table(indiv_merge$generation)

# write file
#fwrite(indiv_merge, "Waldi_indiv_info_withcoord.csv")

# Genotypes (from Waldi_alleles_binned)
genot = fread("Waldi_genotypes.csv")
genot$CirclePlotID = case_when(
  substr(genot$SampleID, 1, 1) == "O" ~ substr(genot$SampleID, 1, 6),
  substr(genot$SampleID, 1, 1) %in% c("B", "R") ~ substr(genot$SampleID, 1, 2),
  TRUE ~ substr(genot$SampleID, 1, 5)
)

### create a file with adults and circle plots coords
#load coordinates of adults only
ad=fread("Waldi_adults_Coordinates_metric.csv")
ad=ad[,-1]

# add Ageclass classification
genot <- genot %>%left_join(indiv_merge %>% select(SampleID, Ageclass_name), by = "SampleID") %>% relocate(Ageclass_name, .after = generation)

# merge adults and circle plots
names(ad)=names(wal_m)
all_coords=rbind(wal_m,ad)
# add coords
genot_merge = merge(genot, all_coords, by = "CirclePlotID", all.x = T)
# correct the generation (split Juveniles from Offspring)
genot_merge[which(genot_merge$Ageclass_name == "Juvenile"), "generation"] = "Juvenile"

#fwrite(genot_merge, "Waldi_genotypes_withcoord.csv")

# ---- Allenwiller -----
# load file with metric coordinates
al_m <- fread("Allenwiller_circlecplots_Coordinates_metric.csv")
al_m=al_m[,-c("V1")]

#check for duplicates
dup_ID=al_m$CirclePlotID[duplicated(al_m$CirclePlotID)]
al_m[al_m$CirclePlotID %in% dup_ID,]

# load circle plots info
al_cp <- fread("Allenwiller_circleplot_info.csv")
# two columns are identical
which(al_cp$SampleID!=al_cp$PlotTree)
#remove one of them and other unnecessary columns
al_cp=al_cp[,-c("SampleID", "SampleIDOriginal", "PlotIDOriginal")]
# rename column as waldi for consistency
names(al_cp)[c(1,2,13)] = colnames(wal_cp)[c(1,2,13)] 

# merge Circular plot info with coordinates of both circle plots and individual adults
al_cp_full=merge(al_cp,al_m,by="CirclePlotID",all=TRUE)

#coordinates of adults only
ad=fread("Allenwiller_adults_Coordinates_metric.csv")
ad=ad[,-1]
names(ad)=c("CirclePlotTreeID","x_mothertree","y_mothertree")

# add the coordinates of the mother tree (CirclePlotTreeID)
al_cp_full2=merge(al_cp_full,ad,by="CirclePlotTreeID",all.x=TRUE,all.y=FALSE)
#fwrite(subset(al_cp_full2, !is.na("CirclePlotTreeID")),"Allenwiller_circleplot_info_withcoord.csv" )

# import individuals file
indiv = fread("Allenwiller_indiv_info.csv")
# create column adult and circle plot id to merge the coordinates
indiv$CirclePlotID = substr(indiv$SampleID ,1,7)
# remove the last digit from VA circle plots (only offspring)
indiv$CirclePlotID[indiv$generation == "Offspring"] <- gsub("([A-Z0-9]+[A-Z])[0-9]+$", "\\1",  indiv$CirclePlotID[indiv$generation == "Offspring"])

# add coordinates of circle plots
indiv_merge = merge(indiv, al_m, by = "CirclePlotID", all.x = T)
# add coordinates of adults
names(ad)=c("CirclePlotID","x_mothertree","y_mothertree")
indiv_merge = merge(indiv_merge, ad, by = "CirclePlotID", all.x = T)
# combine coordinates columns into one
indiv_merge$x <- coalesce(indiv_merge$x, indiv_merge$x_mothertree)
indiv_merge$y <- coalesce(indiv_merge$y, indiv_merge$y_mothertree)
indiv_merge <- indiv_merge[, -c( "x_mothertree", "y_mothertree")]

# add Ageclass classification
levels(unique(factor(indiv_merge$Ageclass)))
indiv_merge <- indiv_merge %>%
  mutate(Ageclass_name = case_when(
    Ageclass == 1  ~ "Seedling",
    Ageclass == 3  ~ "Sapling",
    Ageclass == 4  ~ "Juvenile",
    Ageclass == 10 ~ "Adult",
    TRUE           ~ NA_character_  
  ))

# correct generation based on ageclasses
indiv_merge[which(indiv_merge$Ageclass_name == "Juvenile"), "generation"] = "Juvenile"
table(indiv_merge$generation)

# write file
#fwrite(indiv_merge, "Allenwiller_indiv_info_withcoord.csv")

# reduced dataset without the surrounding individuals
subset <- subset(indiv_merge, x > 4126670 & y > 2840030)
plot(subset$y ~subset$x)
#fwrite(subset, "Allenwiller_indiv_info_withcoord_subset.csv")

# Genotypes
genot = fread("Allenwiller_genotypes.csv")
# create column adult and circle plot id to merge the coordinates
genot$CirclePlotID = substr(genot$SampleID ,1,6)
### create a file with adults and circle plots coords
#load coordinates of adults only
ad=fread("Allenwiller_adults_Coordinates_metric.csv")
ad=ad[,-1]

# add Ageclass classification
genot <- genot %>%left_join(indiv_merge %>% select(SampleID, Ageclass_name), by = "SampleID") %>% relocate(Ageclass_name, .after = generation)
genot[which(genot$Ageclass_name == "Juvenile"), "generation"] = "Juvenile"

# merge adults and circle plots
names(ad)=names(al_m)
all_coords=rbind(al_m,ad)
# add coords to genotypes
genot_merge = merge(genot, all_coords, by = "CirclePlotID", all.x = T)

fwrite(genot_merge, "Allenwiller_genotypes_withcoord.csv")
