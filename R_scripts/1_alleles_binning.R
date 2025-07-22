# Allele binning #

# binning size results for locus DZ447 and merging with Tandem results
library(xlsx)

rm(list=ls())

setwd(dir="~/Data_clean")

# Waldi
waldi_size <- read.xlsx("20230201_AllSize_Waldi.xlsx", sheetIndex = 1)
DZ447 <- c(waldi_size$DZ447, waldi_size$NA..8)

hist(DZ447[! DZ447%in% -1])
plot(sort(DZ447[! DZ447%in% -1]),sequence(table((DZ447[! DZ447%in% -1]))), pch = 16, cex = 0.4)
axis(1, at = seq(183, 253, by = 2), las = 1)

# rounding number down to nearest integer
DZ447b <- floor(DZ447)
plot(sort(DZ447b[! DZ447b%in% -1]),sequence(table((DZ447b[! DZ447b%in% -1]))), pch = 16, cex = 0.4)

## zoom in to check the three alleles 1 bp apart
dz447 <- DZ447[DZ447 > 187 & DZ447 < 193]
plot(sort(dz447[! dz447%in% -1]),sequence(table((dz447[! dz447%in% -1]))), pch = 16, cex = 0.4)    
dz447b <- DZ447b[DZ447b > 187 & DZ447b < 193]
plot(sort(dz447b[! dz447b%in% -1]),sequence(table((dz447b[! dz447b%in% -1]))), pch = 16, cex = 0.4)

# keeping only the needed columns (DZ447 locus) and rounding
waldi_D <- waldi_size[, c("CirclePlot", "SampleID","DZ447", "NA..8")]
waldi_D$DZ447 = floor(waldi_size$DZ447)
waldi_D$NA..8 = floor(waldi_size$NA..8)

# combining with TANDEM results
waldi_tan <- read.table("20230201_AllSize_Waldi_tandem.txt", header = F, fill = T, sep = "\t")
waldi_tan = waldi_tan[-c(1,2,3,5),]
colnames(waldi_tan) = waldi_tan[1,]
colnames(waldi_tan)[19] = "NA..8"
waldi_tan = waldi_tan[-1,]
waldi_tan[, -1] = apply(waldi_tan[, -1], 2, function(x) as.numeric(x))

# replacing DZ447 values from TANDEM with value from mathematical rounding
waldi_binned <- waldi_tan
waldi_binned$DZ447 <- waldi_D$DZ447
waldi_binned$N..8 <- waldi_D$NA..8

#plotting the binned results
plot(sort(waldi_binned$DZ447[! waldi_binned$DZ447%in% -1]),sequence(table((waldi_binned$DZ447[! waldi_binned$DZ447%in% -1]))), pch = 16, cex = 0.4)
axis(1, at = seq(183, 253, by = 2), las = 1)

# renaming
colnames(waldi_binned)[1] = "SampleID"

# saving results
write.csv(waldi_binned, "Waldi_alleles_binned.csv", row.names = F)

# Allenwiller
Allenwiller_size <- read.xlsx("20230201_AllSize_Allenwiller.xlsx", sheetIndex = 1)

DZ447 <- c(Allenwiller_size$DZ447, Allenwiller_size$NA..8)
hist(DZ447[! DZ447%in% -1])
plot(sort(DZ447[! DZ447%in% -1]),sequence(table((DZ447[! DZ447%in% -1]))), pch = 16, cex = 0.4)
axis(1, at = seq(183, 253, by = 2), las = 1)

# rounding number down to nearest integer
DZ447b <- floor(DZ447)
plot(sort(DZ447b[! DZ447b%in% -1]),sequence(table((DZ447b[! DZ447b%in% -1]))), pch = 16, cex = 0.4)

## zoom in to check the three alleles 1 bp apart
dz447 <- DZ447[DZ447 > 187 & DZ447 < 193]
plot(sort(dz447[! dz447%in% -1]),sequence(table((dz447[! dz447%in% -1]))), pch = 16, cex = 0.4)    
dz447b <- DZ447b[DZ447b > 187 & DZ447b < 193]
plot(sort(dz447b[! dz447b%in% -1]),sequence(table((dz447b[! dz447b%in% -1]))), pch = 16, cex = 0.4)

# keeping only the needed columns (DZ447 locus) and rounding
Allenwiller_D <- Allenwiller_size[, c("CirclePlot", "SampleID","DZ447", "NA..8")]
Allenwiller_D$DZ447 = floor(Allenwiller_size$DZ447)
Allenwiller_D$NA..8 = floor(Allenwiller_size$NA..8)

# combining with TANDEM results
Allenwiller_tan <- read.table("20230201_AllSize_Allenwiller_tandem.txt", header = F, fill = T, sep = "\t")
Allenwiller_tan = Allenwiller_tan[-c(1,2,3,5),]
colnames(Allenwiller_tan) <- Allenwiller_tan[1, ]
Allenwiller_tan <- Allenwiller_tan[-1, ]
colnames(Allenwiller_tan)[1] = "SampleID"

# removing VA044N01 and VA044N02
Allenwiller_tan = Allenwiller_tan[!apply(Allenwiller_tan=="VA044N01" | Allenwiller_tan=="VA044N02",1,any),]

# renaming column to distinguish from other NA
colnames(Allenwiller_tan)[19] = "NA..8"

# convert to numeric
Allenwiller_tan[, -1] = apply(Allenwiller_tan[, -1], 2, function(x) as.numeric(x))

# replacing DZ447 values from TANDEM with value from mathematical rounding
Allenwiller_binned <- Allenwiller_tan
Allenwiller_binned$DZ447 <- Allenwiller_D$DZ447
Allenwiller_binned$NA..8<- Allenwiller_D$NA..8

#plotting the binned results
plot(sort(Allenwiller_binned$DZ447[! Allenwiller_binned$DZ447%in% -1]),sequence(table((Allenwiller_binned$DZ447[! Allenwiller_binned$DZ447%in% -1]))), pch = 16, cex = 0.4)
axis(1, at = seq(183, 253, by = 2), las = 1)

# reshaping and renming columns as in waldi for consistency

colnames(Allenwiller_binned) <- colnames(waldi_binned)

# saving results
write.csv(Allenwiller_binned, "Allenwiller_alleles_binned.csv", quote=F, row.names = F)
