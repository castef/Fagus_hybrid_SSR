#' @name nh_subsetR
#' @title Subset NewHybrids format datasets
#' @description Function allows the user to subset NewHybrids format datasets to explore the efficacy of variously sized panels. Note, the function will create a file with your file name and "_subsetted" appended to it (viz. "YourFileName_subsetted.txt"), so this file must be renamed so it is not overwritten if multiple subsetted datasets are desired.
#' @param NHData The file path to the NewHybrids formatted dataset to be subsetted
#' @param loci A vector of Loci names to be retained in the dataset
#' @importFrom dplyr filter summarise ungroup group_by
#' @importFrom stringr str_extract str_split
#' @export



nh_subsetR <- function(NHData, loci = NULL){
  
  
  NHD <- read.table(NHData, header = FALSE, quote = "", sep = "\t", stringsAsFactors = FALSE)
  
  NHD_nometa <- NHD[-(1:5),] ## remove unneeded
  
  temp <- as.data.frame(do.call(rbind, stringr::str_split(NHD_nometa, " ")))
  
  colnamevec <- unlist(stringr::str_split(NHD[5,], " ")) ## there was a typo here
  
  colnames(temp) <- colnamevec
  
  if(length(loci) > 0){ ### May add ability to subset populations later - have if statement here to allow
    
    loci.to.keep <- which(colnamevec %in% loci)
    
    temp_subsetted <- temp[, c(1, loci.to.keep)]
  }
  
  insertNumIndivs <- paste("NumIndivs", nrow(temp_subsetted))
  
  insertNumLoci <- paste("NumLoci", length(temp_subsetted[-1])) ## will probably have to be -1
  
  ### hard coded stuff
  insertDigits <- "Digits 3"
  insertFormat <- "Format Lumped"
  
  insertLociNames <- paste(colnames(temp_subsetted), collapse = " ")
  
  data.out <- do.call(paste, c(data.frame(temp_subsetted[,]), sep = " "))
  
  data.out <- c(insertNumIndivs, insertNumLoci, insertDigits, insertFormat, insertLociNames, data.out)
  
  
  NHsplit <- c(stringr::str_split(string = NHData, pattern = "/"))
  
  outNameHold <- stringr::str_extract(NHsplit, paste0("[:word:]{3,}", ".txt"))
  dat.name.remove <- which(unlist(NHsplit) == outNameHold)
  outNameHold <- gsub(x = outNameHold, pattern = ".txt", replacement = "_subsetted.txt")
  outNameOut <- paste0(paste(unlist(NHsplit)[-dat.name.remove], collapse = "/"), "/", outNameHold)
  
  write(x = data.out, file = outNameOut)
  
}
