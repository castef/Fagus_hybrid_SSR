#' @name hybridPowerComp
#' @title Assignment power comparison among different SNP subsets using NewHybrids simulated datasets3.
#' @description Evaluates the accuracy with which NewHybrids assigns individuals of known hybrid class to the correct hybrid class in simulated datasets at varying levels of stringency (PofZ). The code will write graphical and numerical results to the directory provided by the user.
#' @param dir File path to the directory in which the NewHybrids results (in individual folders as returned by parallelNH_XX) are located.
#' @param save_output A logical indicating whether the plots and plotting data should be saved to the hard drive. The default is TRUE
#' @param return_workspace A logical indicating whether the plots and plotting data should be returned to the workspace as a list object. The default is FALSE
#' @param Thresholds A vector of critical posterior probability thresholds by which to facet plots. The default is c(0.5,0.6,0.7,0.8,0.9)
#' @param samplesize The number of individuals per NewHybrids class. By (default: NULL) this data will be extracted from the "*individuals.txt" output from parallelnewhybrids if present in the same folder as the PofZ file. This can also explicitly defined as a vector (6 values corresponding to # in P1,P2,F1,F2,BC1,BC2) or a path to a *_Individuals.txt.
#' @param CT The threshold posterior probability of assignment (PofZ) to F2 above which Pure Population 1 or Pure Population 2 individuals are flagged to indicate possible non-convergence. The default is 0.1.
#' @param CTI The proportion of individuals in either Pure Population 1 OR Pure Population 2 allowed to exceed the F2 assignment threshold (PofZCutOff). The default is 0.5.
#' @rdname hybridpowercomp
#' @import ggplot2
#' @import magrittr
#' @importFrom dplyr filter summarise ungroup group_by do
#' @importFrom grid arrow unit
#' @importFrom stringr str_extract
#' @importFrom reshape2 melt
#' @importFrom  scales alpha
#' @export
#'


## function was changed to calculate just accuracy, efficiency and power, 
# as Type I error didn't work out


# dir <- paste0(path.hold, "/Allenwiller_panel_locimost/NH.Results/")
# save_output = TRUE
# return_workspace = FALSE
# Thresholds = c(0.5,0.6,0.7,0.8,0.9)
# samplesize = NULL
# CT = 0.1
# CTI = 0.5


my_hybridPowerComp <-function(dir, save_output = TRUE, return_workspace = FALSE, Thresholds = c(0.5,0.6,0.7,0.8,0.9),  samplesize = NULL, CT = 0.1, CTI = 0.5){
  
  #######################################################################################
  ##Make sure people aren't specifying that the function shouldn't output anything######
  #######################################################################################
  if(save_output == FALSE & return_workspace == FALSE){
    stop("You have asked me to not return any results. If you're not going to look at the analyses, I'm not going to bother running them.")
  }
  #######################################################
  
  if(save_output == TRUE){
    ##########################
    ##Create Results Folders##
    ##########################
    
    #set directory for which holds the New Hybrids output folders
    filedir <- dir
    lfiles <- setdiff(list.files(dir), c("Figures and Data", "NewHybrids Plots")) #ignores Figures folder in case this is run more than once and in case plots made
    if(length(which(list.files(dir) == "Figures and Data")) == 0) {dir.create(paste0(dir, "Figures and Data"))} # if there isn't a 'Figures and Data' folder for output create one
    if(length(which(list.files(paste0(dir, "Figures and Data")) == "pdf")) == 0) {dir.create(paste0(dir, "Figures and Data/pdf"))} #create a folder for pdfs
    if(length(which(list.files(paste0(dir, "Figures and Data")) == "jpg")) == 0) {dir.create(paste0(dir, "Figures and Data/jpg"))} #create a folder for jpgs
    if(length(which(list.files(paste0(dir, "Figures and Data")) == "data")) == 0) {dir.create(paste0(dir, "Figures and Data/data"))} #create a folder for data
  }
  #######################################################
  
  #####################
  ##Convergence Check##
  #####################
  
  #Convergence checker - set to initial state
  arethereproblems = "no"
  
  # Collate the output from New Hybrids together ('p of z' files)
  output <- NULL
  for (i in lfiles)
  {
    
    ### ANALYSIS REPs - if used hybriddetective to simulate the data, unique replicates will have different numbers of loci, a unique simulation number, and a replicate nubmer
    tempfiles <- list.files(paste0(filedir, i)) ## Get names of all files in the directory
    pzfile <- tempfiles[grep("PofZ", tempfiles)] ## Which of the files is the PofZ file?
    tempfile <- read.table(paste0(filedir,i, "/", pzfile), head = TRUE) ### Read the PofZ file in
    
    ## Each analysis must have an accompanying LociAndAlleles file <- use this to figure out how many loci there are
    LociandAlleles <- tempfiles[grep("LociAndAlleles", tempfiles)]
    LandAfile <- readChar(paste0(filedir,i,  "/", LociandAlleles), file.info(paste0(filedir,i,  "/", LociandAlleles))$size)
    numLociExt <- stringr::str_extract(string = LandAfile, pattern = paste0("from ", "[:digit:]{1,5}", " loci")) ### find the string with the number of loci, extract
    numLociWorking <- gsub(x = numLociExt, pattern = "from ", replacement = "") ## remove "from"
    numLociWorking <- as.numeric(gsub(x = numLociWorking, pattern = " loci", replacement = "")) ### This is how many loci there are
    
    #identify the simulation and repeat info
    S_ident <- gsub("_", "", stringr::str_extract(pzfile, paste0("_S", "[:digit:]{1}", "R", "[:digit:]{1}", "_"))) ### if used hybriddetective, will have S#_R# - exctract
    tempfile$sim <- substring(S_ident, 1, 2)
    tempfile$rep <- substring(S_ident, 3, 4)
    tempfile$nLoci <- numLociWorking
    
    tempfile <- tempfile[ , -grep("IndivName", colnames(tempfile))] #delete IndivName
    
    #rename the columns
    colnames(tempfile) <- c("Indv", "Pure1", "Pure2", "F1", "F2", "BC1", "BC2", "sim", "rep", "nLoci")
    tempfile <- tempfile[ , c("Indv", "sim", "rep", "nLoci", "Pure1", "Pure2", "F1", "F2", "BC1", "BC2")]# reorder
    
    #Get the samplesize for a given class
    samplesize = NULL
    IndividualsPath <- tempfiles[grep("individuals.txt", tempfiles)]
    if(length(samplesize) == 1 & is.numeric(samplesize)){samplesize <- rep(samplesize, 6)}
    if(length(samplesize) == 1 & !is.numeric(samplesize)){samplesize <- as.vector(n_class(samplesize)[ ,2])}
    
    ###THE ERROR IS HERE
    if(is.null(samplesize)){samplesize <- as.vector(n_class(paste0(filedir, i, "/", IndividualsPath)))$n}   ### changed [ ,2] into $n

    #common order
    if(sum(tempfile[1:samplesize[1], "Pure1"], na.rm = TRUE) < sum(tempfile[1:samplesize[1], "Pure2"], na.rm = TRUE)){
      pure1 <- tempfile$Pure2;pure2 <- tempfile$Pure1
      bc1 <- tempfile$BC2;bc2 <- tempfile$BC1
      
      tempfile$Pure1 <- pure1;tempfile$Pure2 <- pure2
      tempfile$BC1 <- bc1;tempfile$BC2 <- bc2
    }
    
    #Filter for convervence issues. Based on the 'convergence filter (CT) and % of indviduals permited to fail (CTI)
    #Here we look at the "pure 1 and 2 populations for
    if(length(which(tempfile[1:samplesize[1], "F2"] > CT))/length(1:samplesize[1]) > CTI &
       length(which(tempfile[(samplesize[1] + 1):samplesize[2],"F2"] > CT))/length((samplesize[1] + 1):samplesize[2]) > CTI){
      
      tempfile[ ,5:length(tempfile)] = NA #replace data with NAs
      print(paste("Possible non-convergence detected in", pzfile))
      arethereproblems = "Yes"}
    
    output <- rbind(output, tempfile)
    
  }#end of for loop
  
  #If convergence issues were flagged then the process stops here.
  if(arethereproblems == "Yes")
      {
    stop("Please remove, or re-run those results for which non-convergence was detected", call. = F)
  }
  
  
  temp <- output
  tempinds <- hybriddetective::n_class(paste0(filedir, i, "/", IndividualsPath))           # paste0(filedir, lfiles[1], "/", IndividualsPath)) 
  out.inds.class = NULL
  
  ## get the correct numbers of individuals for known categories
  for(j in 1:length(samplesize)){
    make.class <- rep(as.character(tempinds[j,1]), times = tempinds[j,2])
    out.inds.class <- c(out.inds.class, make.class)
  } # END J Loop
  
  temp$known <- out.inds.class
  colnames(temp)[c(5, 6)] <- c("P1", "P2") ## rename these two columns to match column names between all data frames made
  temp$max.class <- NA ## add a column that will be the class to which each indiviudal is assigned
  
  ## calculate what the most probable (highest PofZ) genotype frequency category is for each indivivudal
  for(k in 1:nrow(temp)){
    temp$max.class[k] = names(which.max(temp[k, 5:10]))
  } ## END K Loop
  
  
  
  #######################################################
  
  ########################
  ## Calculate ACCURACY ##
  ########################
  
  #### CHECK HERE - ACCURACY LOOP - Can speed up???
  
  writeLines("Calculating Accuracy
          ")
  
  ## ACCURACY =  number assigned correctly / total number assigned  -> for each category
  
  temp2 <- temp
  tempAccuracy <- temp
  tempAccuracy$isgood = TRUE ### Set all isgood to TRUE - isgood is a check if the assigned PofZ > the critical PofZ - will change to FALSE in function
  tempAccuracy$domatch = FALSE ### Set all domatch to FALSE - domatch is check if the assigned class == known class. Will change to TRUE in function
  
  
  #Create long form data for dplyr loop
  PofZVector_Accuracy <- rep(50:99/100, each = nrow(tempAccuracy)) #vector of PofZs
  
  tempAccuracyLong <- do.call("rbind", replicate(length(50:99), tempAccuracy, simplify = FALSE))
  tempAccuracyLong$group <- as.character(PofZVector_Accuracy)
  tempAccuracyLong$pofz <- PofZVector_Accuracy
  
  AccuracyData <- tempAccuracyLong%>%
    dplyr::group_by(pofz,nLoci,max.class)%>%
    dplyr::do(accuracyfunction(.))%>%  ### the error is here
    dplyr::ungroup()%>%data.frame()
  
  #Set up plotting data
  AccuracyData$nLoci <- factor(x = AccuracyData$nLoci, levels = ordered(unique(as.numeric(AccuracyData$nLoci))))
  AccuracyData$max.class <- factor(x = AccuracyData$max.class, levels = c("P1", "P2", "F1", "F2", "BC1", "BC2"))
  
  
  ### Can now plot Accuracy 6 class boxplot
  
  tempAccuracyLong2 <- tempAccuracyLong
  tempAccuracyLong2$max.class2 <- as.character(tempAccuracyLong2$max.class)
  tempAccuracyLong2$max.class2[tempAccuracyLong2$max.class2 %in% c("F1", "F2", "BC1", "BC2")] = "Hyb"
  
  AccuracyDataBoxPlot <- tempAccuracyLong2%>%
    dplyr::group_by(pofz,nLoci,max.class2)%>%
    dplyr::do(accuracyfunction(.))%>%
    dplyr::ungroup()%>%data.frame()
  
  AccuracyDataBoxPlot$max.class2 <- factor(x = AccuracyDataBoxPlot$max.class2, levels = c("P1", "P2", "Hyb"))
  
  ## Can now plot Accuracy Pure Hyb boxplot
  
  #Summary among simulations --
  SummaryAccuracy <- AccuracyData%>%
    dplyr::group_by(pofz,nLoci,max.class)%>%
    dplyr::summarise(mean=mean(means,na.rm=T),
                     sd=sd(means,na.rm=T))%>%
    dplyr::ungroup()%>%data.frame()
  
  get.y.min.AccuracyLine <- min((SummaryAccuracy$mean - SummaryAccuracy$sd), na.rm = TRUE)
  
  get.y.min.AccuracyThreshold <- min(dplyr::filter(SummaryAccuracy,pofz %in% Thresholds)$mean-dplyr::filter(SummaryAccuracy,pofz %in% Thresholds)$sd)
  
  ## Can now plot Accuracy Line 6 Class Line Plot -> Class as colour, faceted by Loci, with SD
  ## Can now plot Accuracy Line Plot 6 Class -> Loci as colour, faceted by 6 Classes
  ## Can now plot Accuracy Dot Plot 6 Class -> Colour as Class, Faceted by PofZ %in% c(0.05, 0.6, 0.7, 0.8, 0.9)
  
  
  #ComboHybrids ------------
  ComboHybridAccuracy <- AccuracyData
  ComboHybridAccuracy$class <- as.character(AccuracyData$max.class)
  ComboHybridAccuracy$class[ComboHybridAccuracy$class %in% c("F1", "F2", "BC1", "BC2")] = "Hybrid"
  ComboHybridAccuracy$class[ComboHybridAccuracy$class == "P1"] = "Pure1"
  ComboHybridAccuracy$class[ComboHybridAccuracy$class == "P2"] = "Pure2"
  
  ComboHybridAccuracy <- ComboHybridAccuracy%>%
    dplyr::group_by(nLoci,pofz,class)%>%
    dplyr::summarise(mprob = mean(means,na.rm=T),
                     sdprob = sd(means,na.rm=T))%>%
    dplyr::ungroup()%>%data.frame()
  
  ComboHybridAccuracy$class <- factor(ComboHybridAccuracy$class, levels=c("Pure1", "Pure2", "Hybrid")) # set plotting levels
  
  
  ## Can now plot Accuracy Dot Plot 6 Class -> Colour as Class, Faceted by PofZ %in% c(0.05, 0.6, 0.7, 0.8, 0.9) and Pure Hyb
  ## Can now plot Accuracy Line Plot Hyb Pure -> Loci as colour, faceted by Hyb Pure
  
  AccuracyData2 <- AccuracyData
  AccuracyData2$class <- as.character(AccuracyData2$max.class)
  AccuracyData2$class[AccuracyData2$class %in% c("F1", "F2", "BC1", "BC2")] = "Hybrid"
  AccuracyData2$class[AccuracyData2$class == "P1"] = "Pure1"
  AccuracyData2$class[AccuracyData2$class == "P2"] = "Pure2"
  
  SummaryAccuracy2 <- AccuracyData2%>%
    dplyr::group_by(pofz,nLoci,class)%>%
    dplyr::summarise(mean=mean(means,na.rm=T),
                     sd=sd(means,na.rm=T))%>%
    dplyr::ungroup()%>%data.frame()
  
  ## Can now plot Accuracy Line Pure Hyb Line Plot -> Class as colour, faceted by Loci, with SD
  
  
  SummaryAccuracy$Group = NA
  
  SummaryAccuracy$Group[SummaryAccuracy$max.class %in% c("P1", "P2")] = "Pure"
  SummaryAccuracy$Group[SummaryAccuracy$max.class %in% c("F1", "F2")] = "Generational Hybrids"
  SummaryAccuracy$Group[SummaryAccuracy$max.class %in% c("BC1", "BC2")] = "Back-cross"
  
  SummaryAccuracy$Group <- factor(x = SummaryAccuracy$Group, levels = c("Pure", "Generational Hybrids", "Back-cross"))
  
  ## Can now plot Accuracy Line Plot Group Facet -> Class as colour, faceted by LOCI and "hybrid grouping", with SD
  
  
  ##################################################
  
  ###########################
  ### CALCULATE EFFICIENCY ##
  ###########################
  
  writeLines("
            Calculating Efficiency
            ")
  
  ## average and SD the  replicate runs of each simulation in New Hybrids. Filter is just a holder for the dplyr:: call
  sim_data <-dplyr::filter(output)%>%dplyr::group_by(nLoci,sim,Indv)%>%dplyr::summarise(Pure1_sd=sd(Pure1),Pure1=mean(Pure1),
                                                                                        Pure2_sd=sd(Pure2),Pure2=mean(Pure2),
                                                                                        F1_sd=sd(F1),F1=mean(F1),
                                                                                        F2_sd=sd(F2),F2=mean(F2),
                                                                                        BC1_sd=sd(BC1),BC1=mean(BC1),
                                                                                        BC2_sd=sd(BC2),BC2=mean(BC2))%>%dplyr::ungroup()%>%data.frame()
  #pull out just the means of the replicates
  sim_means <- sim_data[ ,-grep("_sd", colnames(sim_data))]
  ## Look at assignment success as a function of threshold probability
  
  num.sim <- length(which(sim_means$sim == "S1"))/6/length(unique(sim_means$nLoci))
  
  ## find he dim of each class in a given sim (length = sum n_class$n
  classvec2 <- rep(c("Pure1", "Pure2", "F1", "F2", "BC1", "BC2"), times = samplesize)
  
  ProbOutput <- NULL
  for (s in unique(sim_means$nLoci)){
    
    lsub <- filter(sim_means,nLoci == s)
    
    for(i in unique(sim_means$sim)){
      tempsub <- filter(lsub,sim == i)
      
      for(q in 50:99/100){ # probability of 50 - 99%
        
        p1.p <- length(which(tempsub[which(classvec2 == "Pure1"), "Pure1"] > q))/samplesize[1]
        p2.p <- length(which(tempsub[which(classvec2 == "Pure2"), "Pure2"] > q))/samplesize[2]
        F1.p <- length(which(tempsub[which(classvec2 == "F1"), "F1"] > q))/samplesize[3]
        F2.p <- length(which(tempsub[which(classvec2 == "F2"), "F2"] > q))/samplesize[4]
        BC1.p <- length(which(tempsub[which(classvec2 == "BC1"), "BC1"] > q))/samplesize[5]
        BC2.p <- length(which(tempsub[which(classvec2 == "BC2"), "BC2"] > q))/samplesize[6]
        tempout <- data.frame(nLoci = s, sim = i, level = q, prob = c(p1.p, p2.p, F1.p, F2.p, BC1.p, BC2.p),
                              class=c("Pure1", "Pure2", "F1", "F2", "BC1", "BC2"))
        ProbOutput <- rbind(ProbOutput,tempout)
        
      } # end q loop
    } # end i loop
  } # end s loop
  
  ### set proper levels for ProbOutput
  ProbOutput$class <- factor(x = ProbOutput$class, levels = c("Pure1", "Pure2", "F1", "F2", "BC1", "BC2"))
  
  #combined hybrid probabilities
  ProbOutput2 <- NULL
  for (s in unique(sim_means$nLoci)){
    
    lsub <- filter(sim_means,nLoci == s)
    
    for(i in unique(sim_means$sim)){
      
      tempsub <- filter(lsub,sim == i)
      tempsub$phyb <- rowSums(tempsub[ ,c("F1", "F2", "BC1", "BC2")])
      
      for(q in 50:99/100){ # probability of 50 - 99%
        
        p1.p <- length(which(tempsub[which(classvec2 == "Pure1"), "Pure1"] > q))/samplesize[1]
        p2.p <- length(which(tempsub[which(classvec2 == "Pure2"), "Pure2"] > q))/samplesize[2]
        Hybrid <- length(which(tempsub[which(classvec2 %in% c("F1", "F2", "BC1", "BC2")), "phyb"] > q))/sum(samplesize[3:6])
        tempout <- data.frame(nLoci = s, sim = i, level = q, prob = c(p1.p, p2.p, Hybrid),
                              class = c("Pure1", "Pure2", "Hybrid"))
        ProbOutput2 <- rbind(ProbOutput2, tempout)
        
      } # end q loop
    } # end i loop
  } # end s loop
  
  ### set proper levels for ProbOutput2
  ProbOutput2$class <- factor(x = ProbOutput2$class, levels = c("Pure1", "Pure2", "Hybrid"))
  
  # get the mean and standard error for the estimates of assignment succes based on NH probabilty among simulations
  FinalData <- data.frame(ProbOutput%>%dplyr::group_by(nLoci,level,class)%>%dplyr::summarise(mprob = mean(prob,na.rm=T),
                                                                                             sdprob = sd(prob,na.rm=T))%>%dplyr::ungroup())
  FinalData$class <- factor(FinalData$class, levels=c("Pure1", "Pure2", "F1", "F2", "BC1", "BC2")) # NH class
  
  # set plotting levels
  FinalData$group <- "Pure"
  FinalData[which(FinalData$class %in% c("BC1", "BC2")), "group"] <- "Back-cross"
  FinalData[which(FinalData$class %in% c("F1", "F2")), "group"] <- "Generational hybrids"
  
  FinalData$group <-  factor(FinalData$group, levels = c("Pure", "Generational hybrids", "Back-cross"))
  FinalData$class <- factor(FinalData$class, levels = c("Pure1", "Pure2", "F1", "F2", "BC1", "BC2"))
  
  #ComboHybrids ------------
  FinalData2 <- data.frame(ProbOutput2%>%group_by(nLoci,level,class)%>%summarise(mprob = mean(prob,na.rm=T),
                                                                                 sdprob = sd(prob,na.rm=T))%>%ungroup())
  
  FinalData2$class <- factor(FinalData2$class, levels=c("Pure1","Pure2","Hybrid")) # set plotting levels
  
  
  
  ##################################################
  
  ####################################
  ### CALCULATE OVERALL POWER!!!!!! ##
  ####################################
  
  writeLines("
                      Calculating Power!!!!
                      ")
  
  ### Create power box plots
  performance_efficiency <- FinalData
  performance_accuracy <- SummaryAccuracy
  
  
  performance_efficiency$class <- as.character(performance_efficiency$class)
  performance_efficiency$class[which(performance_efficiency$class == "Pure1")] = "P1"
  performance_efficiency$class[which(performance_efficiency$class == "Pure2")] = "P2"
  
  performance_accuracy$to.merge <- interaction(performance_accuracy$max.class, performance_accuracy$nLoci, performance_accuracy$pofz)
  performance_efficiency$to.merge <- interaction(performance_efficiency$class, performance_efficiency$nLoci, performance_efficiency$level)
  
  
  performance_merge <- merge(x = performance_accuracy, y = performance_efficiency, by = "to.merge")
  
  performance_merge$performance <- performance_merge$mean*performance_merge$mprob
  
  performance_merge$class <- factor(performance_merge$class, levels = c("P1", "P2", "F1", "F2", "BC1", "BC2"))
  
  
  performance_merge_Hyb <- performance_merge
  performance_merge_Hyb$class <- as.character(performance_merge_Hyb$max.class)
  performance_merge_Hyb$class[performance_merge_Hyb$class %in% c("F1", "F2", "BC1", "BC2")] = "Hybrid"
  performance_merge_Hyb$class[performance_merge_Hyb$class == "P1"] = "Pure1"
  performance_merge_Hyb$class[performance_merge_Hyb$class == "P2"] = "Pure2"
  
  
  
  ## Power with SD plots
  getEfficiency <- ProbOutput
  getAccuracy <- AccuracyData
  getAccuracy$nLoci = as.double(as.character(getAccuracy$nLoci))
  getAccuracy$max.class <- as.character(getAccuracy$max.class)
  getEfficiency$sim <- as.character(getEfficiency$sim)
  getEfficiency$class <- as.character(getEfficiency$class)
  getEfficiency$class[getEfficiency$class == "Pure1"] = "P1"
  getEfficiency$class[getEfficiency$class == "Pure2"] = "P2"
  
  start_Power <- merge(x = getEfficiency, y = getAccuracy, by.y = c("simulation", "max.class", "pofz", "nLoci"), by.x = c("sim", "class", "level", "nLoci"))
  
  start_Power$Power <- start_Power$prob * start_Power$means
  
  Final_Power <- data.frame(start_Power%>%dplyr::group_by(nLoci,level,class)%>%dplyr::summarise(meanPower = mean(Power,na.rm=T),
                                                                                                sdPower = sd(Power,na.rm=T))%>%dplyr::ungroup())
  
  Final_Power$class <- factor(x = Final_Power$class, levels = c("P1", "P2", "F1", "F2", "BC1", "BC2"))
  
  
  
  AccuracyData3 <- AccuracyData
  AccuracyData3$class <- as.character(AccuracyData3$max.class)
  AccuracyData3$class[AccuracyData3$class %in% c("F1", "F2", "BC1", "BC2")] = "Hybrid"
  AccuracyData3$class[AccuracyData3$class == "P1"] = "Pure1"
  AccuracyData3$class[AccuracyData3$class == "P2"] = "Pure2"
  
  ProbOutput2$class <- as.character(ProbOutput2$class)
  getEfficiency2 <- ProbOutput2
  
  getAccuracy2 <- AccuracyData3
  start_Power2 <- merge(x = getEfficiency2, y = getAccuracy2, by.y = c("simulation", "class", "pofz", "nLoci"), by.x = c("sim", "class", "level", "nLoci"))
  start_Power2$Power = start_Power2$prob * start_Power2$means
  Final_Power2 <- data.frame(start_Power2%>%dplyr::group_by(nLoci,level,class)%>%dplyr::summarise(meanPower = mean(Power,na.rm=T),
                                                                                                  sdPower = sd(Power,na.rm=T))%>%dplyr::ungroup())
  Final_Power2$class <- factor(x = Final_Power2$class, levels = c("Pure1", "Pure2", "Hybrid"))
  
  Final_Power$group = NA
  Final_Power$group[Final_Power$class %in% c("P1", "P2")] = "Pure"
  Final_Power$group[Final_Power$class %in% c("F1", "F2")] = "Generational Hybrids"
  Final_Power$group[Final_Power$class %in% c("BC1", "BC2")] = "Back-cross"
  Final_Power$group <- factor(x = Final_Power$group, levels = c("Pure", "Generational Hybrids", "Back-cross"))
  
  ##################################################
  
  
  ############################
  ## MEAN AND SD PofZ SCORE ##
  ############################
  
  writeLines("
                  Calculating Mean Posterior Probabilities
                  ")
  
  ## average and SD the  replicate runs of each simulation in New Hybrids. Filter is just a holder for the dplyr:: call
  sim_data <-dplyr::filter(output)%>%dplyr::group_by(nLoci,sim,Indv)%>%dplyr::summarise(Pure1_sd=sd(Pure1),Pure1=mean(Pure1),
                                                                                        Pure2_sd=sd(Pure2),Pure2=mean(Pure2),
                                                                                        F1_sd=sd(F1),F1=mean(F1),
                                                                                        F2_sd=sd(F2),F2=mean(F2),
                                                                                        BC1_sd=sd(BC1),BC1=mean(BC1),
                                                                                        BC2_sd=sd(BC2),BC2=mean(BC2))%>%dplyr::ungroup()%>%data.frame()
  
  #pull out just the means of the replicates
  sim_means <- sim_data[ ,-grep("_sd", colnames(sim_data))]
  ######### There was an error here when Nick ran it
  ## assign the classes to the data
  classvec <- rep(c("Pure1", "Pure2", "F1", "F2", "BC1", "BC2"), times = samplesize)
  classvec <- rep(classvec, times = nrow(sim_means)/length(classvec))
  sim_means$class <- classvec
  
  #Compare the simulations using boxplots
  boxdata <- NULL
  for (i in unique(sim_means$nLoci))
  {
    for (j in unique(sim_means$class))
    {
      bdat_temp <- filter(sim_means, class == j, nLoci == i)
      tout <- data.frame(sim = bdat_temp$sim, Indv = bdat_temp$Indv, class = j, nLoci = i, value = bdat_temp[ ,j])
      boxdata <- rbind(boxdata, tout)
    } ## END J Loop
  } ## END I LOOP
  
  boxdata$nLoci <- factor(boxdata$nLoci)
  
  
  #Combined loci
  sim_means2 <- sim_means
  sim_means2$hybrid <- rowSums(sim_means2[ ,c("F1", "F2", "BC1", "BC2")])
  sim_means2[which(sim_means2$class == "Pure1"),"hybrid"] = sim_means2[which(sim_means2$class =="Pure1"), "Pure1"] #add values of the Pure
  sim_means2[which(sim_means2$class == "Pure2"),"hybrid"] = sim_means2[which(sim_means2$class =="Pure2"), "Pure2"] #add values of the Pure
  
  sim_means2$hclass <- "Hybrid"
  sim_means2[which(sim_means$class == "Pure1"), "hclass"] <- "Pure1"
  sim_means2[which(sim_means$class == "Pure2"), "hclass"] <- "Pure2"
  
  sim_means2$hclass <- factor(sim_means2$hclass, levels = c("Pure1", "Pure2", "Hybrid"))
  


  
  ########################
  
  ########################################################################
  ##Prepare a list of the plots that will be made and their descriptions##
  ########################################################################
  
  Plot.Names = c("Plot_1", "Plot_2", "Plot_3", "Plot_4", "Plot_5", "Plot_6", "Plot_7", "Plot_8", "Plot_9", "Plot_10", "Plot_11", "Plot_12", "Plot_13", "Plot_14", "Plot_15",
                 "Plot_16", "Plot_17", "Plot_18", "Plot_19", "Plot_20", "Plot_21", "Plot_22", "Plot_23", "Plot_24", "Plot_25", "Plot_26", "Plot_27",  "Plot_30", "Plot_31")
  
  Plot.Description = c(
    "Accuracy Boxplot - by hybrid class, and critical posterior probability", #1
    "Efficiency Boxplot - by hybrid class, and critical posterior probability", #2
    "Power Boxplot - by hybrid class, and critical posterior probability", #3
    "Accuracy Boxplot - by pure classes and all hybrids, and critical posterior probability", #4
    "Efficiency Boxplot - by pure classes and all hybrids, and critical posterior probability", #5
    "Power Boxplot - by pure classes and all hybrids, and critical posterior probability", #6
    "Accuracy Lineplot - by hybrid class, faceted by panel size", #7
    "Efficiency Lineplot - by hybrid class, faceted by panel size", #8
    "Power Lineplot - by hybrid class, faceted by panel size", #9
    "Accuracy Lineplot - by pure classes and all hybrids, faceted by panel size", #10
    "Efficiency Lineplot - by pure classes and all hybrids, faceted by panel size", #11
    "Power Lineplot - by pure classes and all hybrids, faceted by panel size", #12
    "Accuracy Lineplot - by hybrid class, faceted by panel size and \"hybrid grouping\"", #13
    "Efficiency Lineplot - by hybrid class, faceted by panel size and \"hybrid grouping\"", #14
    "Power Lineplot - by hybrid class, faceted by panel size and \"hybrid grouping\"", #15
    "Accuracy Lineplot - by panel size, faceted by hybrid class", #16
    "Efficiency Lineplot - by panel size, faceted by hybrid class", #17
    "Power Lineplot - by panel size, faceted by hybrid class", #18
    "Accuracy Lineplot - by panel size, faceted by pure classes and all hybrids", #19
    "Efficiency Lineplot - by panel size, faceted by pure classes and all hybrids", #20
    "Power Lineplot - by panel size, faceted by pure classes and all hybrids", #21
    "Accuracy Dotplot - by hybrid class, faceted by critical posterior probability %in% Thresholds", #22
    "Efficiency Dotplot - by hybrid class, faceted by critical posterior probability %in% Thresholds", #23
    "Power Dotplot - by hybrid class, faceted by critical posterior probability %in% Thresholds", #24
    "Accuracy Dotplot - by pure classes and all hybrids, faceted by critical posterior probability %in% Thresholds", #25
    "Efficiency Dotplot - by pure classes and all hybrids, faceted by critical posterior probability %in% Thresholds", #26
    "Power Dotplot - by pure classes and all hybrids, faceted by critical posterior probability %in% Thresholds", #27
    "Mean Posterior Probability of Assignment - per simulation, faceted by hybrid class", #30
    "Mean Posterior Probability of Assignment - per simulation, faceted by  pure classes and all hybrids" #31
  )
  
  ### make the legend a dataframe - this will then be added to the list that is composed of all plots
  Plot.Legend <- data.frame(Plot = Plot.Names, Description = Plot.Description)
  
  
  ########################
  
  ######################################
  #### Do the darn plotting already ####
  ######################################
  
  writeLines("
                      Makin' you some plots
                      ")
  
  # Accuracy Boxplot - by hybrid class, and critical posterior probability
  Plot_1 <-
    ggplot2::ggplot(filter(AccuracyData,pofz %in% c(0.5,0.75,0.9)), aes(x = max.class, y = means, fill = max.class)) +
    geom_boxplot() +
    facet_grid(pofz~nLoci) +
    labs(x = "Genotype Frequency Class", y = "Accuracy") +
    scale_fill_brewer(palette = "Dark2") +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position = "none", , legend.key = element_blank(),
          strip.background = element_rect(, colour = "black", fill = "white"),
          strip.text.x = element_text(colour = "black"), strip.text.y = element_text(colour = "black"))
  
  
  
  # Efficiency Boxplot - by hybrid class, and critical posterior probability
  Plot_2 <-
    ggplot2::ggplot(filter(ProbOutput,level %in% c(0.5,0.75,0.9)), aes(x = class, y = prob, fill = class)) +
    geom_boxplot() +
    facet_grid(level~nLoci) +
    labs(x = "Genotype Frequency Class", y = "Efficiency") +
    scale_fill_brewer(palette = "Dark2") +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position = "none", strip.background = element_rect(, colour = "black", fill = "white"),
          strip.text.x = element_text(colour = "black"), strip.text.y = element_text(colour = "black")) +
    coord_cartesian(ylim = c(min(ProbOutput$prob), max(ProbOutput$prob)))
  
  
  
  # Power Boxplot - by hybrid class, and critical posterior probability
  Plot_3 <-
    ggplot2::ggplot(filter(performance_merge,level %in% c(0.5,0.75,0.9)), aes(x = class, y = performance, fill = class)) +
    geom_boxplot() +
    facet_grid(level~nLoci.y) +
    labs(x = "Genotype Frequency Class", y = "Power") +
    scale_fill_brewer(palette = "Dark2") +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position = "none", , legend.key = element_blank(),
          strip.background = element_rect(, colour = "black", fill = "white"),
          strip.text.x = element_text(colour = "black"), strip.text.y = element_text(colour = "black")) +
    coord_cartesian(ylim = c(min(performance_merge$performance), max(performance_merge$performance)))
  
  
  
  
  # Accuracy Boxplot - by pure classes and all hybrids, and critical posterior probability
  Plot_4 <-
    ggplot2::ggplot(filter(AccuracyDataBoxPlot,pofz %in% c(0.5,0.75,0.9)), aes(x = max.class2, y = means, fill = max.class2)) +
    geom_boxplot() +
    facet_grid(pofz~nLoci) +
    labs(x = "Genotype Frequency Class", y = "Accuracy") +
    scale_fill_brewer(palette = "Dark2") +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position = "none", , legend.key = element_blank(),
          strip.background = element_rect(, colour = "black", fill = "white"),
          strip.text.x = element_text(colour = "black"), strip.text.y = element_text(colour = "black")) +
    coord_cartesian(ylim = c(min(AccuracyDataBoxPlot$means), max(AccuracyDataBoxPlot$means)))
  
  
  
  # Efficiency Boxplot - by pure classes and all hybrids, and critical posterior probability
  Plot_5 <-
    ggplot2::ggplot(filter(ProbOutput2,level %in% c(0.5,0.75,0.9)), aes(x = class, y = prob, fill = class)) +
    geom_boxplot() +
    facet_grid(level~nLoci) +
    labs(x = "Genotype Frequency Class", y = "Efficiency") +
    scale_fill_brewer(palette = "Dark2") +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position = "none", strip.background = element_rect(, colour = "black", fill = "white"),
          strip.text.x = element_text(colour = "black"), strip.text.y = element_text(colour = "black")) +
    coord_cartesian(ylim = c(min(ProbOutput2$prob), max(ProbOutput2$prob)))
  
  
  
  # Power Boxplot - by pure classes and all hybrids, and critical posterior probability
  Plot_6 <-
    ggplot2::ggplot(filter(performance_merge_Hyb,pofz %in% c(0.5,0.75,0.9)), aes(x = class, y = performance, fill = class)) +
    geom_boxplot() +
    facet_grid(level~nLoci.y) +
    labs(x = "Genotype Frequency Class", y = "Power") +
    scale_fill_brewer(palette = "Dark2") +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position = "none", , legend.key = element_blank(),
          strip.background = element_rect(, colour = "black", fill = "white"),
          strip.text.x = element_text(colour = "black"), strip.text.y = element_text(colour = "black")) +
    coord_cartesian(ylim = c(min(performance_merge_Hyb$performance), max(performance_merge_Hyb$performance)))
  
  
  
  # Accuracy Lineplot - by hybrid class, faceted by panel size
  Plot_7 <-
    ggplot2::ggplot(SummaryAccuracy) +
    geom_line(aes(x = pofz, y = mean, colour = max.class), lwd = 1.25) +
    geom_line(aes(y = mean+sd, x = pofz, colour = max.class), linetype = 2) +
    geom_line(aes(y = mean-sd, x = pofz, colour = max.class), linetype = 2) +
    facet_wrap(~nLoci, ncol = 3) +
    theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position="bottom", legend.key = element_blank(),
          strip.background = element_rect(colour = "black", fill = "white")) +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Critical Posterior Probability Threshold", y = expression("Accuracy "%+-%"sd"), col="Genotype Frequency Class") +
    coord_cartesian(ylim = c(get.y.min.AccuracyLine, 1))
  
  
  
  # Efficiency Lineplot - by hybrid class, faceted by panel size
  Plot_8 <-
    ggplot2::ggplot(FinalData, aes(x = level, y = mprob, col = class)) +
    geom_line(lwd = 1.25) +
    geom_line(aes(y = mprob + sdprob, x = level, colour = class), linetype = 2) +
    geom_line(aes(y = mprob - sdprob, x = level, colour = class), linetype = 2) +
    facet_wrap(~nLoci, ncol = 3) +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position="bottom", legend.key = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"), text = element_text(colour = "black")) +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Critical Posterior Probability Threshold", y = expression("Efficiency "%+-%"sd"), col = "Genotype Frequency Class") +
    coord_cartesian(ylim = c(min(FinalData$mprob - FinalData$sdprob), 1))
  
  
  
  # Power Lineplot - by hybrid class, faceted by panel size
  Plot_9 <-
    ggplot2::ggplot(Final_Power, aes(x = level, y = meanPower, col = class)) +
    geom_line(lwd = 1.25) +
    geom_line(aes(y = meanPower + sdPower, x = level, colour = class), linetype = 2) +
    geom_line(aes(y = meanPower - sdPower, x = level, colour = class), linetype = 2) +
    facet_wrap(~nLoci, ncol = 3) +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position="bottom", legend.key = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"), text = element_text(colour = "black")) +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Critical Posterior Probability Threshold", y = expression("Power "%+-%"sd"), col = "Genotype Frequency Class") +
    coord_cartesian(ylim = c(min(Final_Power$meanPower - Final_Power$sdPower, na.rm = TRUE), 1))
  
  
  min(Final_Power[!is.na(Final_Power$meanPower) & !is.na(Final_Power$sdPower),]$meanPower - Final_Power[!is.na(Final_Power$meanPower) & !is.na(Final_Power$sdPower),]$sdPower)
  
  # Accuracy Lineplot - by pure classes and all hybrids, faceted by panel size
  Plot_10 <-
    ggplot2::ggplot(SummaryAccuracy2) +
    geom_line(aes(x = pofz, y = mean, colour = class), lwd = 1.25) +
    geom_line(aes(y = mean+sd, x = pofz, colour = class), linetype = 2) +
    geom_line(aes(y = mean-sd, x = pofz, colour = class), linetype = 2) +
    facet_wrap(~nLoci, ncol = 3) +
    theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position="bottom", legend.key = element_blank(),
          strip.background = element_rect(colour = "black", fill = "white")) +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Critical Posterior Probability Threshold", y = expression("Accuracy "%+-%"sd"), col="Genotype Frequency Class") +
    coord_cartesian(ylim = c(get.y.min.AccuracyLine, 1))
  
  
  
  # Efficiency Lineplot - by pure classes and all hybrids, faceted by panel size
  Plot_11 <-
    ggplot2::ggplot(FinalData2) +
    geom_line(aes(x = level, y = mprob, col = class), lwd = 1.25) +
    geom_line(aes(x = level, y = mprob + sdprob, col = class), lty = 2) +
    geom_line(aes(x = level, y = mprob - sdprob, col = class), lty = 2) +
    facet_wrap(~nLoci, ncol = 3)+
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position = "bottom", legend.key = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"), text = element_text(colour = "black")) +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Critical Posterior Probability Threshold", y = expression("Efficiency "%+-%"sd"), col = "Genotype Frequency Class") +
    coord_cartesian(ylim = c(min(FinalData2$mprob - FinalData2$sdprob, na.rm = TRUE), 1))
  
  
  
  # Power Lineplot - by pure classes and all hybrids, faceted by panel size
  Plot_12 <-
    ggplot2::ggplot(Final_Power2, aes(x = level, y = meanPower, col = class)) +
    geom_line(lwd = 1.25) +
    geom_line(aes(y = meanPower + sdPower, x = level, colour = class), linetype = 2) +
    geom_line(aes(y = meanPower - sdPower, x = level, colour = class), linetype = 2) +
    facet_wrap(~nLoci, ncol = 3) +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position="bottom", legend.key = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"), text = element_text(colour = "black")) +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Critical Posterior Probability Threshold", y = expression("Power "%+-%"sd"), col = "Genotype Frequency Class") +
    coord_cartesian(ylim = c(min(Final_Power2$meanPower - Final_Power2$sdPower, na.rm = TRUE), 1))
  
  
  
  # Accuracy Lineplot - by hybrid class, faceted by panel size and \"hybrid grouping\"
  Plot_13 <-
    ggplot2::ggplot(SummaryAccuracy) +
    geom_line(aes(x = pofz, y = mean, colour = max.class), lwd = 1.25) +
    geom_line(aes(y = mean+sd, x = pofz, colour = max.class), linetype = 2) +
    geom_line(aes(y = mean-sd, x = pofz, colour = max.class), linetype = 2) +
    facet_grid( Group ~ nLoci, scales = "free_y") +
    theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position="bottom", legend.key = element_blank(),
          strip.background = element_rect(colour = "black", fill = "white")) +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Critical Posterior Probability Threshold", y = expression("Accuracy "%+-%"sd"), col="Genotype Frequency Class") +
    coord_cartesian(ylim = c(get.y.min.AccuracyLine, 1))
  
  
  
  # Efficiency Lineplot - by hybrid class, faceted by panel size and \"hybrid grouping\"
  Plot_14 <-
    ggplot2::ggplot(FinalData, aes(x = level, y = mprob, col = class)) +
    geom_line(lwd = 1.25) +
    geom_line(aes(y = mprob + sdprob, x = level, colour = class), linetype = 2) +
    geom_line(aes(y = mprob - sdprob, x = level, colour = class), linetype = 2) +
    facet_grid( group ~ nLoci, scales = "free_y") +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position="bottom", legend.key = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"), text = element_text(colour = "black")) +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Critical Posterior Probability Threshold", y = expression("Efficiency "%+-%"sd"), col = "Genotype Frequency Class") +
    coord_cartesian(ylim = c(min(FinalData$mprob - FinalData$sdprob, na.rm = TRUE), 1))
  
  
  
  # Power Lineplot - by hybrid class, faceted by panel size and \"hybrid grouping\"
  Plot_15 <-
    ggplot2::ggplot(Final_Power, aes(x = level, y = meanPower, col = class)) +
    geom_line(lwd = 1.25) +
    geom_line(aes(y = meanPower + sdPower, x = level, colour = class), linetype = 2) +
    geom_line(aes(y = meanPower - sdPower, x = level, colour = class), linetype = 2) +
    facet_grid( group ~ nLoci, scales = "free_y") +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position="bottom", legend.key = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"), text = element_text(colour = "black")) +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Critical Posterior Probability Threshold", y = expression("Power "%+-%"sd"), col = "Genotype Frequency Class") +
    coord_cartesian(ylim = c(min(Final_Power$meanPower - Final_Power$sdPower, na.rm = TRUE), 1))
  
  min(Final_Power$meanPower - Final_Power$sdPower, na.rm = TRUE)
  
  # Accuracy Lineplot - by panel size, faceted by hybrid class"
  Plot_16 <-
    ggplot2::ggplot(SummaryAccuracy) +
    geom_line(aes(x = pofz, y = mean, colour = nLoci), lwd = 1.25) +
    geom_line(aes(y = mean - sd, x = pofz, colour = nLoci), linetype = 2) +
    geom_line(aes(y = mean + sd, x = pofz, colour = nLoci), linetype = 2) +
    facet_wrap(~max.class, nrow = 3) +
    theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position="bottom", legend.key = element_blank(),
          strip.background = element_rect(colour = "black", fill = "white")) +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Critical Posterior Probability Threshold", y = expression("Accuracy "%+-%"sd"), col="Panel Size (Loci)") +
    coord_cartesian(ylim = c(get.y.min.AccuracyLine, 1))
  
  
  
  # Efficiency Lineplot - by panel size, faceted by hybrid class"
  Plot_17 <-
    ggplot2::ggplot(data = FinalData) +
    geom_line(aes(x = level, y = mprob, col = factor(nLoci)), lwd = 1.25) +
    geom_line(aes(x = level, y = mprob + sdprob, col = factor(nLoci)), lty = 2) +
    geom_line(aes(x = level, y = mprob - sdprob, col = factor(nLoci)), lty = 2) +
    facet_wrap(~class, nrow = 3) +
    scale_color_brewer(palette = "Dark2") +
    theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"),
          panel.grid.major = element_line(colour = "grey90"), legend.position = "bottom", legend.key = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"), text = element_text(colour = "black")) +
    labs(x = "Critical Posterior Probability Threshold", y = expression("Efficiency "%+-%"sd"), col = "Panel Size (Loci)")  +
    coord_cartesian(ylim = c(min(FinalData$mprob - FinalData$sdprob, na.rm = TRUE), 1))
  
  
  
  
  # Power Lineplot - by panel size, faceted by hybrid class
  Plot_18 <-
    ggplot2::ggplot(Final_Power, aes(x = level, y = meanPower, colour = as.factor(nLoci))) +
    geom_line(lwd = 1.25) +
    geom_line(aes(y = meanPower + sdPower, x = level, colour = as.factor(nLoci)), linetype = 2) +
    geom_line(aes(y = meanPower - sdPower, x = level, colour = as.factor(nLoci)), linetype = 2) +
    facet_wrap(~class, nrow = 3) +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position="bottom", legend.key = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"), text = element_text(colour = "black")) +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Critical Posterior Probability Threshold", y = expression("Power "%+-%"sd"), col = "Panel Size (Loci)") +
    coord_cartesian(ylim = c(min(Final_Power$meanPower - Final_Power$sdPower, na.rm = TRUE), 1))
  
  
  
  # Accuracy Lineplot - by panel size, faceted by pure classes and all hybrids
  Plot_19 <-
    ggplot2::ggplot(ComboHybridAccuracy) +
    geom_line(aes(x = pofz, y = mprob, colour = nLoci), lwd = 1.25) +
    geom_line(aes(y = mprob - sdprob, x = pofz, colour = nLoci), linetype = 2) +
    geom_line(aes(y = mprob + sdprob, x = pofz, colour = nLoci), linetype = 2) +
    facet_wrap(~class, ncol = 3) +
    theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position="bottom", legend.key = element_blank(),
          strip.background = element_rect(colour = "black", fill = "white")) +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Critical Posterior Probability Threshold", y = expression("Accuracy "%+-%"sd"), col="Panel Size (Loci)") +
    coord_cartesian(ylim = c(get.y.min.AccuracyLine, 1))
  
  
  
  # Efficiency Lineplot - by panel size, faceted by pure classes and all hybrids
  Plot_20 <-
    ggplot2::ggplot(data = FinalData2) +
    geom_line(aes(x = level, y = mprob, col = factor(nLoci)), lwd = 1.25) +
    geom_line(aes(x = level, y = mprob + sdprob,col = factor(nLoci)), lty = 2) +
    geom_line(aes(x = level, y = mprob - sdprob, col = factor(nLoci)), lty = 2) +
    facet_grid(~class) +
    scale_color_brewer(palette = "Dark2") +
    theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"),
          panel.grid.major = element_line(colour = "grey90"), legend.position = "bottom", legend.key = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"), text = element_text(colour = "black")) +
    labs(x = "Critical Posterior Probability Threshold", y = expression("Efficiency "%+-%"sd"), col = "Panel Size (Loci)") +
    coord_cartesian(ylim = c(min(FinalData2$mprob, na.rm = TRUE), 1))
  
  
  
  # Power Lineplot - by panel size, faceted by pure classes and all hybrids
  Plot_21 <-
    ggplot2::ggplot(Final_Power2, aes(x = level, y = meanPower, colour = as.factor(nLoci))) +
    geom_line(lwd = 1.25) +
    geom_line(aes(y = meanPower + sdPower, x = level, colour = as.factor(nLoci)), linetype = 2) +
    geom_line(aes(y = meanPower - sdPower, x = level, colour = as.factor(nLoci)), linetype = 2) +
    facet_wrap(~class, ncol = 3) +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position="bottom", legend.key = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"), text = element_text(colour = "black")) +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Critical Posterior Probability Threshold", y = expression("Power "%+-%"sd"), col = "Panel Size (Loci)") +
    coord_cartesian(ylim = c(min(Final_Power$meanPower, na.rm = TRUE), 1))
  
  
  
  # Accuracy Dotplot - by hybrid class, faceted by critical posterior probability %in% Thresholds
  Plot_22 <-
    ggplot2::ggplot(dplyr::filter(SummaryAccuracy,pofz %in% Thresholds), aes(x = factor(nLoci), y = mean, col = max.class, group = max.class)) +
    geom_point(size = 2.5, position = position_dodge(0.5)) + geom_path(lwd = 0.9, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(0.5), width = 0.5) +
    facet_grid(~pofz) +
    labs(x="Panel Size (Loci)",y=expression("Accuracy "%+-%"sd"), col = "Genotype Frequency Class", group = "") +
    scale_color_brewer(palette = "Dark2") +
    theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"),
          panel.grid.major = element_line(colour = "grey90"),
          legend.position = "bottom", legend.key = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"), text = element_text(colour = "black")) +
    coord_cartesian(ylim = c(get.y.min.AccuracyThreshold, 1))
  
  
  # Efficiency Dotplot - by hybrid class, faceted by critical posterior probability %in% Thresholds
  Plot_23 <-
    ggplot2::ggplot(dplyr::filter(FinalData, level %in% Thresholds), aes(x = factor(nLoci), y = mprob, col = class, group = class)) +
    geom_point(size = 2.5, position = position_dodge(0.5)) +
    geom_path(lwd = 0.9, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymin = mprob - sdprob, ymax = mprob + sdprob), width = 0.5, position = position_dodge(0.5)) +
    facet_grid(~level) +
    labs(x = "Panel Size (Loci)", y = expression("Efficiency "%+-%"sd"), col = "Genotype Frequency Class", group = "") +
    scale_color_brewer(palette = "Dark2") +
    theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"),
          panel.grid.major = element_line(colour = "grey90"),
          legend.position = "bottom", legend.key = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"), text = element_text(colour = "black")) +
    coord_cartesian(ylim = c(min(FinalData$mprob - FinalData$sdprob, na.rm = TRUE), 1))
  
  
  
  
  # Power Dotplot - by hybrid class, faceted by critical posterior probability %in% Thresholds
  Plot_24 <-
    ggplot2::ggplot(dplyr::filter(Final_Power,level %in% Thresholds), aes(x = factor(nLoci), y = meanPower, col = class, group = class)) +
    geom_point(size = 2.5, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymin = meanPower - sdPower, ymax = meanPower + sdPower), position = position_dodge(0.5), width = 0.5) +
    facet_grid(~level) +
    labs(x="Panel Size (Loci)",y=expression("Power "%+-%"sd"), col = "Genotype Frequency Class", group = "") +
    scale_color_brewer(palette = "Dark2") +
    theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"),
          panel.grid.major = element_line(colour = "grey90"),
          legend.position = "bottom", legend.key = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"), text = element_text(colour = "black")) +
    coord_cartesian(ylim = c(min(Final_Power$meanPower - Final_Power$sdPower, na.rm = TRUE), 1))
  
  
  
  # Accuracy Dotplot - by pure classes and all hybrids, faceted by critical posterior probability %in% Thresholds
  Plot_25 <-
    ggplot2::ggplot(dplyr::filter(ComboHybridAccuracy, pofz %in% Thresholds), aes(x = factor(nLoci), y = mprob, col = class, group = class)) +
    geom_point(size = 2.5, position = position_dodge(0.5)) +
    geom_path(lwd = 0.9, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymin = (mprob - sdprob), ymax = (mprob + sdprob)), width = 0.5, position = position_dodge(0.5)) +
    facet_grid(~pofz) +
    labs(x = "Panel Size (Loci)", y = expression("Accuracy "%+-%"sd"), col = "Genotype Frequency Class", group = "") +
    scale_color_brewer(palette = "Dark2") +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position = "bottom", legend.key = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"),
          text = element_text(colour = "black")) +
    coord_cartesian(ylim = c(min(ComboHybridAccuracy$mprob - ComboHybridAccuracy$sdprob, na.rm = TRUE), 1))
  
  
  
  # Efficiency Dotplot - by pure classes and all hybrids, faceted by critical posterior probability %in% Thresholds
  Plot_26 <-
    ggplot2::ggplot(dplyr::filter(FinalData2, level %in% Thresholds), aes(x = factor(nLoci), y = mprob, col = class, group = class)) +
    geom_point(size = 2.5, position = position_dodge(0.5)) + geom_path(lwd = 0.9, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymin = mprob - sdprob, ymax = mprob + sdprob), width = 0.5, position = position_dodge(0.5)) +
    facet_grid(~level) +
    labs(x = "Panel Size (Loci)", y = expression("Efficiency "%+-%"sd"), col = "Genotype Frequency Class", group = "")+
    scale_color_brewer(palette = "Dark2") +
    theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"),
          panel.grid.major = element_line(colour = "grey90"),
          legend.position = "bottom", legend.key = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"), text = element_text(colour = "black")) +
    coord_cartesian(ylim = c(min(FinalData$mprob - FinalData$sdprob, na.rm = TRUE), 1))
  
  
  
  # Power Dotplot - by pure classes and all hybrids, faceted by critical posterior probability %in% Thresholds
  Plot_27 <-
    ggplot2::ggplot(dplyr::filter(Final_Power2,level %in% Thresholds), aes(x = factor(nLoci), y = meanPower, col = class, group = class)) +
    geom_point(size = 2.5, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymin = meanPower - sdPower, ymax = meanPower + sdPower), position = position_dodge(0.5), width = 0.5) +
    facet_grid(~level) +
    labs(x="Panel Size (Loci)",y=expression("Power "%+-%"sd"), col = "Genotype Frequency Class", group = "") +
    scale_color_brewer(palette = "Dark2") +
    theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"),
          panel.grid.major = element_line(colour = "grey90"),
          legend.position = "bottom", legend.key = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"), text = element_text(colour = "black")) +
    coord_cartesian(ylim = c(min(Final_Power2$meanPower - Final_Power2$sdPower, na.rm = TRUE), 1))
  
  
  
  
  # Mean Posterior Probability of Assignment - per simulation, faceted by hybrid class
  Plot_30 <-
    ggplot2::ggplot(boxdata, aes(x = nLoci, y = value, fill = sim, group = sim)) +
    geom_boxplot(alpha = 0.8, outlier.size = 0) +
    facet_wrap(~class, nrow = 3, scales = "free_y") +
    labs(y = "Posterior Probability", x = "Panel Size (Loci)") +
    scale_fill_manual(values = c("grey75", "grey75", "grey75")) +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position = "none", strip.background = element_rect(colour = "black", fill = "white"),
          text = element_text(colour = "black"))
  
  
  
  # Mean Posterior Probability of Assignment - per simulation, faceted by  pure classes and all hybrids
  Plot_31 <-
    ggplot2::ggplot(sim_means2, aes(x = factor(nLoci), y = hybrid, fill = sim)) +
    geom_boxplot(alpha = 0.8, outlier.size = 0) +
    facet_wrap(~hclass, nrow = 3, scales = "free_y") +
    labs(y = "Posterior Probability", x = "Panel Size (Loci)") +
    scale_fill_manual(values = c("grey75", "grey75", "grey75")) +
    theme(panel.background = element_rect(fill = "white", colour = "black"), plot.background = element_rect(colour = "white"), panel.grid.major = element_line(colour = "grey90"),
          legend.position = "none", strip.background = element_rect(colour = "black", fill = "white"), text = element_text(colour = "black"))
  
  
  
  
  HybridPower_Plots <- list(Plot_1, Plot_2, Plot_3, Plot_4, Plot_5, Plot_6, Plot_7, Plot_8, Plot_9, Plot_10, Plot_11, Plot_12, Plot_13, Plot_14, Plot_15, Plot_16, Plot_17, Plot_18, Plot_19, Plot_20, Plot_21, Plot_22, Plot_23, Plot_24, Plot_25, Plot_26, Plot_27,   Plot_30, Plot_31)
  
  
  
  
  
  ########################
  
  ############################################################
  #### Get the data that was plotted for people and stuff ####
  ############################################################
  
  
  
  Plot_1_Data <- dplyr::filter(AccuracyData,pofz %in% c(0.5,0.75,0.9))
  Plot_1_Data <- data.frame(x = Plot_1_Data$max.class, y = Plot_1_Data$means, fill = Plot_1_Data$max.class, facet.1 = Plot_1_Data$pofz, facet.2 = Plot_1_Data$nLoci)
  
  Plot_2_Data <- dplyr::filter(ProbOutput,level %in% c(0.5,0.75,0.9))
  Plot_2_Data = data.frame(x = Plot_2_Data$class, y = Plot_2_Data$prob, fill = Plot_2_Data$class, facet.1 = Plot_2_Data$nLoci)
  
  Plot_3_Data <- dplyr::filter(performance_merge,level %in% c(0.5,0.75,0.9))
  Plot_3_Data <- data.frame(x = Plot_3_Data$class, y = Plot_3_Data$performance, fill = Plot_3_Data$class, facet.1 = Plot_3_Data$level, facet.2 = Plot_3_Data$nLoci.y)
  
  Plot_4_Data <- dplyr::filter(AccuracyDataBoxPlot,pofz %in% c(0.5,0.75,0.9))
  Plot_4_Data <- data.frame(x = Plot_4_Data$max.class2, y = Plot_4_Data$means, fill = Plot_4_Data$max.class2, facet.1 = Plot_4_Data$pofz, facet.2 = Plot_4_Data$nLoci)
  
  Plot_5_Data <- dplyr::filter(ProbOutput2,level %in% c(0.5,0.75,0.9))
  Plot_5_Data <- data.frame(x = Plot_5_Data$class, y = Plot_5_Data$prob, fill = Plot_5_Data$class, facet.1 = Plot_5_Data$level, facet.2 = Plot_5_Data$nLoci)
  
  Plot_6_Data <- dplyr::filter(performance_merge_Hyb,pofz %in% c(0.5,0.75,0.9))
  Plot_6_Data <- data.frame(x = Plot_6_Data$class, y = Plot_6_Data$performance, fill = Plot_6_Data$class, facet.1 = Plot_6_Data$level, facet.2 = Plot_6_Data$nLoci.y)
  
  Plot_7_Data <- SummaryAccuracy
  Plot_7_Data <- data.frame(x = Plot_7_Data$pofz, y = Plot_7_Data$mean, sd.upper = (Plot_7_Data$mean + Plot_7_Data$sd), sd.lower = (Plot_7_Data$mean - Plot_7_Data$sd), colour = Plot_7_Data$max.class, facet = Plot_7_Data$nLoci)
  
  Plot_8_Data <- FinalData
  Plot_8_Data <- data.frame(x = Plot_8_Data$level, y = Plot_8_Data$mprob, sd.upper = (Plot_8_Data$mprob + Plot_8_Data$sdprob), sd.lower = (Plot_8_Data$mprob - Plot_8_Data$sd), colour = Plot_8_Data$class, facet = as.factor(Plot_8_Data$nLoci))
  
  Plot_9_Data <- Final_Power
  Plot_9_Data <- data.frame(x = Plot_9_Data$level, y = Plot_9_Data$meanPower, sd.upper = (Plot_9_Data$meanPower + Plot_9_Data$sdPower), sd.lower = (Plot_9_Data$meanPower - Plot_9_Data$sdPower), colour = Plot_9_Data$class, facet = as.factor(Plot_9_Data$nLoci))
  
  Plot_10_Data <- SummaryAccuracy2
  Plot_10_Data <- data.frame(x = Plot_10_Data$pofz, y = Plot_10_Data$mean, sd.upper = (Plot_10_Data$mean + Plot_10_Data$sd), sd.lower = (Plot_10_Data$mean - Plot_10_Data$sd), colour = Plot_10_Data$class, facet = as.factor(Plot_10_Data$nLoci))
  
  Plot_11_Data <- FinalData2
  Plot_11_Data <- data.frame(x = Plot_11_Data$level, y = Plot_11_Data$mprob, sd.upper = (Plot_11_Data$mprob + Plot_11_Data$sdprob), sd.lower = (Plot_11_Data$mprob - Plot_11_Data$sdprob), colour = Plot_11_Data$class, facet = as.factor(Plot_11_Data$nLoci))
  
  Plot_12_Data <- Final_Power2
  Plot_12_Data <- data.frame(x = Plot_12_Data$level, y = Plot_12_Data$meanPower, sd.upper = (Plot_12_Data$meanPower + Plot_12_Data$sdPower), sd.lower = (Plot_12_Data$meanPower - Plot_12_Data$sdPower), colour = Plot_12_Data$class, facet = as.factor(Plot_12_Data$nLoci))
  
  Plot_13_Data <- SummaryAccuracy
  Plot_13_Data <- data.frame(x = Plot_13_Data$pofz, y = Plot_13_Data$mean, sd.upper = (Plot_13_Data$mean + Plot_13_Data$sd), sd.lower = (Plot_13_Data$mean - Plot_13_Data$sd), colour = Plot_13_Data$max.class, facet.1 = Plot_13_Data$Group, facet.2 = as.factor(Plot_13_Data$nLoci))
  
  Plot_14_Data <- FinalData
  Plot_14_Data <- data.frame(x = Plot_14_Data$level, y = Plot_14_Data$mprob, sd.upper = (Plot_14_Data$mprob + Plot_14_Data$sdprob), sd.lower = (Plot_14_Data$mprob - Plot_14_Data$sdprob), colour = Plot_14_Data$class, facet.1 = Plot_14_Data$group, facet.2 = as.factor(Plot_14_Data$nLoci))
  
  Plot_15_Data <- Final_Power
  Plot_15_Data <- data.frame(x = Plot_15_Data$level, y = Plot_15_Data$meanPower, sd.upper = (Plot_15_Data$meanPower + Plot_15_Data$sdPower), sd.lower = (Plot_15_Data$meanPower - Plot_15_Data$sdPower), colour = Plot_15_Data$class, facet.1 = Plot_15_Data$group, facet.2 = as.factor(Plot_15_Data$nLoci))
  
  Plot_16_Data <- SummaryAccuracy
  Plot_16_Data <- data.frame(x = Plot_16_Data$pofz, y = Plot_16_Data$mean, sd.upper = (Plot_16_Data$mean + Plot_16_Data$sd), sd.lower = (Plot_16_Data$mean - Plot_16_Data$sd), colour = as.factor(Plot_16_Data$nLoci), facet = Plot_16_Data$max.class)
  
  Plot_17_Data <- FinalData
  Plot_17_Data <- data.frame(x = Plot_17_Data$level, y = Plot_17_Data$mprob, sd.upper = (Plot_17_Data$mprob + Plot_17_Data$sdprob), sd.lower = (Plot_17_Data$mprob - Plot_17_Data$sdprob), colour = as.factor(Plot_17_Data$nLoci), facet = Plot_17_Data$class)
  
  Plot_18_Data <- Final_Power
  Plot_18_Data <- data.frame(x = Plot_18_Data$level, y = Plot_18_Data$meanPower, sd.upper = (Plot_18_Data$meanPower + Plot_18_Data$sdPower), sd.lower = (Plot_18_Data$meanPower - Plot_18_Data$sdPower), colour = as.factor(Plot_18_Data$nLoci), facet = Plot_18_Data$class)
  
  Plot_19_Data <- ComboHybridAccuracy
  Plot_19_Data <- data.frame(x = Plot_19_Data$pofz, y = Plot_19_Data$mprob, sd.upper = (Plot_19_Data$mprob + Plot_19_Data$sdprob), sd.lower = (Plot_19_Data$mprob - Plot_19_Data$sdprob), colour = as.factor(Plot_19_Data$nLoci), facet = Plot_19_Data$class)
  
  Plot_20_Data <- FinalData2
  Plot_20_Data <- data.frame(x = Plot_20_Data$level, y = Plot_20_Data$mprob, sd.upper = (Plot_20_Data$mprob + Plot_20_Data$sdprob), sd.lower = (Plot_20_Data$mprob - Plot_20_Data$sdprob), colour = as.factor(Plot_20_Data$nLoci), facet = Plot_20_Data$class)
  
  Plot_21_Data <- Final_Power2
  Plot_21_Data <- data.frame(x = Plot_21_Data$level, y = Plot_21_Data$meanPower, sd.upper = (Plot_21_Data$meanPower + Plot_21_Data$sdPower), sd.lower = (Plot_21_Data$meanPower - Plot_21_Data$sdPower), colour = as.factor(Plot_21_Data$nLoci), facet = Plot_21_Data$class)
  
  Plot_22_Data <- dplyr::filter(SummaryAccuracy,pofz %in% Thresholds)
  Plot_22_Data <- data.frame(x = as.factor(Plot_22_Data$nLoci), y = Plot_22_Data$mean, sd.upper = (Plot_22_Data$mean + Plot_22_Data$sd), sd.lower = (Plot_22_Data$mean - Plot_22_Data$sd), colour = Plot_22_Data$max.class, group = Plot_22_Data$max.class, facet = Plot_22_Data$pofz)
  
  Plot_23_Data <- dplyr::filter(FinalData, level %in% Thresholds)
  Plot_23_Data <- data.frame(x = as.factor(Plot_23_Data$nLoci), y = Plot_23_Data$mprob, sd.upper = (Plot_23_Data$mprob + Plot_23_Data$sdprob), sd.lower = (Plot_23_Data$mprob - Plot_23_Data$sdprob), colour = Plot_23_Data$class, group = Plot_23_Data$class, facet = Plot_23_Data$level)
  
  Plot_24_Data <- dplyr::filter(Final_Power,level %in% Thresholds)
  Plot_24_Data <- data.frame(x = as.factor(Plot_24_Data$nLoci), y = Plot_24_Data$meanPower, sd.upper = (Plot_24_Data$meanPower + Plot_24_Data$sdPower), sd.lower = (Plot_24_Data$meanPower - Plot_24_Data$sdPower), colour = Plot_24_Data$class, group = Plot_24_Data$class, facet = Plot_24_Data$level)
  
  Plot_25_Data <- dplyr::filter(ComboHybridAccuracy, pofz %in% Thresholds)
  Plot_25_Data <- data.frame(x = as.factor(Plot_25_Data$nLoci), y = Plot_25_Data$mprob, sd.upper = (Plot_25_Data$mprob + Plot_25_Data$sdprob), sd.lower = (Plot_25_Data$mprob - Plot_25_Data$sdprob), colour = Plot_25_Data$class, group = Plot_25_Data$class, facet = Plot_25_Data$pofz)
  
  Plot_26_Data <- dplyr::filter(FinalData2, level %in% Thresholds)
  Plot_26_Data <- data.frame(x = as.factor(Plot_26_Data$nLoci), y = Plot_26_Data$mprob, sd.upper = (Plot_26_Data$mprob + Plot_26_Data$sdprob), sd.lower = (Plot_26_Data$mprob - Plot_26_Data$sdprob), colour = Plot_26_Data$class, group = Plot_26_Data$class, facet = Plot_26_Data$level)
  
  Plot_27_Data <- dplyr::filter(Final_Power2,level %in% Thresholds)
  Plot_27_Data <- data.frame(x = as.factor(Plot_27_Data$nLoci), y = Plot_27_Data$meanPower, sd.upper = (Plot_27_Data$meanPower + Plot_27_Data$sdPower), sd.lower = (Plot_27_Data$meanPower - Plot_27_Data$sdPower), colour = Plot_27_Data$class, group = Plot_27_Data$class, facet = Plot_27_Data$level)
  
  
  Plot_30_Data <- boxdata
  Plot_30_Data <- data.frame(x = as.factor(Plot_30_Data$nLoci), y = Plot_30_Data$value, group = Plot_30_Data$sim, facet = Plot_30_Data$class)
  
  Plot_31_Data <- sim_means2
  Plot_31_Data <- data.frame(x = as.factor(Plot_31_Data$nLoci), y = Plot_31_Data$hybrid, fill = Plot_31_Data$sim, facet = Plot_31_Data$hclass)
  
  
  
  
  HybridPower_Data <- list(Plot_1_Data, Plot_2_Data, Plot_3_Data, Plot_4_Data, Plot_5_Data, Plot_6_Data, Plot_7_Data, Plot_8_Data, Plot_9_Data, Plot_10_Data, Plot_11_Data, Plot_12_Data, Plot_13_Data, Plot_14_Data, Plot_15_Data, Plot_16_Data, Plot_17_Data, Plot_18_Data, Plot_19_Data, Plot_20_Data, Plot_21_Data, Plot_22_Data, Plot_23_Data, Plot_24_Data, Plot_25_Data, Plot_26_Data, Plot_27_Data, Plot_30_Data, Plot_31_Data)
  
  
  HybridPower_Output <- list(legend = Plot.Legend, plots = HybridPower_Plots, data = HybridPower_Data)
  
  ########################
  
  ####################################
  #### Save the sweet sweet plots ####
  ####################################
  
  writeLines("
                      I'm saving your plots for you over here
                      ")
  
  
  writeLines("
                      I'm savin' the data too
                      ")
  
  ### Save data and plots is save_output == TRUE
  if(save_output == TRUE){
    ## save the plot legends as well
    write.csv(x = Plot.Legend, file = paste0(dir, "Figures and Data/pdf/Plot_Legends.csv"),row.names = FALSE, quote = FALSE)
    write.csv(x = Plot.Legend, file = paste0(dir, "Figures and Data/jpg/Plot_Legends.csv"),row.names = FALSE, quote = FALSE)
    for(j in 1:length(Plot.Names)){
      # Save .pdf version of plots
      ggsave(paste0(dir, "Figures and Data/pdf/", Plot.Names[j], ".pdf"), HybridPower_Plots[[j]], height = 10, width = 10)
      # Save .jpg version of plots
      ggsave(paste0(dir, "Figures and Data/jpg/", Plot.Names[j], ".jpg"), HybridPower_Plots[[j]], height = 10, width = 10)
      # Save plot data
      write.csv(x = HybridPower_Data[[j]], file = paste0(dir, "Figures and Data/data/", Plot.Names[j], ".csv"), row.names = FALSE, quote = FALSE)
    }
  }
  
  ## Return plots and data
  if(return_workspace == TRUE){
    return(HybridPower_Output)
  }
  
  
  
  
}
