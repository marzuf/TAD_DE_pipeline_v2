startTime <- Sys.time()

options(scipen=100)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script8_name <- "8c_runAllDown"
script_name <- "14i2split_cumulAllDown_limited_AUC_randomTADsShuffle"
cat(paste0(pipScriptDir, "/", script_name, ".R"), "\n")
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))


source("main_settings.R")
#source("run_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(RColorBrewer, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 

registerDoMC(ifelse(SSHFS,2, nCpu)) # loaded from main_settings.R

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

# ADDED 27.11.2018 to check using other files
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)

#****************************************************************************************************************************

# INPUT DATA
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)


#### DRAW THE LEFT PLOT OF THE INITIAL 14E BUT THIS TIME PLOT THE DIFFERENT LINES AND AREA FOR
# -> FIX SIZE
# -> GAUSSIAN SIZE

# !!! WHAT IS DIFFERENT HERE: THE PERMUTATIONS DATA ARE NOT STORED IN A DATAFRAME BUT IN A LIST:
# WHEN SHUFFLING THE GENES, I ALWAYS KEEP THE SAME SET OF TADs BUT NOW THAT I SHUFFLE TAD
# THE NUMBER OF TADs MIGHT BE DIFFERENT (REGIONS OF THE GENOME WHERE NO GENES -> TAD DISCARDED)

# UPDATE SELECT THE GENES ACCORDING TO THE SETTINGS PREPARED IN 0_PREPGENEDATA
initList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_geneList.Rdata"))))
geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))
txt <- paste0(toupper(script_name), "> Start with # genes: ", length(geneList), "/", length(initList), "\n")
printAndLog(txt, pipLogFile)

stopifnot(!any(duplicated(names(geneList))))

gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% geneList,]
geneNbr <- setNames(as.numeric(table(gene2tadDT$region)), names(table(gene2tadDT$region)))

### SET OUTPUT
plotType <- "svg"
myHeight <- ifelse(plotType == "png", 480 , 7)
myWidth <- ifelse(plotType == "png", 600, 10)

# if plotSeparated == TRUE => 1 plot per ratio, otherwise all on the same figure (#x2 plots)
plotSeparated <- F

# "permThresh" the quantile of permutations to take is loaded from main_settings.R

# retrieved from main_settings.R

# stopifnot(any(c(step14_for_randomTADsShuffle  )))
# this version is only implemented for step14_for_randomTADsShuffle !!!
cat("... warning: this version is only implemented for step14_for_randomTADsShuffle !!!\n")
step14_for_randomTADsShuffle <- TRUE
step14_for_randomTADsFix <- FALSE
step14_for_randomTADsGaussian <- FALSE

ratioSizeFix <- percentTADsizeFix
ratioSizeGaussian <- percentTADsizeGaussian 

# for each of the ratioSize -> define the colors
# same colors for gaussian and fix
# but transparency for gaussian (because polygon overlap)

# polygonPermutCol =  rgb(0/255,76/255,153/255, 0.3)
myPalette1 <- brewer.pal(9, "YlOrRd")
names(myPalette1) <- seq(0.1,0.9,by=0.1)

myPalette2 <- rev(brewer.pal(9, "YlGn"))
names(myPalette2) <- seq(1.1,1.9,by=0.1)

myPalette <- c(myPalette1, myPalette2)

fixColors_rgb <- lapply(myPalette,  function(x) as.numeric(col2rgb(x)))

gaussianColors <- lapply(fixColors_rgb, function(x) {
  x <- x/255
  rgb(red = x[1], green=x[2], blue=x[3], alpha=0.3)
})

fixColors <- myPalette

shufflePolygonCol <-  rgb(0/255,76/255,153/255, 0.3)

plotTrueType <- "l"
# "l" to plot the true data as line, "p" to plot them as points

###################################################################################################
# get the direction of up/down
###################################################################################################
# retrieve the direction of up/down
DE_topTable <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_topTable.Rdata"))))
DE_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_geneList.Rdata"))))
exprDT <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_rnaseqDT.Rdata"))))
samp1 <- eval(parse(text=load(paste0(setDir, "/", sample1_file))))
samp2 <- eval(parse(text=load(paste0(setDir, "/", sample2_file))))
stopifnot(all(DE_topTable$genes %in% names(DE_geneList)))
stopifnot(!any(duplicated(names(DE_geneList))))
stopifnot(all(colnames(exprDT) %in% c(samp1, samp2)))
stopifnot(all(samp1 %in% colnames(exprDT)))
stopifnot(all(samp2 %in% colnames(exprDT)))
maxDownGene <- DE_topTable$genes[which.min(DE_topTable$logFC)]
stopifnot(maxDownGene %in% rownames(exprDT))
mean_expr1 <- mean(unlist(c(exprDT[maxDownGene, samp1])), na.rm=T)
mean_expr2 <- mean(unlist(c(exprDT[maxDownGene, samp2])), na.rm=T)

if(mean_expr1 > mean_expr2) {
  subtitDir <- paste0("down: ", toupper(cond1), " > ", toupper(cond2))
} else{
  subtitDir <- paste0("down: ", toupper(cond2), " > ", toupper(cond1))
}

# allDown retrieved from main_settings.R
# allDown_limited <- c("ratioDown", "rescWeighted", "rescWeightedQQ", "prodSignedRatio")

if(! plotSeparated) {
  nColPlot <- 3  # here I plot only with the area
  # nRowPlot <- length(allDown_limited)*2/nColPlot # only 1 plot/row in this version
  nRowPlot <- length(allDown_limited)*3/nColPlot
  outFile <- paste0(curr_outFold, "/allRatios_cumsum_obs_permut.", plotType)
  do.call(plotType, list(outFile, height=myHeight*nRowPlot, width=myWidth*nColPlot))
  par(mfrow=c(nRowPlot, nColPlot))
}

all_ratio_split <- c()

for(curr_ratio_type in allDown_limited) {
  cat(paste0("*** START ", curr_ratio_type, "\n"))
  
  obs_curr_down <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/all_obs_", curr_ratio_type, ".Rdata"))))
  sort_obs_curr_down <- sort(obs_curr_down, decreasing = TRUE)
  
  # FOR ratioDown => plot ratioConcord, departure from 0.5
  if(curr_ratio_type == "ratioDown") {
    my_stat_curr_ratio <- "ratioDown_Concord"
    departureFromValue <- 0.5
    # => Concord, departure 0.5
    # Transform ratioDown -> ratioConcord
    # change so that the ratioDown ranges between 0.5 and 1 (-> e.g. treats 0.1 as 0.9)
    # transf. ratioDown -> ratioConcord
    sort_obs_curr_down <- abs(sort_obs_curr_down - 0.5) + 0.5
  } else if(curr_ratio_type == "rescWeightedQQ" | curr_ratio_type == "rescWeighted" ) {
    my_stat_curr_ratio <- paste0(curr_ratio_type, "_Concord")
    departureFromValue <- 0.5
    # => Concord, departure 0.5
    # Transform rescWeightedQQ -> rescWeightedQQConcord
    # change so that the ratioDown ranges between 0.5 and 1 (-> e.g. treats 0.1 as 0.9)
    # transf. ratioDown -> ratioConcord
    sort_obs_curr_down <- abs(sort_obs_curr_down - 0.5) + 0.5
  } else if(curr_ratio_type == "prodSignedRatio") {
    my_stat_curr_ratio <- "prodSignedRatio"
    departureFromValue <- 0
    # => raw (departure 0)
    # prodSignedRatio -> does not need to be transformed
    sort_obs_curr_down <- sort_obs_curr_down
  }
  sort_obs_curr_down <- sort(sort_obs_curr_down, decreasing = TRUE)
  observ_vect <- sort(sort_obs_curr_down, decreasing = TRUE)
  observ_vect <- cumsum(abs(observ_vect - departureFromValue))
  
  ### PREPARE THE RANDOM DATA HERE 
  shuff_curr_down <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/", curr_ratio_type, "_randomShuffleList.Rdata"))))
  sort_shuff_curr_down <- lapply(shuff_curr_down, function(x) {
    x <- sort(x, decreasing=TRUE)
    if(curr_ratio_type == "ratioDown" | curr_ratio_type == "rescWeightedQQ" | curr_ratio_type == "rescWeighted" ) {
      x <- abs(x - 0.5) + 0.5
    }
    x <- sort(x, decreasing=TRUE)
  })
  cumsum_sort_shuff_curr_down <- lapply(sort_shuff_curr_down, function(x) {
    cumsum(abs(x - departureFromValue))
  })

  ############################################################################################ PLOT THE 1ST PLOT
  # maxTADs <- max(c(unlist(lapply(all_random_List, function(x) lapply(x, length))), length(observ_vect) ))
  maxTADs <- max(c(unlist(lapply(shuff_curr_down, length))), length(observ_vect) )
  x_val <- 1:maxTADs
  y_range <- range(c(unlist(cumsum_sort_shuff_curr_down), observ_vect))
  stopifnot(!is.na(x_val))
  stopifnot(!is.na(y_range))
  
  # create list of equal length for the shuffled data
  maxShuffTADs <- max(c(unlist(lapply(shuff_curr_down, length))))
  sort_shuff_curr_down_equalsize <- lapply(sort_shuff_curr_down, function(x) x[1:maxShuffTADs])
  cumsum_sort_shuff_curr_down_equalsize <- lapply(cumsum_sort_shuff_curr_down, function(x) x[1:maxShuffTADs])
  stopifnot(0 == diff(unlist(lapply(cumsum_sort_shuff_curr_down_equalsize, length))))
  stopifnot(0 == diff(unlist(lapply(sort_shuff_curr_down_equalsize, length))))
  stopifnot(unique(unlist(lapply(sort_shuff_curr_down_equalsize, length))) == maxShuffTADs)
  stopifnot(unique(unlist(lapply(cumsum_sort_shuff_curr_down_equalsize, length))) == maxShuffTADs)
  sort_shuff_dt <- do.call(cbind, sort_shuff_curr_down_equalsize)
  rownames(sort_shuff_dt) <- NULL  # cannot compare TAD names across shufflings...
  stopifnot(nrow(sort_shuff_dt) == maxShuffTADs)
  cumsum_sort_shuff_dt <- apply(sort_shuff_dt,2 , function(x) cumsum(abs(x - departureFromValue)))
  rownames(cumsum_sort_shuff_dt) <- NULL  # cannot compare TAD names across shufflings...
  stopifnot(nrow(sort_shuff_dt) == nrow(cumsum_sort_shuff_dt))
  stopifnot(ncol(sort_shuff_dt) == ncol(cumsum_sort_shuff_dt))

  my_main <- paste0(curr_ratio_type, ": cumsum departure from ", departureFromValue)
  my_ylab <- paste0("cumsum(abs(", curr_ratio_type, " - ", departureFromValue,"))")
  my_xlab <- paste0("regions ranked by decreasing ", curr_ratio_type)
  
  if(plotSeparated){
    outFile <- paste0(curr_outFold, "/", curr_ratio_type, "_departure", departureFromValue,"_cumsum_obs_permut_withRatio_fullCurve.", plotType)
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  }
  # START BY PLOTTING EMPTY PLOT
  plot(NULL, 
       main = my_main,
       xlab=my_xlab,
       ylab=my_ylab,
       xlim = range(x_val),
       ylim = y_range,
       bty="l")
  
  # ADD THE TRUE DATA
  if(plotTrueType == "l") {
    lines(x = 1:length(observ_vect),
          y=observ_vect)
  } else if(plotTrueType == "p") {
    points(x = 1:length(observ_vect),
           y=observ_vect)
  } else {
    stop("error\n")
  }

  ## add the shuffle data:
  polygonPermutCol <- shufflePolygonCol
  stopifnot(!is.na(maxShuffTADs))
  stopifnot(nrow(cumsum_sort_shuff_dt) == maxShuffTADs)
  minLine <- apply(cumsum_sort_shuff_dt,1,min,na.rm=TRUE)
  maxLine <- apply(cumsum_sort_shuff_dt,1,max,na.rm=TRUE)
    
  stopifnot(length(minLine) == maxShuffTADs)
  stopifnot(length(maxLine) == maxShuffTADs)
  stopifnot(all(maxLine>=minLine))
  polygon(x = c(1:maxShuffTADs, rev(1:maxShuffTADs)), 
          y = c( minLine, rev(maxLine)),
          border=polygonPermutCol,
          col = polygonPermutCol)

  mtext(subtitDir, font=3)
  
  # add the AUC ratio
  # convert back the cumsumdt to list to input to the function # input is the ratio, not cumsum !
  auc_ratio <- calc_aucObsPerm_cumsum_list(observ_vect = sort_obs_curr_down,
                                                   shuff_List=sort_shuff_curr_down,
                                                   thresh_perm =permThresh,
                                                   doPlot=F,
                                                   departureValue=departureFromValue)
  legend("bottomright", legend= paste0("AUC ratio  = ", sprintf("%.2f", auc_ratio)), bty="n")

  all_ratio_split <- append(all_ratio_split, setNames(auc_ratio, paste0(curr_ratio_type, "_auc_shuffTADs_all")))
  
  if(plotSeparated){
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
  
  
  ############################################################################################ PLOT THE 2ND PLOT
  # => plot the first part of the curve where the observed data == 1
  if(plotSeparated){
    outFile <- paste0(curr_outFold, "/", curr_ratio_type, "_departure05_cumsum_obs_permut__withRatio_1stPart.", plotType)
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  }
  # take only the TADs that have a ratio == 1
  onlyOne_TADs <- names(sort_obs_curr_down)[sort_obs_curr_down == 1]
  onlyOne_sort_obs_curr_down <- sort_obs_curr_down[onlyOne_TADs]
  onlyOne_observ_vect <- observ_vect[names(observ_vect) %in% onlyOne_TADs]
  stopifnot(sort_obs_curr_down[1:length(onlyOne_observ_vect)] == 1)

  # draw the observed
  # START BY PLOTTING EMPTY PLOT
  plot(NULL, 
       main = my_main,
       xlab=my_xlab,
       ylab=my_ylab,
       xlim = c(0, length(onlyOne_observ_vect)),
       ylim = c(0, max(onlyOne_observ_vect)),
       bty="l")
  
  
  if(length(onlyOne_TADs) > 0) {
    # ADD THE TRUE DATA
    if(plotTrueType == "l") {
      lines(x = 1:length(onlyOne_observ_vect),
            y=onlyOne_observ_vect)
    } else if(plotTrueType == "p") {
      points(x = 1:length(onlyOne_observ_vect),
             y=onlyOne_observ_vect)
    } else {
      stop("error\n")
    }
    # ADD THE BLUE POLYGON FOR THE SHUFFTADS RANDOM DATA
    # draw the random
    # for the random data, take the first as many as 1 in the observed (ensure that not more TADs than in the shuffling)
    max_onlyOne <- min(length(onlyOne_observ_vect), nrow(sort_shuff_dt))
    onlyOne_sort_shuff_dt <- sort_shuff_dt[1:max_onlyOne,]
    onlyOne_cumsum_sort_shuff_dt <- apply(onlyOne_sort_shuff_dt,2 , function(x) cumsum(abs(x - departureFromValue)))
    stopifnot(nrow(onlyOne_cumsum_sort_shuff_dt) <= length(onlyOne_observ_vect))
    onlyOne_minLine <- apply(onlyOne_cumsum_sort_shuff_dt,1,min,na.rm=TRUE)
    onlyOne_maxLine <- apply(onlyOne_cumsum_sort_shuff_dt,1,max,na.rm=TRUE)
    stopifnot(all(onlyOne_maxLine>=onlyOne_minLine))
    polygon(x = c(1:length(onlyOne_maxLine), rev(1:length(onlyOne_maxLine))), 
            y = c( onlyOne_minLine, rev(onlyOne_maxLine)),
            border=polygonPermutCol,
            col = polygonPermutCol)
    
    mtext(subtitDir, font=3)
    # add the AUC ratio
    # convert back the cumsumdt to list to input to the function -> should pass the not cumsum to the function
    onlyOne_shuffList <- lapply(1:ncol(onlyOne_sort_shuff_dt), function(x) onlyOne_sort_shuff_dt[,x] )
    stopifnot(unlist(lapply(onlyOne_shuffList,  length)) == max_onlyOne)
    auc_ratio_onlyOne <- calc_aucObsPerm_cumsum_list(observ_vect = onlyOne_sort_obs_curr_down,
                                                     shuff_List=onlyOne_shuffList,
                                                thresh_perm =permThresh,
                                                doPlot=F,
                                                departureValue=departureFromValue)
    legend("bottomright", legend= paste0("AUC ratio  = ", sprintf("%.2f", auc_ratio_onlyOne)), bty="n")
    
    all_ratio_split <- append(all_ratio_split, setNames(auc_ratio_onlyOne, paste0(curr_ratio_type, "_auc_shuffTADs_one")))
    
  } else {
    legend("bottomright", legend= paste0("undef - no obs. TADs with ratio equal 1"), bty="n")
    
    all_ratio_split <- append(all_ratio_split, setNames(NA, paste0(curr_ratio_type, "_auc_shuffTADs_one")))
    
  }
  if(plotSeparated){
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
  
  ############################################################################################ PLOT THE 3D PLOT
  # => plot the second part of the curve where the observed data < 1
  if(plotSeparated){
    outFile <- paste0(curr_outFold, "/", curr_ratio_type, "_departure05_cumsum_obs_permut_withRatio_2ndPart.", plotType)
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  } 
  # take only the TADs with ratios < 1
  keepFrom <- min(which(sort_obs_curr_down < 1))
  stopifnot(keepFrom == length(onlyOne_observ_vect) + 1)
  noOne_observ_vect <- observ_vect[keepFrom:length(observ_vect)]
  stopifnot(sort_obs_curr_down[keepFrom:length(sort_obs_curr_down)] < 1)
  noOne_sort_obs_curr_down <- sort_obs_curr_down[keepFrom:length(sort_obs_curr_down)]
  
  # draw the observed
  # START BY PLOTTING EMPTY PLOT
  plot(NULL, 
       main = my_main,
       xlab=my_xlab,
       ylab=my_ylab,
       xlim = c(0, length(noOne_observ_vect)),
       ylim = c(0, max(noOne_observ_vect)),
       bty="l")
  # ADD THE TRUE DATA
  if(plotTrueType == "l") {
    lines(x = 1:length(noOne_observ_vect),
          y=noOne_observ_vect)
  } else if(plotTrueType == "p") {
    points(x = 1:length(noOne_observ_vect),
           y=noOne_observ_vect)
  } else {
    stop("error\n")
  }
  # ADD THE BLUE POLYGON FOR THE SHUFFTADS RANDOM DATA
  if(keepFrom <= nrow(cumsum_sort_shuff_dt)) {
    # draw the random
    # for the random data, take the first as many as 1 in the observed (ensure that not more TADs than in the shuffling)
    noOne_sort_shuff_dt <- sort_shuff_dt[keepFrom:nrow(sort_shuff_dt),]
    noOne_cumsum_sort_shuff_dt <- apply(noOne_sort_shuff_dt,2 , function(x) cumsum(abs(x - departureFromValue)))
    stopifnot(nrow(noOne_cumsum_sort_shuff_dt) <= length(noOne_observ_vect))
    noOne_minLine <- apply(noOne_cumsum_sort_shuff_dt,1,min,na.rm=TRUE)
    noOne_maxLine <- apply(noOne_cumsum_sort_shuff_dt,1,max,na.rm=TRUE)
    stopifnot(all(noOne_maxLine>=noOne_minLine))
    polygon(x = c(1:length(noOne_maxLine), rev(1:length(noOne_maxLine))), 
            y = c( noOne_minLine, rev(noOne_maxLine)),
            border=polygonPermutCol,
            col = polygonPermutCol)
    
    mtext(subtitDir, font=3)
    # add the AUC ratio
    # convert back the cumsumdt to list to input to the function # for the function to calculate, do not use the cumsum
    noOne_shuffList <- lapply(1:ncol(noOne_sort_shuff_dt), function(x) noOne_sort_shuff_dt[,x] )
    stopifnot(length(unique(unlist(lapply(noOne_shuffList,  length)))) == 1 )
    auc_ratio_noOne <- calc_aucObsPerm_cumsum_list(observ_vect = noOne_sort_obs_curr_down,
                                                     shuff_List=noOne_shuffList,
                                                     thresh_perm =permThresh,
                                                     doPlot=F,
                                                     departureValue=departureFromValue)
    legend("bottomright", legend= paste0("AUC ratio  = ", sprintf("%.2f", auc_ratio_noOne)), bty="n")
    
    all_ratio_split <- append(all_ratio_split, setNames(auc_ratio_noOne, paste0(curr_ratio_type, "_auc_shuffTADs_2ndPart")))
    
    
  } else {
    legend("bottomright", legend= paste0("undef. (no random data available)"), bty="n")
    
    all_ratio_split <- append(all_ratio_split, setNames(NA, paste0(curr_ratio_type, "_auc_shuffTADs_2ndPart")))
    
  }

  
  if(plotSeparated){
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
  
  
} # end-for iterating over ratios

if(!plotSeparated){
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

outFile <- file.path(curr_outFold, "all_ratio_split.Rdata")
save(all_ratio_split, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))


