startTime <- Sys.time()

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
script_name <- "14i_cumulAllDown_limited_randomTADsShuffle"
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

stopifnot(any(c(step14_for_randomTADsShuffle  )))
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
  nColPlot <- 1  # here I plot only with the area
  # nRowPlot <- length(allDown_limited)*2/nColPlot # only 1 plot/row in this version
  nRowPlot <- length(allDown_limited)*1/nColPlot
  outFile <- paste0(curr_outFold, "/allRatios_cumsum_obs_permut.", plotType)
  do.call(plotType, list(outFile, height=myHeight*nRowPlot, width=myWidth*nColPlot))
  par(mfrow=c(nRowPlot, nColPlot))
}

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
    }
    observ_vect <- sort(sort_obs_curr_down, decreasing = TRUE)
    observ_vect <- cumsum(abs(observ_vect - departureFromValue))
    
    # load and store all the permut to have the limits for the plot limits
    # sort when adding them to the list !
    all_random_List <- list()
    
    if(step14_for_randomTADsFix){
      for(ratioSize in ratioSizeFix) {
        tmpDT <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/", curr_ratio_type, "_randomFixSizeList_", gsub("\\.", "",ratioSize), ".Rdata"))))
        tmpDT <- lapply(tmpDT, function(x) {
          x <- sort(x, decreasing=TRUE)
          if(curr_ratio_type == "ratioDown" | curr_ratio_type == "rescWeightedQQ" | curr_ratio_type == "rescWeighted" ) {
            x <- abs(x - 0.5) + 0.5
          }
          x <- sort(x, decreasing=TRUE)
          cumsum(abs(x - departureFromValue))
          })
        all_random_List[["fix"]][[paste0(ratioSize)]] <- tmpDT
      }  
    }
    
    if(step14_for_randomTADsGaussian){
      for(ratioSize in ratioSizeGaussian) {
        tmpDT <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/", curr_ratio_type, "_randomGaussianList_", gsub("\\.", "",ratioSize), ".Rdata"))))
        tmpDT <- lapply(tmpDT, function(x) {
          x <- sort(x, decreasing=TRUE)
          if(curr_ratio_type == "ratioDown" | curr_ratio_type == "rescWeightedQQ" | curr_ratio_type == "rescWeighted" ) {
            x <- abs(x - 0.5) + 0.5
          }
          x <- sort(x, decreasing=TRUE)
          cumsum(abs(x - departureFromValue))
        })
        all_random_List[["gauss"]][[paste0(ratioSize)]] <- tmpDT
      }  
    }
    
    if(step14_for_randomTADsShuffle){

      tmpDT <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/", curr_ratio_type, "_randomShuffleList.Rdata"))))
      tmpDT <- lapply(tmpDT, function(x) {
        x <- sort(x, decreasing=TRUE)
        if(curr_ratio_type == "ratioDown" | curr_ratio_type == "rescWeightedQQ" | curr_ratio_type == "rescWeighted" ) {
          x <- abs(x - 0.5) + 0.5
        }
        x <- sort(x, decreasing=TRUE)
        cumsum(abs(x - departureFromValue))
      })
      all_random_List[["shuffle"]] <- tmpDT
    }
    
    
    # range of the plots:
    maxTADs <- max(c(unlist(lapply(all_random_List, function(x) lapply(x, length))), length(observ_vect) ))
    x_val <- 1:maxTADs
    y_range <- range(c(unlist(all_random_List), observ_vect))
    stopifnot(!is.na(x_val))
    stopifnot(!is.na(y_range))
    
    my_main <- paste0(curr_ratio_type, ": cumsum departure from ", departureFromValue)
    my_ylab <- paste0("cumsum(abs(", curr_ratio_type, " - ", departureFromValue,"))")
    my_xlab <- paste0("regions ranked by decreasing ", curr_ratio_type)
    
    if(plotSeparated){
      outFile <- paste0(curr_outFold, "/", curr_ratio_type, "_departure", departureFromValue,"_cumsum_obs_permut_AUC.", plotType)
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

    # NOW FOR FIX SIZE AND GAUSSIAN, AND FOR EACH OF THE SIZE RATIOS, ADD THE CURVES/AREAS FOR THE RANDOM

    if(step14_for_randomTADsFix){
      subList <- all_random_List[["fix"]]
      foo <- lapply(1:length(subList), function(i) {
        sizeList <- subList[[i]]
        curr_ratioSize <- names(subList)[i]
        polygonPermutCol <- fixColors[curr_ratioSize]
        maxElements <- max(unlist(lapply(sizeList, length)))  
        stopifnot(!is.na(maxElements))
        sizeList2 <- lapply(sizeList, function(x) x[1:maxElements] <- x[1:maxElements])
        # minLine <- do.call(pmin, list(sizeList2, na.rm=TRUE)) # ??? why is this not working ???
        # maxLine <- do.call(pmax, list(sizeList2, na.rm=TRUE)) # ??? why is this not working ???
        sizeList2dt <- do.call(cbind,sizeList2)
        stopifnot(nrow(sizeList2dt) == maxElements)
        minLine <- apply(sizeList2dt,1,min,na.rm=TRUE)
        maxLine <- apply(sizeList2dt,1,max,na.rm=TRUE)
        stopifnot(length(minLine) == maxElements)
        stopifnot(length(maxLine) == maxElements)
        stopifnot(all(maxLine>=minLine))
        polygon(x = c(1:maxElements, rev(1:maxElements)), 
                y = c( minLine, rev(maxLine)),
                border=polygonPermutCol,
                col = polygonPermutCol)
      })
      
    }
    if(step14_for_randomTADsGaussian){
      subList <- all_random_List[["gauss"]]
      # for each of the sizeRatio
      # - extend the vectors so that they have all the same number of elements (max #)
      # - for each position find the min and max value
      foo <- lapply(1:length(subList), function(i) {
        sizeList <- subList[[i]]
        curr_ratioSize <- names(subList)[i]
        polygonPermutCol <- gaussianColors[[curr_ratioSize]]
        maxElements <- max(unlist(lapply(sizeList, length)))  
        stopifnot(!is.na(maxElements))
        sizeList2 <- lapply(sizeList, function(x) x[1:maxElements] <- x[1:maxElements])
        # minLine <- do.call(pmin, list(sizeList2, na.rm=TRUE)) # ??? why is this not working ???
        # maxLine <- do.call(pmax, list(sizeList2, na.rm=TRUE)) # ??? why is this not working ???
        sizeList2dt <- do.call(cbind,sizeList2)
        stopifnot(nrow(sizeList2dt) == maxElements)
        minLine <- apply(sizeList2dt,1,min,na.rm=TRUE)
        maxLine <- apply(sizeList2dt,1,max,na.rm=TRUE)
        stopifnot(length(minLine) == maxElements)
        stopifnot(length(maxLine) == maxElements)
        stopifnot(all(maxLine>=minLine))
        polygon(x = c(1:maxElements, rev(1:maxElements)), 
                y = c( minLine, rev(maxLine)),
                border=polygonPermutCol,
                col = polygonPermutCol)
      })
    }
    
    if(step14_for_randomTADsShuffle){
      sizeList <- all_random_List[["shuffle"]]
      polygonPermutCol <- shufflePolygonCol
      maxElements <- max(unlist(lapply(sizeList, length)))  
      stopifnot(!is.na(maxElements))
      sizeList2 <- lapply(sizeList, function(x) x[1:maxElements] <- x[1:maxElements])
      # minLine <- do.call(pmin, list(sizeList2, na.rm=TRUE)) # ??? why is this not working ???
      # maxLine <- do.call(pmax, list(sizeList2, na.rm=TRUE)) # ??? why is this not working ???
      sizeList2dt <- do.call(cbind,sizeList2)
      stopifnot(nrow(sizeList2dt) == maxElements)
      minLine <- apply(sizeList2dt,1,min,na.rm=TRUE)
      maxLine <- apply(sizeList2dt,1,max,na.rm=TRUE)
      stopifnot(length(minLine) == maxElements)
      stopifnot(length(maxLine) == maxElements)
      stopifnot(all(maxLine>=minLine))
      polygon(x = c(1:maxElements, rev(1:maxElements)), 
              y = c( minLine, rev(maxLine)),
              border=polygonPermutCol,
              col = polygonPermutCol)
    }

    if(step14_for_randomTADsFix | step14_for_randomTADsGaussian) {
      ratioLeg <- union(ratioSizeGaussian, ratioSizeFix)
      legend("topleft", lty=1, legend=ratioLeg, col = fixColors[as.character(ratioLeg)], bty="n")
    }
    
    mtext(subtitDir, font=3)
      if(plotSeparated){
        foo <- dev.off()
        cat(paste0("... written: ", outFile, "\n"))
      }
}

if(!plotSeparated){
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))


