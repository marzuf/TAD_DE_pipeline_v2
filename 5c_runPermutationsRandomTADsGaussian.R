#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS

# - script0: pipeline_regionList.Rdata
# - script0: rangeTADgenes.Rdata
# - script1: DE_topTable.Rdata
# - script1: DE_geneList.Rdata
################################################################################

################  OUTPUT
# - permutationsList_gaussian_<sizeRatio>.Rdata
################################################################################

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))


args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script_name <- "5c_runPermutationsRandomTADsGaussian"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

# retrieved from main_settings.R or settingF
stopifnot(exists("gene2tadAssignMethod"))
if(gene2tadAssignMethod == "maxOverlap")
  suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

registerDoMC(ifelse(SSHFS, 2, nCpu))

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

# ADDED 27.11.2018 to check using other files
txt <- paste0("inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)

nRandom <- ifelse(SSHFS, 5, nRandomPermutGaussian) # retrieved from main_settings.R

### generate randomTADs by simply partitioning the genome in Gaussian distributed size regions
percentTADsize <- percentTADsizeGaussian # retrieve from main_settings.R

#******************************************************************************* !! HARD CODED
# shuffling start positions of the TADs
set.seed(20180126)
# !!! WARNING !!!
# compared to the other way of permutations, here the shuffled chr1_TAD1 does not correspond to the true TAD1 !!! (different start, end, size)
# the positions of the random TADs is not recorded
# in the shuffled partition, there is no gap between TADs

# similar to 5b but here the size of the random TADs is not fix to have equal size partition
# but drawn from a Gaussian distribution of the "true" initial TADs
#*******************************************************************************

################################****************************************************************************************
####################################################### PREPARE INPUT DATA
################################****************************************************************************************

rangeTADsize <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rangeTADgenes.Rdata"))))
stopifnot(rangeTADsize[1] >= minNbrGeneTAD)

DE_topTable <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_topTable.Rdata"))))
DE_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_geneList.Rdata"))))
stopifnot(nrow(DE_topTable) == length(DE_geneList))

entrezList <- unlist(sapply(DE_topTable$genes, function(x) DE_geneList[x]))
DE_topTable$genes <- entrezList
stopifnot(setequal(DE_topTable$genes, DE_geneList))

entrezDT <- read.delim(entrezDT_file, header=TRUE, stringsAsFactors = F)
entrezDT$entrezID <- as.character(entrezDT$entrezID)
stopifnot(all(DE_topTable$genes %in% entrezDT$entrezID))

gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
# don't do this here !
# gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% as.character(geneList),]

### PREPARE THE GENES I WILL ASSIGN TO THE RANDOM TADs
# do not take the genes on chrY
gene2assignDT <- entrezDT[entrezDT$entrezID %in% gene2tadDT$entrezID,]
stopifnot(nrow(gene2tadDT) == nrow(gene2assignDT))
gene2assignDT <- gene2assignDT[,c("entrezID", "chromo", "start", "end", "strand")]

# take only those for which I have expression data !!!
gene2assignDT <- gene2assignDT[as.character(gene2assignDT$entrezID) %in% DE_topTable$genes,]

# INPUT DATA -TADs
tadDT <- read.delim(TADpos_file, header=F, col.names = c("chromo", "region",  "start", "end"), stringsAsFactors = F)
maxEndDT <- aggregate(end~chromo, data=tadDT, FUN=max, na.rm=TRUE)
tadDT <- tadDT[grepl("_TAD", tadDT$region),]
tadDT$size <- tadDT$end - tadDT$start + 1
meanTADdt <- aggregate(size~chromo, data=tadDT, FUN=mean, na.rm=TRUE)
sdTADdt <- aggregate(size~chromo, data=tadDT, FUN=sd, na.rm=TRUE)
# for DI issue with chr19 1 value -> sd is NA
sdTADdt$size[is.na(sdTADdt$size)] <- 0

### take only the filtered data according to initial settings
pipeline_regionList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_regionList.Rdata"))))
if(useTADonly) {
  if(any(grepl("_BOUND", pipeline_regionList))) {
    stop("! data were not prepared for \"useTADonly\" !")
  }
}
initLen <- length(unique(gene2tadDT$region))
gene2tadDT <- gene2tadDT[gene2tadDT$region %in% pipeline_regionList,]
txt <- paste0(toupper(script_name), "> Take only filtered regions: ", length(unique(gene2tadDT$region)), "/", initLen, "\n")
printAndLog(txt, pipLogFile)

all_chromo <- unique(as.character(meanTADdt$chromo))
stopifnot(length(all_chromo)  >= 22)
stopifnot(all_chromo %in% paste0("chr", c(1:22, "X")))

################################****************************************************************************************
####################################################### RUN PERMUTATIONS FOR ALL SIZE RATIOS
################################****************************************************************************************

for(randSizeRatio in percentTADsize) {
  
  all_permDT <- foreach(perm = 1:nRandom) %dopar%{
    cat("... ratio", randSizeRatio, "- start random TADs:", perm, "/", nRandom, "\n")
    allChr_randTADdt <- foreach(chromo = all_chromo, .combine="rbind") %do% {
      chrEnd <- maxEndDT$end[maxEndDT$chromo == chromo]
      stopifnot(length(chrEnd) > 0)
      meanTADsize <- meanTADdt$size[meanTADdt$chromo == chromo]
      stopifnot(length(meanTADsize) > 0)
      sdTADsize <- sdTADdt$size[meanTADdt$chromo == chromo]
      stopifnot(length(sdTADsize) > 0)
      # > change here for the Gaussian TAD size
      allSizes <- numeric()
      while(TRUE) {
        newSize <- round(rnorm(n=1,mean=meanTADsize*randSizeRatio, sd=sdTADsize*randSizeRatio),0)
        if(newSize <= 0) next
        if(sum(allSizes)+newSize > chrEnd)
          break
        allSizes <- c(allSizes, newSize)
      }
      # add the remaining bp
      rmd_all <-  (chrEnd-sum(allSizes))%/%length(allSizes)
      stopifnot(rmd_all >= 0)
      if(rmd_all > 0){
        allSizes <- allSizes + rmd_all
      } 
      rmd_part <- chrEnd - sum(allSizes)
      stopifnot(rmd_part >=0)
      if(rmd_part > 0) {
        allSizes[1:rmd_part] <- allSizes[1:rmd_part] + 1  
      }
      stopifnot(sum(allSizes) == chrEnd)
      # shuffle the set of sizes -> not needed when drawn from Gaussian
      newEnd <- cumsum(allSizes)
      stopifnot(newEnd[length(newEnd)] == chrEnd)
      newStart <- c(1, newEnd[-length(newEnd)]+1)
      
      randTADdt <- data.frame(chromo=chromo,
                              start=newStart, 
                              end = newEnd,
                              region=paste0(chromo, "_TAD", 1:length(newStart)),
                              stringsAsFactors = FALSE)
      # ensure not overlapping
      if(nrow(randTADdt) > 1) {
        for(i in 2:nrow(randTADdt))
          stopifnot(randTADdt$start[i] > randTADdt$end[i-1])
      }
      # ensure all starts smaller than ends
      stopifnot(randTADdt$start < randTADdt$end)
      # stopifnot(abs(diff(range(randTADdt$start - randTADdt$end)))<=1) # not true when drawn from gaussian
      randTADdt
    } # end building new set of TADs for all chromosomes
    # now assign the genes 
    # no need to filter inter-TADs: when partitioning, I generate only TADs
    all_gene_posDT <- assignGene2TADs(regionDT = allChr_randTADdt, geneDT = gene2assignDT, assignMethod = gene2tadAssignMethod)
    
    # if needed, select the TADs containing appropriate number of genes
    if(useFilterSizeData){
      nGenesByTAD <- setNames(as.numeric(table(all_gene_posDT$region)), names(table(all_gene_posDT$region)))
      toKeep <- names(nGenesByTAD)[nGenesByTAD >= minNbrGeneTAD & nGenesByTAD <= rangeTADsize[2]]
      all_gene_posDT <- all_gene_posDT[as.character(all_gene_posDT$region) %in% toKeep,]
      stopifnot(nrow(all_gene_posDT) > 0)
    }  
    stopifnot(!any(duplicated(all_gene_posDT$entrezID)))
    assignedGenes <- setNames(all_gene_posDT$region, all_gene_posDT$entrezID)
    ##### TEMPORARY SAVE FOR CHECK
    # save(allChr_randTADdt, file = paste0(curr_outFold, "/", "allChr_randTADdt_", gsub("\\.","", randSizeRatio), "_", perm,".Rdata"))
    # save(all_gene_pos, file = paste0(curr_outFold, "/", "all_gene_pos_", gsub("\\.","", randSizeRatio), "_", perm,".Rdata"))
    ##############################
    assignedGenes
  }
  names(all_permDT) <- paste0("region", 1:nRandom)
  ####################################################### WRITE OUTPUT
  permutationsList <- all_permDT
  stopifnot(length(permutationsList) == (nRandom))
  save(permutationsList, file = paste0(curr_outFold, "/permutationsList_gaussian_", gsub("\\.","", randSizeRatio), ".Rdata"))
  cat(paste0("... written: ", paste0(curr_outFold, "/permutationsList_gaussian_", gsub("\\.","", randSizeRatio), ".Rdata"), "\n"))
} # end-for iterating over TAD size ratio

################################****************************************************************************************
####################################################### WRITE OUTPUT
################################****************************************************************************************

txt <- paste0(toupper(script_name), "> Number of permutations: ", nRandom, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> WARNING: genes that map to duplicated entrezID are removed ! \n")
printAndLog(txt, pipLogFile)

cat("!!! WARNING: with this shuffling, the name of the TADs do not correspond with the names of the initial \"true\" TADs !!! \n")

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))

