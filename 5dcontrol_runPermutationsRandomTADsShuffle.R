#!/usr/bin/Rscript

startTime <- Sys.time()

set.seed(20180202)


################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: pipeline_regionList.Rdata
# - script0: rangeTADgenes.Rdata
# - script1: DE_topTable.Rdata
# - script1: DE_geneList.Rdata
################################################################################

################  OUTPUT
# - permutationsList_<sizeRatio>.Rdata
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
script_name <- "5d_runPermutationsRandomTADsShuffle"
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
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)

nRandom <- ifelse(SSHFS, 5, nRandomPermutShuffle) # retrieved from main_settings.R

### generate randomTADs using the shuffle randomization
shuffle_sourceFile <- paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/ezh2_utils_fct.R")

shuffle_chromoPartition_v1 <- local({
  source(shuffle_sourceFile, local = TRUE)
  environment(shuffle_chromoPartition_v1) <- .GlobalEnv
  shuffle_chromoPartition_v1
})
# shuffle_chromoPartition_v1 <- function(domainDT, chrSize, preservePattern=TRUE, seed = NULL)

#******************************************************************************* !! HARD CODED
# shuffling start positions of the TADs
set.seed(20180126)
# !!! WARNING !!!
# compared to the other way of permutations, here the shuffled chr1_TAD1 does not correspond to the true TAD1 !!! (different start, end, size)
# the positions of the random TADs is not recorded
# in the shuffled partition, there is no gap between TADs
# similar to 5b: TAD randomization, but here randomization by shuffling TADs
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
head(gene2tadDT)

### PREPARE THE GENES I WILL ASSIGN TO THE RANDOM TADs
# do not take the genes on chrY
gene2assignDT <- entrezDT[entrezDT$entrezID %in% gene2tadDT$entrezID,]
head(gene2assignDT)
stopifnot(nrow(gene2tadDT) == nrow(gene2assignDT))
# gene2assignDT <- gene2assignDT[,c("entrezID", "chromo", "start", "end", "strand")]
stopifnot(c("chromo", "start", "end", "entrezID") %in% colnames(gene2assignDT))
if(gene2tadAssignMethod == "startPos") {
  stopifnot("strand" %in% colnames(gene2assignDT))
}

# take only those for which I have expression data !!!
gene2assignDT <- gene2assignDT[as.character(gene2assignDT$entrezID) %in% DE_topTable$genes,]

# INPUT DATA -TADs
tadDT <- read.delim(TADpos_file, header=F, col.names = c("chromo", "region",  "start", "end"), stringsAsFactors = F)
maxEndDT <- aggregate(end~chromo, data=tadDT, FUN=max, na.rm=TRUE)
tadDT <- tadDT[grepl("_TAD", tadDT$region),]
tadDT$size <- tadDT$end - tadDT$start + 1
domainDT <- tadDT[,c("chromo", "start", "end")] 

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

# for the given number of permutations
# for each chromosome
all_chromo <- unique(as.character(domainDT$chromo))
#stopifnot(length(all_chromo)  >= 22)
stopifnot(all_chromo %in% paste0("chr", c(1:22, "X")))

################################****************************************************************************************
####################################################### RUN PERMUTATIONS FOR ALL SIZE RATIOS
################################****************************************************************************************
  
all_permDT <- foreach(perm = 1:nRandom) %dopar%{
  cat("... start random TADs:", perm, "/", nRandom, "\n")
  allChr_randTADdt <- foreach(chromo = all_chromo, .combine="rbind") %do% {
    chrEnd <- maxEndDT$end[maxEndDT$chromo == chromo]
    # select the initial set of TADs for this chromosome
    chromo_domainDT <- domainDT[domainDT$chromo == chromo,]
    stopifnot(nrow(chromo_domainDT) > 0)
    randTADdt <- shuffle_chromoPartition_v1(domainDT=chromo_domainDT, chrSize = chrEnd , preservePattern=FALSE)
    # cat("perm", perm, " - ", chromo, " - after shuffle\n")
    randTADdt$region <- paste0(chromo, "_TAD", 1:nrow(randTADdt))
    stopifnot(nrow(randTADdt) > 0)
    # the area of the chromo covered by TADs should be the same before/after shuffling !
    stopifnot(abs(sum(chromo_domainDT$end-chromo_domainDT$start) - sum(randTADdt$end-randTADdt$start)) < 1e-10)
    ## ADD THE GAPS IF NOT ONLY RUNNING PIPELINE FOR TADs
    
    # save(chromo_domainDT, file=paste0("chromo_domainDT_", perm, ".Rdata"))
    
    if(!useTADonly) {
      gapsDT <- data.frame(chromo = chromo, 
                 start = c(1, randTADdt$end+1),
                 end = c(randTADdt$start-1, chrEnd),
                 stringsAsFactors = FALSE)
      gapsDT <- gapsDT[gapsDT$end > gapsDT$start,]
      if(nrow(gapsDT) > 1) {
        stopifnot(gapsDT$start[1] == 1)
        gapsDT$region <- paste0(chromo, "_BOUND", 1:nrow(gapsDT))
        filled_randTADdt <- rbind(randTADdt, gapsDT)  
      } else {
        filled_randTADdt <- randTADdt  
      }
      filled_randTADdt <- filled_randTADdt[order(filled_randTADdt$start),]
      filled_randTADdt <- unique(filled_randTADdt)
      stopifnot(nrow(filled_randTADdt) >= nrow(randTADdt))  
      stopifnot(gapsDT$end[nrow(gapsDT)] == chrEnd)
      if(nrow(filled_randTADdt) > 1) stopifnot(filled_randTADdt$end[1:(nrow(filled_randTADdt)-1)] + 1 == filled_randTADdt$start[2:nrow(filled_randTADdt)])

      randTADdt <- filled_randTADdt
    }
    
    # save(randTADdt, file=paste0("randTADdt_", perm, ".Rdata"))
    
    # ensure not overlapping
    if(nrow(randTADdt) > 1) {
      for(i in 2:nrow(randTADdt))
        stopifnot(randTADdt$start[i] > randTADdt$end[i-1])
    }
    # ensure all starts smaller than ends
    stopifnot(randTADdt$start < randTADdt$end)
    # save(randTADdt, file=paste0("randTADdt_", perm, ".Rdata"))
    randTADdt
    
  } # end building new set of TADs for all chromosomes
  
  # now assign the genes to the shuffled TADs
  # GENES HAVE ALREADY BEEN FILTERED TO CONSIDER ONLY THOSE WITH MAPPING POSITION + EXPRESSION DATA AVAILABLE
  # ASSIGN THE GENES TO TADs
  # warning("!!! use this script if assigning genes to TADs based on start position !!!\n")
  if(useTADonly) {
    stopifnot(!any(grepl("_BOUND", allChr_randTADdt$region)))
  }
  # save(randTADdt, file=paste0("allChr_randTADdt_", perm, ".Rdata"))
  # save(gene2assignDT, file=paste0("gene2assignDT_", perm, ".Rdata"))
  all_gene_posDT <- assignGene2TADs(regionDT = allChr_randTADdt, geneDT = gene2assignDT, assignMethod = gene2tadAssignMethod)
  
  # save(all_gene_posDT, file=paste0("all_gene_posDT_", perm, ".Rdata"))
  
  stopifnot(!any(duplicated(all_gene_posDT$entrezID)))
  # this should be only genes for which I have expression data !
  stopifnot(all_gene_posDT$entrezID %in% DE_topTable$genes)
  
  ##### TEMPORARY SAVE FOR CHECK
  if(perm %in% c(1:10)){
    save(allChr_randTADdt, file = paste0(curr_outFold, "/", "allChr_randTADdt_shuffle_", perm,".Rdata"))
    save(all_gene_posDT, file = paste0(curr_outFold, "/", "all_gene_posDT_shuffle_", perm,".Rdata"))
  }
  ##############################
  if(useFilterSizeData){
    nGenesByTAD <- setNames(as.numeric(table(all_gene_posDT$region)), names(table(all_gene_posDT$region)))
    toKeep <- names(nGenesByTAD)[nGenesByTAD >= minNbrGeneTAD & nGenesByTAD <= rangeTADsize[2]]
    all_gene_posDT <- all_gene_posDT[as.character(all_gene_posDT$region) %in% toKeep,]
    stopifnot(nrow(all_gene_posDT) > 0)
  }
  
  assignedGenes <- setNames(all_gene_posDT$region, all_gene_posDT$entrezID)
  assignedGenes
}
names(all_permDT) <- paste0("region", 1:nRandom)

################################****************************************************************************************
####################################################### WRITE OUTPUT
################################****************************************************************************************

permutationsList <- all_permDT
stopifnot(length(permutationsList) == (nRandom))
save(permutationsList, file = paste0(curr_outFold, "/permutationsList_shuffle.Rdata"))
cat(paste0("... written: ", paste0(curr_outFold, "/permutationsList_shuffle.Rdata"), "\n"))

txt <- paste0(toupper(script_name), "> Number of permutations: ", nRandom, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> WARNING: genes that map to duplicated entrezID are removed ! \n")
printAndLog(txt, pipLogFile)

cat("!!! WARNING: with this shuffling, the name of the TADs do not correspond with the names of the initial \"true\" TADs !!! \n")


txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))

