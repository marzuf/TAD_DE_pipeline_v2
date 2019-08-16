#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

### !!! NB 02.07.2019 THE cor() FUNCTION IS USED WITHOUT PASSING THE CORRELATION METHOD PARAMETER !!! -> DEFAULT "pearson"

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: pipeline_regionList.Rdata
# - script0: pipeline_geneList.Rdata
# - script0: rna_madnorm_rnaseqDT.Rdata or rna_qqnorm_rnaseqDT.Rdata
# - script5: permutationsDT.Rdata
################################################################################

################  OUTPUT
# - meanCorr_permDT.Rdata
################################################################################

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script5_name <- "5_runPermutationsMedian"
script_name <- "7_runPermutationsMeanTADCorr"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

source(paste0(pipScriptDir, "/", "my_save_pigz.R")) # to use customed fastSave save.pigz()

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(ifelse(SSHFS, 2, nCpu)) # from main_settings.R

# if microarray was not set in the settings file -> by default set to  FALSE
if(!exists("microarray")) microarray <- FALSE

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
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


########################***************************************** HARD-CODED !!!
withDiago <- FALSE
txt <- paste0(toupper(script_name), "> Take the diagonal when computing the mean of lower triangle: ", as.character(withDiago), "\n")
printAndLog(txt, pipLogFile)
#*******************************************************************************

################################****************************************************************************************
####################################################### PREPARE INPUT
################################****************************************************************************************
# INPUT DATA
# UPDATE SELECT THE GENES ACCORDING TO THE SETTINGS PREPARED IN 0_PREPGENEDATA
if(microarray) {
  norm_rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_madnorm_rnaseqDT.Rdata"))))
} else {
  norm_rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_qqnorm_rnaseqDT.Rdata")))) 
}
initList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_geneList.Rdata"))))
geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))

txt <- paste0(toupper(script_name), "> Start with # genes: ", length(geneList), "/", length(initList), "\n")
printAndLog(txt, pipLogFile)

norm_rnaseqDT <- norm_rnaseqDT[names(geneList),]    
stopifnot(all(rownames(norm_rnaseqDT) == names(geneList)))
stopifnot(!any(duplicated(names(geneList))))

gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% names(geneList),] 
# !!! 02.07.2019: THIS IS WRONG => gene2tadDT should be subset using geneList not names(geneList) !!!
# => but this is ok because gene2tadDT is not used later in the script

cat("... load permutation data ...\n")
permutationsDT <- eval(parse(text = load(paste0(pipOutFold, "/", script5_name, "/permutationsDT.Rdata"))))

if(ncol(permutationsDT) != nRandomPermut)
  stop("! NEED TO CHECK: different settings were used for running the permutations !\n")

all_regions <- sort(unique(as.character(permutationsDT[,2])))
# stopifnot(all_regions %in% rownames(permutationsDT) )

pipeline_regionList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_regionList.Rdata"))))
if(useTADonly) {
  if(any(grepl("_BOUND", pipeline_regionList))) {
    stop("! data were not prepared for \"useTADonly\" !")
  }
}
if(!setequal(pipeline_regionList, permutationsDT[,2])) {
  txtWarningRegion <- paste0(toupper(script_name), "> Not the same set of regions in permutDT and pipeline_regionList\n")
  printAndLog(txtWarningRegion, pipLogFile)    
  stop(txtWarningRegion)
} else {
  txtWarningRegion <- ""
}

if(!setequal(pipeline_geneList, rownames(permutationsDT))) {
  txtWarningGene <- paste0(toupper(script_name), "> Not the same set of genes in permutDT and pipeline_geneList\n")
  printAndLog(txtWarningGene, pipLogFile)    
  stop(txtWarningGene)
} else{
  txtWarningGene <- ""
}


if(useTADonly) {
  if(! all (regexpr("_TAD",  gene2tadDT$region[gene2tadDT$entrezID %in% rownames(permutationsDT)]) > 0 )) {
    stop("make not sense to filter TAD genes after permutations if permutations were run with genes belonging to non-TAD regions\n")
  }
  initLen <- length(all_regions)
  all_regions <- all_regions[grep("_TAD", all_regions)]
  # if want to use only TADs, make more sens if permutations were run without the TADs
  txt <- paste0(toupper(script_name), "> Take only TAD regions: ", length(all_regions), "/", initLen, "\n")
  printAndLog(txt, pipLogFile)    
  if(length(all_regions) < initLen){
    stop("make not sense to filter TAD regions after permutations if permutations were run with genes belonging to non-TAD regions\n")
  }
}

################################****************************************************************************************
####################################################### COMPUTE MEAN INTRA CORR BY TAD FOR THE PERMUTATIONS
################################****************************************************************************************

cat("... start intraCorr permutDT \n")

intraTADcorr_permDT_allReg <- foreach(i_col = 1:ncol(permutationsDT), .combine='cbind') %dopar% {
#intraTADcorr_permDT_allReg <- foreach(i_col = 1:2, .combine='cbind') %dopar% {
  
  cat(paste0("... intraTAD correlation for permutation: ", i_col, "/", ncol(permutationsDT), "\n"))
  
  g2t_permDT <- data.frame(entrezID = rownames(permutationsDT), 
                           region = permutationsDT[,i_col], stringsAsFactors = F)
  g2t_permDT$entrezID <- as.character(g2t_permDT$entrezID)
  g2t_permDT$region <-  as.character(g2t_permDT$region)
  
  permutCorr <- sapply(unique(all_regions), function(reg) {
    reg_genes <- g2t_permDT$entrezID[g2t_permDT$region == reg]
    stopifnot( reg_genes %in% geneList)
    subData <- as.data.frame(t(norm_rnaseqDT[which(geneList %in% reg_genes),,drop=F]))
    
    # THE REGIONS HERE ARE UNFILTERED, SO IT IS POSSIBLE THAT THERE ARE SOME REGIONS WITH ONLY 1 GENE
    # (these regions have not been used for the logReg_meanTAD etc.)
    if(ncol(subData) < 2)
      return(NA)
    # columns => the genes
    # rows => the samples
    stopifnot(nrow(subData) == ncol(norm_rnaseqDT))
    ################################################# BECAUSE OF POSSIBLE DUPLICATED GENE ENTREZ ID
    # UPDATE:  now should not have duplicates !!!
    stopifnot(ncol(subData) == length(reg_genes))
    # stopifnot(ncol(subData) == length(geneList[which(geneList %in% reg_genes)]))
    
    #### ALL CORRELATION
    corrMatrix_all <- cor(subData)
    # should be correlation of the genes
    ################################################# BECAUSE OF POSSIBLE DUPLICATED GENE ENTREZ ID
    # UPDATE:  now should not have duplicates !!!
    stopifnot(nrow(corrMatrix_all) == length(reg_genes))
    stopifnot(ncol(corrMatrix_all) == length(reg_genes))
    # stopifnot(ncol(corrMatrix_all) == length(geneList[which(geneList %in% reg_genes)]))
    # stopifnot(nrow(corrMatrix_all) == length(geneList[which(geneList %in% reg_genes)]))
    
    meanCorr_all <- mean(corrMatrix_all[lower.tri(corrMatrix_all, diag = withDiago)], na.rm=T)
    
  })
  # permutCorr <- rbindlist(permutCorr)
  curr_permutDT <- data.frame(permutCorr)
  stopifnot(ncol(curr_permutDT) == 1)
  colnames(curr_permutDT) <- paste0("result", i_col-1)
  stopifnot(all(rownames(curr_permutDT) == all_regions))
  curr_permutDT
}
cat("... end intraCorr permutDT \n")

cat(paste0("*** DONE: ", script_name, "\n"))
#stop("-- ok\n")

meanCorr_permDT <- as.data.frame(intraTADcorr_permDT_allReg)
stopifnot(ncol(meanCorr_permDT) == ncol(permutationsDT))  
colnames(meanCorr_permDT) <- paste0("permutation",  c(1:ncol(permutationsDT)))
stopifnot(nrow(meanCorr_permDT) == length(all_regions))
rownames(meanCorr_permDT) <- all_regions

################################****************************************************************************************
####################################################### WRITE OUTPUT
################################****************************************************************************************

txt <- paste0(toupper(script_name), "> Number of permutations for which mean logFC computed: ", ncol(meanCorr_permDT), "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> Number of regions for which mean logFC computed: ", nrow(meanCorr_permDT), "\n")
printAndLog(txt, pipLogFile)

#save(meanCorr_permDT, file= paste0(curr_outFold, "/meanCorr_permDT.Rdata"))
# update 16.08.2019 => faster save version
my_save.pigz(meanCorr_permDT, pigz_exec_path = pigz_exec_path, file= paste0(curr_outFold, "/meanCorr_permDT.Rdata"))

cat(paste0("... written: ", curr_outFold, "/meanCorr_permDT.Rdata", "\n"))

cat(paste0(txtWarningGene))
cat(paste0(txtWarningRegion))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))



cat("dim(meanCorr_permDT)\n")
cat(dim(meanCorr_permDT))
cat("\n")





