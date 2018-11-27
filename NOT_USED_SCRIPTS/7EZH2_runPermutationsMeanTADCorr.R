SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

startTime <- Sys.time()

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script5_name <- "5_runPermutationsMedian"
script_name <- "7_runPermutationsMeanTADCorr"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
#source("run_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(ifelse(SSHFS, 2, nCpu)) # from main_settings.R

# if microarray was not set in the settings file -> by default set to  FALSE
if(!exists("microarray")) microarray <- FALSE

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

#****************************************************************************************************************************

withDiago <- FALSE
txt <- paste0(toupper(script_name), "> Take the diagonal when computing the mean of lower triangle: ", as.character(withDiago), "\n")
printAndLog(txt, pipLogFile)

#*******************************************************************************
# geneList -> entrezID, names original rownames
# UNCOMMENT ONE OF THE TWO !!!!
#if(useFilterData){
#    ### VERSION 1: script1_name => use only the genes used in DE
#     geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_geneList.Rdata"))))
#     qqnorm_rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_qqnorm_rnaseqDT.Rdata"))))
#     stopifnot(all(rownames(qqnorm_rnaseqDT) == names(geneList)))
#     txt <- paste0(toupper(script_name), "> Use filtered qqnorm_rnaseqDT from DE\n")
#     printAndLog(txt, pipLogFile)
#} else{
#    ### VERSION 0: script0_name => use all genes
#    geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_geneList.Rdata"))))
#    qqnorm_rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_qqnorm_rnaseqDT.Rdata"))))
#    stopifnot(all(rownames(qqnorm_rnaseqDT) == names(geneList)))
#    txt <- paste0(toupper(script_name), "> Use unfiltered qqnorm_rnaseqDT\n")
#    printAndLog(txt, pipLogFile)
#}


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


#*******************************************************************************

#****************************************************************************************************************************

# INPUT DATA
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% names(geneList),]

cat("... load permutation data ...\n")
permutationsDT <- eval(parse(text = load(paste0(pipOutFold, "/", script5_name, "/permutationsDT.Rdata"))))
all_regions <- sort(unique(as.character(permutationsDT[,2])))
# stopifnot(all_regions %in% rownames(permutationsDT) )


### ===>> change EZH2 data, gene can map to boundary and tads, so change from stop() to warning()
if(useTADonly) {
  if(! all (regexpr("_TAD",  gene2tadDT$region[gene2tadDT$entrezID %in% rownames(permutationsDT)]) > 0 )) {
    warning("make not sense to filter TAD genes after permutations if permutations were run with genes belonging to non-TAD regions\n")
  }
  initLen <- length(all_regions)
  all_regions <- all_regions[grep("_TAD", all_regions)]
  # if want to use only TADs, make more sens if permutations were run without the TADs
  txt <- paste0(toupper(script_name), "> Take only TAD regions: ", length(all_regions), "/", initLen, "\n")
  printAndLog(txt, pipLogFile)    
  if(length(all_regions) < initLen){
    warning("make not sense to filter TAD regions after permutations if permutations were run with genes belonging to non-TAD regions\n")
  }
}


#****************************************************************************************************************************
#****************************************************************************************************************************
# COMPUTE MEAN INTRA CORR BY TAD FOR THE PERMUTATIONS
#****************************************************************************************************************************
#****************************************************************************************************************************

cat("... start intraCorr permutDT \n")

intraTADcorr_permDT_allReg <- foreach(i_col = 1:ncol(permutationsDT), .combine='cbind') %dopar% {
  
  cat(paste0("... intraTAD correlation for permutation: ", i_col, "/", ncol(permutationsDT), "\n"))
  
  g2t_permDT <- data.frame(entrezID = rownames(permutationsDT), 
                           region = permutationsDT[,i_col], stringsAsFactors = F)
  g2t_permDT$entrezID <- as.character(g2t_permDT$entrezID)
  g2t_permDT$region <-  as.character(g2t_permDT$region)
  
  permutCorr <- sapply(unique(all_regions), function(reg) {
    reg_genes <- g2t_permDT$entrezID[g2t_permDT$region == reg]
    subData <- as.data.frame(t(norm_rnaseqDT[which(geneList %in% reg_genes),,drop=F]))
    
    # THE REGIONS HERE ARE UNFILTERED, SO IT IS POSSIBLE THAT THERE ARE SOME REGIONS WITH ONLY 1 GENE
    # (these regions have not been used for the logReg_meanTAD etc.)
    if(ncol(subData) < 2)
      return(NA)
    # columns => the genes
    # rows => the samples
    stopifnot(nrow(subData) == ncol(norm_rnaseqDT))
    ################################################# BECAUSE OF POSSIBLE DUPLICATED GENE ENTREZ ID
    # stopifnot(ncol(subData) == length(reg_genes))
    stopifnot(ncol(subData) == length(geneList[which(geneList %in% reg_genes)]))
    
    #### ALL CORRELATION
    corrMatrix_all <- cor(subData)
    # should be correlation of the genes
    ################################################# BECAUSE OF POSSIBLE DUPLICATED GENE ENTREZ ID
    # stopifnot(nrow(corrMatrix_all) == length(reg_genes))
    # stopifnot(ncol(corrMatrix_all) == length(reg_genes))
    stopifnot(ncol(corrMatrix_all) == length(geneList[which(geneList %in% reg_genes)]))
    stopifnot(nrow(corrMatrix_all) == length(geneList[which(geneList %in% reg_genes)]))
    
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

meanCorr_permDT <- as.data.frame(intraTADcorr_permDT_allReg)
stopifnot(ncol(meanCorr_permDT) == ncol(permutationsDT))  
colnames(meanCorr_permDT) <- paste0("permutation",  c(1:ncol(permutationsDT)))
stopifnot(nrow(meanCorr_permDT) == length(all_regions))
rownames(meanCorr_permDT) <- all_regions


txt <- paste0(toupper(script_name), "> Number of permutations for which mean logFC computed: ", ncol(meanCorr_permDT), "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> Number of regions for which mean logFC computed: ", nrow(meanCorr_permDT), "\n")
printAndLog(txt, pipLogFile)

save(meanCorr_permDT, file= paste0(curr_outFold, "/meanCorr_permDT.Rdata"))
cat(paste0("... written: ", curr_outFold, "/meanCorr_permDT.Rdata", "\n"))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))







