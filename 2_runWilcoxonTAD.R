#!/usr/bin/Rscript

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: pipeline_regionList.Rdata
# - script0: rna_geneList.Rdata
# - script0: pipeline_geneList.Rdata
# - script0: rna_madnorm_rnaseqDT.Rdata or rna_qqnorm_rnaseqDT.Rdata
################################################################################

################  OUTPUT
# - wilcox_ttest_meanTAD_qq.Rdata
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
script_name <- "2_runWilcoxonTAD"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
#source("run_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

# if microarray was not set in the settings file -> by default set to  FALSE
if(!exists("microarray")) microarray <- FALSE

suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

registerDoMC(ifelse(SSHFS, 2, nCpu))

# ADDED 16.11.2018 to check using other files
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)

#*******************************************************************************
if(microarray) {
  norm_rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_madnorm_rnaseqDT.Rdata"))))
}else {
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

# INPUT DATA
regionList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_regionList.Rdata"))))
if(useTADonly) {
  if(any(grepl("_BOUND", pipeline_regionList))) {
    stop("! data were not prepared for \"useTADonly\" !")
  }
}
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% as.character(geneList),]

########## LOAD THE DATA
samp1 <- eval(parse(text=load(paste0(setDir, "/", sample1_file))))
samp2 <- eval(parse(text=load(paste0(setDir, "/", sample2_file))))

stopifnot(all(samp1 %in% colnames(norm_rnaseqDT)))
stopifnot(all(samp2 %in% colnames(norm_rnaseqDT)))

# filter for the size (UPDATE: done in 0_prepGeneData)
regionsFilter <- regionList

txt <- paste0(toupper(script_name), "> Start with # regions: ", length(regionList), "/", length(unique(gene2tadDT$region)), "\n")
printAndLog(txt, pipLogFile)

if(useTADonly) {
    initLen <- length(regionsFilter)
    regionsFilter <- regionsFilter[grep("_TAD", regionsFilter)]
    txt <- paste0(toupper(script_name), "> Take only the TAD regions: ", length(regionsFilter),"/", initLen, "\n")
    printAndLog(txt, pipLogFile)
}

cat("... start Wilcoxon tests\n")

wilcox_ttest_meanTAD_qq <- foreach(i_reg = 1:length(regionsFilter)) %dopar% {
  cat(paste0("... wilcox tests, region ", i_reg, "/", length(regionsFilter), "\n"))
  reg <- regionsFilter[i_reg]
  reg_genes <- gene2tadDT$entrezID[gene2tadDT$region == reg]
  subData <- as.data.frame(norm_rnaseqDT[which(geneList[rownames(norm_rnaseqDT)] %in% reg_genes),,drop=F])
  subDataFoo <- subData
  subData <- as.data.frame(matrix(colMeans(subData), nrow=1))
  colnames(subData) <- colnames(subDataFoo)
  
  exprValues <- as.vector(as.matrix(subData))
  # unlist column by column (e.g. the first 5 values are the first 5 values of the 1st column)
  stopifnot(exprValues[1:5] == subData[1,1:5])
  stopifnot(exprValues[(length(exprValues)-4):length(exprValues)] == subData[1,(ncol(subData)-4):ncol(subData)])
  stopifnot(nrow(subData) == 1)
  sampleValues <- sapply(colnames(subData), function(x) rep(x, nrow(subData)))
  stopifnot(all(sampleValues[1:5] == colnames(subData)[1:5]))
  stopifnot(all(sampleValues[(length(sampleValues)-4):length(sampleValues)] == colnames(subData)[(ncol(subData)-4):ncol(subData)]))
  # 1 for samp1, 0 for samp2
  sampleBinary <- sapply(sampleValues, function(x) {
    ifelse(x %in% samp1, 1,
           ifelse(x %in% samp2, 0, NA))
  })
  stopifnot(all(!is.na(sampleBinary)))
  stopifnot(length(sampleBinary) == length(exprValues))
  samp1_values <- exprValues[sampleBinary == 1] 
  samp2_values <- exprValues[sampleBinary == 0]
  stopifnot( length(samp1_values) + length(samp2_values) == length(sampleBinary))
  ttest_stat <- t.test(samp1_values, samp2_values)$statistic
  names(ttest_stat) <- "ttest_stat"
  ttest_pval <- t.test(samp1_values, samp2_values)$p.value
  names(ttest_pval) <- "ttest_pval"
  wilcox_stat <- wilcox.test(samp1_values, samp2_values)$statistic
  names(wilcox_stat) <- "wilcox_stat"
  wilcox_pval <- wilcox.test(samp1_values, samp2_values)$p.value
  names(wilcox_pval) <- "wilcox_pval"
  c(ttest_stat, ttest_pval, wilcox_stat, wilcox_pval)
}

stopifnot(length(wilcox_ttest_meanTAD_qq) == length(regionsFilter))
names(wilcox_ttest_meanTAD_qq) <- regionsFilter

cat(paste0("... end Wilcoxon tests\n"))

save(wilcox_ttest_meanTAD_qq, file= paste0(curr_outFold, "/wilcox_ttest_meanTAD_qq.Rdata"))
cat(paste0("... written: ", paste0(curr_outFold, "/wilcox_ttest_meanTAD_qq.Rdata"), "\n"))


txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))


