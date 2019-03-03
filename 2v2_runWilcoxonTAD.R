#!/usr/bin/Rscript

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: pipeline_regionList.Rdata
# - script0: rna_geneList.Rdata
# - script0: pipeline_geneList.Rdata
# - script0: rna_madnorm_rnaseqDT.Rdata or rna_qqnorm_rnaseqDT.Rdata
################################################################################

################  OUTPUT
# - wilcox_pairedTAD_meanExpr_fpkm.Rdata
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
script_name <- "2v2_runWilcoxonTAD"
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
  ### CHANGED HERE FOR THE V2 VERSION !!! 02.03.2019
  # norm_rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_madnorm_rnaseqDT.Rdata"))))
  norm_rnaseqDT <- eval(parse(text = load(file.path(pipOutFold, script0_name, "rna_fpkmDT.Rdata"))))
  
}else {
  ### CHANGED HERE FOR THE V2 VERSION !!! 02.03.2019
  # norm_rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_qqnorm_rnaseqDT.Rdata")))) 
  norm_rnaseqDT <- eval(parse(text = load(file.path(pipOutFold, script0_name, "rna_fpkmDT.Rdata"))))
  
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

wilcox_pairedTAD_meanExpr_fpkm <- foreach(i_reg = 1:length(regionsFilter)) %dopar% {

  cat(paste0("... wilcox tests, region ", i_reg, "/", length(regionsFilter), "\n"))
  reg <- regionsFilter[i_reg]
  reg_genes <- gene2tadDT$entrezID[gene2tadDT$region == reg]
  subData <- as.data.frame(norm_rnaseqDT[which(geneList[rownames(norm_rnaseqDT)] %in% reg_genes),,drop=F])
  
  cond1_DT <-  as.data.frame(norm_rnaseqDT[which(geneList[rownames(norm_rnaseqDT)] %in% reg_genes), samp1,drop=F])
  cond2_DT <- as.data.frame(norm_rnaseqDT[which(geneList[rownames(norm_rnaseqDT)] %in% reg_genes), samp2,drop=F])
  
  reg_genes_avgExpr_cond1 <- rowMeans(cond1_DT)
  reg_genes_avgExpr_cond2 <- rowMeans(cond2_DT)
  stopifnot(!is.na(reg_genes_avgExpr_cond1))
  stopifnot(!is.na(reg_genes_avgExpr_cond2))
  stopifnot(names(reg_genes_avgExpr_cond1) == names(reg_genes_avgExpr_cond2) )
  
  wTest <- wilcox.test(reg_genes_avgExpr_cond1, reg_genes_avgExpr_cond2, paired=TRUE, alternative="two.sided") 
  
  list(wilcoxTest_pval=as.numeric(wTest$p.value), 
       wilcoxTest_stat=as.numeric(wTest$statistic))
}

stopifnot(length(wilcox_pairedTAD_meanExpr_fpkm) == length(regionsFilter))
names(wilcox_pairedTAD_meanExpr_fpkm) <- regionsFilter

cat(paste0("... end Wilcoxon tests\n"))

outFile <- file.path(curr_outFold, "wilcox_pairedTAD_meanExpr_fpkm.Rdata")
save(wilcox_pairedTAD_meanExpr_fpkm, file= outFile)
cat(paste0("... written: ", outFile, "\n"))


txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))


