#!/usr/bin/Rscript

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: pipeline_geneList.Rdata
# - script1: DE_topTable.Rdata
# - script11: emp_pval_combined.Rdata
################################################################################

################  OUTPUT
# - [...].Rdata + plots
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
script11_name <- "11_runEmpPvalCombined"

script_name <- "23_funcEnrich_TopTableGenes_TADsGenes"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
#source("run_settings.R")
source(settingF)
# settingF = "SETTING_FILES_NOVOOM/run_settings_GSE71119_dediffSM_MFSM.R"
# settingF="../TAD_DE_pipeline/SETTING_FILES_cleanInput/run_settings_GSE101521_control_mdd.R"
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
# suppressPackageStartupMessages(library(clusterProfiler, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
# suppressPackageStartupMessages(library(org.Hs.eg.db, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 

# library(foreach)
# library(doMC)

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 480, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

### RETRIEVE FROM MAIN_SETTINGS.R
p_adj_DE_thresh <- 0.05
stopifnot(exists("p_adj_DE_thresh"))

nTopTADs <- 10
stopifnot(exists("nTopTADs"))

cat(paste0("... found setting: nTopTADs = ", nTopTADs, "\n"))
cat(paste0("... found setting: p_adj_DE_thresh = ", p_adj_DE_thresh, "\n"))

#****************************** FOR DEBUG
# nTopTADs <- 10
# p_adj_DE_thresh <- 0.05
#******************************

nShowPlot <- 10
################################****************************************************************************************
####################################################### PREPARE INPUT
################################****************************************************************************************

# INPUT DATA
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))

DE_topTable <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_topTable.Rdata"))))
stopifnot(DE_topTable$genes == rownames(DE_topTable))
stopifnot(names(geneList) %in% DE_topTable$genes)
curr_topTable_DT <- DE_topTable[DE_topTable$genes %in% names(geneList),]
stopifnot(curr_topTable_DT$genes %in% names(geneList)) # not true
stopifnot(nrow(curr_topTable_DT) > 0)
rm(DE_topTable)
curr_topTable_DT$entrezID <- unlist(sapply(curr_topTable_DT$genes, function(x) geneList[names(geneList) == x]))
stopifnot(!any(is.na(curr_topTable_DT$entrezID)))

curr_g2t_DT <- gene2tadDT[gene2tadDT$entrezID %in% geneList,]

pval_combined <- eval(parse(text = load(paste0(pipOutFold, "/", script11_name, "/", "emp_pval_combined.Rdata"))))
adj_pval_combined <- sort(p.adjust(pval_combined, method="BH"))

# retain the top nTopTADs TADs based on the combined pval
# use rank because of the ties !
topTADs <- names(rank(adj_pval_combined, ties="min")[rank(adj_pval_combined, ties="min") <= nTopTADs])
topTADs_genes <- as.character(c(unlist(sapply(topTADs, function(x) curr_g2t_DT$entrezID[curr_g2t_DT$region == x]))))
topTADs_genes_symb <- as.character(unlist(sapply(topTADs_genes, function(x) names(geneList[geneList == x]))))

cat(paste0("... # of genes derived from topTADs: ", length(topTADs_genes_symb), "\n"))

# select the DE genes
curr_topTable_DT <- curr_topTable_DT[order(curr_topTable_DT$adj.P.Val, curr_topTable_DT$P.Value),]
TopTable_DE_genes_symb <- curr_topTable_DT$genes[curr_topTable_DT$adj.P.Val <= p_adj_DE_thresh]
TopTable_DE_genes <- unlist(sapply(TopTable_DE_genes_symb, function(x) geneList[x]))

TopTable_notDE_genes_symb <- curr_topTable_DT$genes[curr_topTable_DT$adj.P.Val >= p_adj_DE_thresh]
TopTable_notDE_genes <- unlist(sapply(TopTable_notDE_genes_symb, function(x) geneList[x]))


### 
#                    DE gene | not DE gene
#  -----------------------------------
#  inTopTADs     |    #      |   #
# -------------------------------------
#  notInTopTADs  |    #      |   #
# 
#
# filled by column
DE_inTop <- sum(TopTable_DE_genes %in% topTADs_genes)
DE_notInTop <- sum(!TopTable_DE_genes %in% topTADs_genes)
stopifnot(sum(topTADs_genes %in% TopTable_DE_genes) == DE_inTop)
notDE_inTop <- sum(!topTADs_genes %in% TopTable_DE_genes)
stopifnot(sum(TopTable_notDE_genes %in% topTADs_genes) == notDE_inTop)
notDE_notInTop <- sum(! TopTable_notDE_genes %in% topTADs_genes)

testMatrix <- data.frame(DE_genes = c(DE_inTop, DE_notInTop),
                         notDE_genes = c(notDE_inTop, notDE_notInTop),
                         row.names=c("inTopTADs", "notInTopTADs"))
print(testMatrix)

FT_result <- fisher.test(testMatrix, alternative = "greater")

FT_result$p.value
FT_result$estimate

################################################################################################################
######################################################## HYPERGEO tests
################################################################################################################

# iterate over the topTADs and compute 1 test/topTAD

for(currTAD in topTADs) {
  
  currTADs_genes <- as.character(c(unlist(sapply(currTAD, function(x) curr_g2t_DT$entrezID[curr_g2t_DT$region == x]))))
  
  # filled by column
  DE_inCurrTAD <- sum(TopTable_DE_genes %in% currTADs_genes)
  DE_notInCurrTAD <- sum(!TopTable_DE_genes %in% currTADs_genes)
  stopifnot(sum(currTADs_genes %in% TopTable_DE_genes) == DE_inCurrTAD)
  notDE_inCurrTAD <- sum(!currTADs_genes %in% TopTable_DE_genes)
  stopifnot(sum(TopTable_notDE_genes %in% currTADs_genes) == notDE_inCurrTAD)
  notDE_notInCurrTAD <- sum(! TopTable_notDE_genes %in% currTADs_genes)
  
  testMatrix <- data.frame(DE_genes = c(DE_inCurrTAD, DE_notInCurrTAD),
                           notDE_genes = c(notDE_inCurrTAD, notDE_notInCurrTAD),
                           row.names=c("inTopTADs", "notInTopTADs"))
  print(testMatrix)
  
  FT_result <- fisher.test(testMatrix, alternative = "greater")
  
  FT_result$p.value
  FT_result$estimate
  
}

############################################################################################################################
############################################################################################################################
############################################################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))



