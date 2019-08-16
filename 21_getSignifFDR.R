#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script11: emp_pval_combined.Rdata
################################################################################

################  OUTPUT
# - emp_pval_meanCorr.Rdata + plots
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

script_name <- "21_getSignifFDR"

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script6_name <- "6_runPermutationsMeanLogFC"
script7_name <- "7_runPermutationsMeanTADCorr"
script8_name <- "8c_runAllDown"

stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
#source("run_settings.R")
source(settingF)
# settingF = "SETTING_FILES_NOVOOM/run_settings_GSE71119_dediffSM_MFSM.R"
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

registerDoMC(ifelse(SSHFS, 2, nCpu))

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

# empFDR_thresh should be loaded from main_settings.R
stopifnot(exists("empFDR_tresh"))
stopifnot(empFDR_tresh <= 1 & empFDR_tresh >= 0)

tadDT <- read.delim(TADpos_file, header=F, col.names = c("chromo", "region",  "start", "end"), stringsAsFactors = F)
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))

################################****************************************************************************************
####################################################### FIND SIGNIF TADs INTERSECT FDR LOGFC INTRACORR THRESHOLDS
################################****************************************************************************************

###### 1) load data

############# load logFC data

cat("... load logFC data\n")

# obs_logFC_file <- file.path(curr_outFold, script3_name, "all_meanLogFC_TAD.Rdata")
obs_logFC_file <- file.path(pipOutFold, script3_name, "all_meanLogFC_TAD.Rdata")
stopifnot(file.exists(obs_logFC_file))
obs_vect_logFC <- eval(parse(text = load(obs_logFC_file)))

# shuff_logFC_file <- file.path(curr_outFold, script6_name, "meanLogFC_permDT.Rdata")
shuff_logFC_file <- file.path(pipOutFold, script6_name, "meanLogFC_permDT.Rdata")
stopifnot(file.exists(shuff_logFC_file))
permutDT_logFC <- eval(parse(text = load(shuff_logFC_file)))

# because this is absolute logFC
cut_off_seq_logFC <- seq(0,5,0.01)

############# load intraCorr data

cat("... load intraCorr data\n")

# obs_intraCorr_file <- file.path(curr_outFold, script4_name, "all_meanCorr_TAD.Rdata")
obs_intraCorr_file <- file.path(pipOutFold, script4_name, "all_meanCorr_TAD.Rdata")
stopifnot(file.exists(obs_intraCorr_file))
obs_vect_intraCorr <- eval(parse(text = load(obs_intraCorr_file)))

# shuff_intraCorr_file <- file.path(curr_outFold, script7_name, "meanCorr_permDT.Rdata")
shuff_intraCorr_file <- file.path(pipOutFold, script7_name, "meanCorr_permDT.Rdata")
stopifnot(file.exists(shuff_intraCorr_file))
permutDT_intraCorr <- eval(parse(text = load(shuff_intraCorr_file)))

cut_off_seq_intraCorr <- seq(0,1,0.01)

stopifnot(setequal(names(obs_vect_logFC), names(obs_vect_intraCorr)))


###### 2) compute FDR for a sequence of logFC and intraCorr thresholds

##### compute empFDR for logFC

cat("... Start computing empFDR for logFC\n")

empFDR_seq_logFC <- foreach(x = cut_off_seq_logFC , .combine='c') %dopar% {
  get_SAM_FDR(obs_vect_logFC, permutDT_logFC, cut_off = x, symDir = "symmetric", withPlot = F)
}
names(empFDR_seq_logFC) <- cut_off_seq_logFC

outFile <- file.path(curr_outFold, paste0("empFDR_seq_logFC.Rdata"))
save(empFDR_seq_logFC, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


##### compute empFDR for intraCorr

cat("... Start computing empFDR for intraCorr\n")

empFDR_seq_intraCorr <- foreach(x = cut_off_seq_intraCorr , .combine='c') %dopar% {
  get_SAM_FDR(obs_vect_intraCorr, permutDT_intraCorr, cut_off = x, symDir = "higher", withPlot = F)
}
names(empFDR_seq_intraCorr) <- cut_off_seq_intraCorr

outFile <- file.path(curr_outFold, paste0("empFDR_seq_intraCorr.Rdata"))
save(empFDR_seq_intraCorr, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

###### 3) transform the empFDR to make them <= 1 and monotonically decreasing; and select the smallest threshold that pass the desired empFDR threshold
#### transform the empFDR so that they are 1) not > 1 and 2) monotonically decreasing
# (code inspired from: csaw::empiricalFDR)
# forced monotonically decreasing example:
# if a higher FDR comes after, replaces current FDR with the higher FDR
# empFDR_seq_logFC_toy <- c(1, 0.97, 0.90, 0.92, 0.89, 0.75, 0.88, 0.86, 0.67, 0.66, 0.5, 0.6, 0.4, 0.2, 0.3)
# rev(cummax(rev(empFDR_seq_logFC_toy)))
# 1.00 0.97 0.90 0.92 0.89 0.75 0.88 0.86 0.67 0.66 0.50 0.60 0.40 0.20 0.30
# 1.00 0.97 0.92 0.92 0.89 0.88 0.88 0.86 0.67 0.66 0.60 0.60 0.40 0.30 0.30
# # (-> start from the right and retain the highest encounter value)
# empFDR_seq_logFC_toy <- c(1, 0.97, 0.4, 0.3, 0.7)
# rev(cummax(rev(empFDR_seq_logFC_toy)))
# 1.00 0.97 0.70 0.70 0.70

empFDR_seq_logFC <- empFDR_seq_logFC[!is.na(empFDR_seq_logFC) & ! is.infinite(empFDR_seq_logFC)]
empFDR_seq_logFC[empFDR_seq_logFC > 1] <- 1
empFDR_seq_logFC <- rev(cummax(rev(empFDR_seq_logFC)))
# which is the logFC threshold that leads to empFDR_thresh FDR ?
belowThresh_logFC <- empFDR_seq_logFC[empFDR_seq_logFC <= empFDR_tresh]

if(length(belowThresh_logFC) > 0) {
  # take the less stringent threshold that leads to this FDR
  logFC_thresh <- as.numeric(names(belowThresh_logFC[1]))
  stopifnot(!is.na(logFC_thresh))
  logFC_thresh_empFDR <- as.numeric(belowThresh_logFC[1])
  stopifnot(!is.na(logFC_thresh_empFDR))
  
  cat(paste0("> For emp. FDR = ", empFDR_tresh, " found: \n"))
  cat(paste0("... logFC threshold: ", logFC_thresh, " (empFDR = ", round(logFC_thresh_empFDR,2), ")\n"))
  
} else {
  logFC_thresh <- NA
  logFC_thresh_empFDR <- NA
  cat(paste0("> For emp. FDR = ", empFDR_tresh, " found: \n"))
  cat(paste0("... none logFC threshold for the desired FDR \n"))
}

empFDR_seq_intraCorr <- empFDR_seq_intraCorr[!is.na(empFDR_seq_intraCorr) & ! is.infinite(empFDR_seq_intraCorr)]
empFDR_seq_intraCorr[empFDR_seq_intraCorr > 1] <- 1
empFDR_seq_intraCorr <- rev(cummax(rev(empFDR_seq_intraCorr)))
belowThresh_intraCorr <- empFDR_seq_intraCorr[empFDR_seq_intraCorr <= empFDR_tresh]

if(length(belowThresh_intraCorr) > 0) {
  # take the less stringent threshold that leads to this FDR
  intraCorr_thresh <- as.numeric(names(belowThresh_intraCorr[1]))
  stopifnot(!is.na(intraCorr_thresh))
  intraCorr_thresh_empFDR <- as.numeric(belowThresh_intraCorr[1])
  stopifnot(!is.na(intraCorr_thresh_empFDR))
  cat(paste0("> For emp. FDR = ", empFDR_tresh, " found: \n"))
  cat(paste0("... intraCorr threshold: ", intraCorr_thresh, " (empFDR = ", round(intraCorr_thresh_empFDR,2), ")\n"))
  
} else{
  intraCorr_thresh <- NA
  intraCorr_thresh_empFDR <- NA
  cat(paste0("> For emp. FDR = ", empFDR_tresh, " found: \n"))
  cat(paste0("... none intraCorr threshold for the desired FDR \n"))
}

###### 4) retrieve signif TADs for logFC and intraCorr based on the thresholds found at step 3

stopifnot(setequal(names(obs_vect_intraCorr), names(obs_vect_logFC)))
stopifnot(length(obs_vect_intraCorr) == length(obs_vect_logFC))
nTotTADs <- length(obs_vect_logFC)

if(!is.na(logFC_thresh)) {
  signifTADs_logFC <- names(obs_vect_logFC[abs(obs_vect_logFC) >= logFC_thresh])
  cat("> Found # TADs with signif. logFC:\t", length(signifTADs_logFC), "/", nTotTADs, "\n")
} else {
  signifTADs_logFC <- character()
}

if(!is.na(intraCorr_thresh)) {
  signifTADs_intraCorr <- names(obs_vect_intraCorr[abs(obs_vect_intraCorr) >= intraCorr_thresh])
  cat("> Found # TADs with signif. intraCorr:\t", length(signifTADs_intraCorr), "/", nTotTADs, "\n")
} else{
  signifTADs_intraCorr <- character()
}

###### 5) final set of signif TADs = intersect signif. intraCorr signif. logFC
signifTADs <- intersect(signifTADs_intraCorr, signifTADs_logFC)

cat("> Found # TADs with signif. intraCorr + logFC:\t", length(signifTADs), "/", nTotTADs, "\n")

outFile <- file.path(curr_outFold, paste0("signifTADs.Rdata"))
save(signifTADs, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(curr_outFold, paste0("signifTADs_logFC.Rdata"))
save(signifTADs_logFC, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(curr_outFold, paste0("signifTADs_intraCorr.Rdata"))
save(signifTADs_intraCorr, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


###### 6) build output table
if(length(signifTADs) > 0) {
g2t_DT <- gene2tadDT[gene2tadDT$entrezID %in% geneList & gene2tadDT$region %in% signifTADs, c("entrezID", "region"),]
stopifnot(g2t_DT$entrezID %in% geneList)
g2t_DT$symbol <- unlist(sapply(g2t_DT$entrezID, function(x) names(geneList[geneList == x])))
g2t_DT_agg <- aggregate(.~region, function(x) paste0(x, collapse=","), data=g2t_DT)
tadDT <- tadDT[tadDT$region %in% g2t_DT_agg$region,]
out_DT <- left_join(tadDT, g2t_DT_agg, by="region")
out_DT$logFC <- unlist(sapply(out_DT$region, function(x) round(obs_vect_logFC[x], 4))) 
out_DT$intraCorr <- unlist(sapply(out_DT$region, function(x) round(obs_vect_intraCorr[x], 4))) 
out_DT$chromo <- factor(out_DT$chromo, levels = paste0("chr", c(1:22, "X")))
out_DT <- out_DT[order(abs(out_DT$logFC), out_DT$intraCorr, as.numeric(out_DT$chromo), out_DT$start),]
stopifnot(!is.na(out_DT))
} else{
  out_DT <- data.frame()
}
outFile <- file.path(curr_outFold, "signifTADs_empFDR_logFC_intraCorr_intersect.txt")
write.table(out_DT, file = outFile, sep="\t", quote=F, row.names=F, col.names=T)
cat(paste0("... written: ", outFile, "\n"))

all_empFDR_signifThresh <- c(
  desired_empFDR_thresh = empFDR_tresh,
  effective_logFC_empFDR_thresh = logFC_thresh_empFDR,
  logFC_thresh = logFC_thresh,
  effective_intraCorr_empFDR_thresh = intraCorr_thresh_empFDR,
  intraCorr_thresh = intraCorr_thresh
)
outFile <- file.path(curr_outFold, paste0("all_empFDR_signifThresh.Rdata"))
save(all_empFDR_signifThresh, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


all_empFDR_signifThresh <- round(all_empFDR_signifThresh, 4)

nbrSignif <- c(
  nSignifTADs_logFC = length(signifTADs_logFC),
  nSignifTADs_intraCorr = length(signifTADs_intraCorr),
  nSignifTADs_intersect = length(signifTADs)
)

out_thresh_DT <- t(data.frame(c(dataset=basename(pipOutFold), all_empFDR_signifThresh, nbrSignif)))

outFile <- file.path(curr_outFold, paste0("FDR_results_summary.txt"))
write.table(out_thresh_DT, file = outFile, col.names=T, row.names=F, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))


######################
################################################
txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))

