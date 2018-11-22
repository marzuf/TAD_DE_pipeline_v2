#!/usr/bin/Rscript

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script11: emp_pval_combined.Rdata
################################################################################

################  OUTPUT
# - emp_pval_meanCorr.Rdata + plots
################################################################################
SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

script_name <- "20_getSignifCombined"
script0_name <- "0_prepGeneData"
script11_name <- "11_runEmpPvalCombined"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
#source("run_settings.R")
source(settingF)
# settingF = "SETTING_FILES_NOVOOM/run_settings_GSE71119_dediffSM_MFSM.R"
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

stopifnot(exists("adjCombinedPval_thresh"))
stopifnot(adjCombinedPval_thresh <= 1 & adjCombinedPval_thresh >= 0)

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

tadDT <- read.delim(TADpos_file, header=F, col.names = c("chromo", "region",  "start", "end"), stringsAsFactors = F)
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))

################################****************************************************************************************
####################################################### FIND SIGNIF TADs BASED ON ADJ. COMBINED P-VAL
################################****************************************************************************************

# LOAD INPUT DATA
all_combined_pval <- eval(parse(text = load(paste0(pipOutFold, "/", script11_name, "/", "emp_pval_combined.Rdata"))))
stopifnot(names(all_combined_pval) %in% gene2tadDT$region)

# 1) correct the combined p-val
adj_all_combined_pval <- p.adjust(all_combined_pval, method="BH")

# 2) find the TADs that have adj. combined p-val smaller than the desired alpha
signifTADs <- names(adj_all_combined_pval[adj_all_combined_pval <= adjCombinedPval_thresh ])
cat(paste0("... found signif. TADs: ", length(signifTADs), "/", length(adj_all_combined_pval), "\n"))

# 3) build the output table
if(length(signifTADs) == 0) {
  out_DT <- data.frame()  
} else{
  g2t_DT <- gene2tadDT[gene2tadDT$entrezID %in% geneList & gene2tadDT$region %in% signifTADs, c("entrezID", "region"),]
  stopifnot(g2t_DT$entrezID %in% geneList)
  g2t_DT$symbol <- unlist(sapply(g2t_DT$entrezID, function(x) names(geneList[geneList == x])))
  g2t_DT_agg <- aggregate(.~region, function(x) paste0(x, collapse=","), data=g2t_DT)
  tadDT <- tadDT[tadDT$region %in% g2t_DT_agg$region,]
  out_DT <- left_join(tadDT, g2t_DT_agg, by="region")
  out_DT$adj_combinedPval <- unlist(sapply(out_DT$region, function(x) round(adj_all_combined_pval[x], 4))) 
  out_DT$chromo <- factor(out_DT$chromo, levels = paste0("chr", c(1:22, "X")))
  out_DT <- out_DT[order(out_DT$adj_combinedPval, as.numeric(out_DT$chromo), out_DT$start),]
  stopifnot(!is.na(out_DT))
}

outFile <- file.path(curr_outFold, "signifTADs_combinedPval.txt")
write.table(out_DT, file = outFile, sep="\t", quote=F, row.names=F, col.names=T)
cat(paste0("... written: ", outFile, "\n"))

if(nrow(out_DT) > 1) {
outFile <- file.path(curr_outFold, paste0("adjCombinedPval_density.", plotType))
plot(density(out_DT$adj_combinedPval), main="Adj. combnined Pval distribution", xlab="adj. combined Pval")
abline(v = adjCombinedPval_thresh, lty=2)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))
}

out_thresh_DT <- t(data.frame(c(dataset=basename(pipOutFold), adjPvalThresh = adjCombinedPval_thresh, nbrSignif = length(signifTADs))))
outFile <- file.path(curr_outFold, paste0("adjPvalCombined_results_summary.txt"))
write.table(out_thresh_DT, file = outFile, col.names=T, row.names=F, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(curr_outFold, "signifTADs.Rdata")
save(signifTADs, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(curr_outFold, "adjCombinedPval_thresh.Rdata")
save(adjCombinedPval_thresh, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

######################################################################
txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))

