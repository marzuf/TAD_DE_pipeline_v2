#!/usr/bin/Rscript

# modified 02.08.2019 -> to accomodate with the 100'000 without need of computing intraTADcorr and 8c script-> do the emp. FDR stuff only for the FC

options(scipen=100)

startTime <- Sys.time()

suppressPackageStartupMessages(library(data.table, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# empirical FDR
# logFC as effect size (measure of actual difference I want)

# rank by logFC
# set a cut-off k above which I say as significant
# N = number of real solutions with logFC > k

# for each permutation, compute how many solutions greater than k
# R = average number of random solutions with logFC > k
# => represents the number of false positives
# => FDR = R/N

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

source(paste0(setDir, "/mnt/ed4/marie/scripts/RNA_seq_v2_before0405/RNAseq_fct.R"))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script6_name <- "6_runPermutationsMeanLogFC"


script_name <- "19onlyFC_SAM_emp_measurement"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

# cat(paste0("setDir = ", setDir, "\n"))
source("main_settings.R") # setDir is the main_settings not in run_settings
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))


# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

# caller <- "DI"
# curr_dataset <- "TCGAbrca_lum_bas"
# outFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2", "SAM_emp_measurement", caller, curr_dataset)
# system(paste0("mkdir -p ", outFold))

# TADpos_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_", caller, "_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")    
# # file with assignment from entrez to all regions
# gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_", caller, "_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
# curr_outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER", curr_dataset)

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 480, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

fixCutOffSeq <- TRUE

# ADDED 15.08.2019 to check using other files
txt <- paste0("inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)


#**************************************************************************************************** COMBINED LOG FC AND INTRA TAD CORR

TADpos_DT <- read.delim(TADpos_file, header=F, stringsAsFactors = F, col.names=c("chromo", "region", "start", "end"))
gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))


############# logFC
# obs_logFC_file <- file.path(curr_outFold, script3_name, "all_meanLogFC_TAD.Rdata")
obs_logFC_file <- file.path(pipOutFold, script3_name, "all_meanLogFC_TAD.Rdata")
stopifnot(file.exists(obs_logFC_file))
obs_vect_logFC <- eval(parse(text = load(obs_logFC_file)))

# shuff_logFC_file <- file.path(curr_outFold, script6_name, "meanLogFC_permDT.Rdata")
shuff_logFC_file <- file.path(pipOutFold, script6_name, "meanLogFC_permDT.Rdata")
stopifnot(file.exists(shuff_logFC_file))
permutDT_logFC <- eval(parse(text = load(shuff_logFC_file)))

# because this is absolute logFC
cut_off_seq_logFC <- seq(0,5,0.5)
cut_off_seq_logFC <- seq(0,5,0.05) # MZ: UPDATE 16.07.2019



###########################################################################################


stopifnot(setequal(names(obs_vect_logFC), rownames(permutDT_logFC)))


commonReg <- sort(names(obs_vect_logFC))

# geneList_file <- file.path(curr_outFold, script0_name, "pipeline_geneList.Rdata")
geneList_file <- file.path(pipOutFold, script0_name, "pipeline_geneList.Rdata")
stopifnot(file.exists(geneList_file))
geneList <- eval(parse(text = load(geneList_file)))

g2t_DT <- gene2tad_DT[gene2tad_DT$entrezID %in% geneList,]
stopifnot(nrow(g2t_DT) > 0)
nbrGenes <- setNames(as.numeric(table(g2t_DT$region)), names(table(g2t_DT$region)))


nbrGenes <- nbrGenes[commonReg]
stopifnot(length(nbrGenes) > 0)

ii <- cut(nbrGenes, breaks = seq(min(nbrGenes), max(nbrGenes), len = 100),
          include.lowest = TRUE)
colors <- colorRampPalette(c("blue","white", "red"))(99)[ii]



empFDR_list <- list()

#**************************************************************************************************** EMP. FDR FROM LOG. FC
curr_variable <- "logFC"

obs_vect <- obs_vect_logFC
permutDT <- permutDT_logFC

curr_variable_plotName <- "abs(logFC)"

interRegion <- intersect(names(obs_vect), rownames(permutDT))
obs_vect <- obs_vect[interRegion]
permutDT <- permutDT[interRegion,]

if(fixCutOffSeq) {
  cut_off_seq <- cut_off_seq_logFC
} else{
  cut_off_seq <- round(seq(0, max(obs_vect[-which.max(obs_vect)]), length.out = 10),1) 
}
# symmetric: observ_N <- sum(abs(obs_vect) >= abs(cut_off))
empFDR_seq <- unlist(sapply(cut_off_seq, function(x) get_SAM_FDR(obs_vect, permutDT, cut_off = x, symDir = "symmetric", withPlot = F)))

toKeep <- !is.infinite(empFDR_seq) & !is.na(empFDR_seq)
slopeFDR <- as.numeric(coef(lm(empFDR_seq[toKeep] ~ cut_off_seq[toKeep]))["cut_off_seq[toKeep]"])

# symmetric ! - abs. logFC 
nbrObservedSignif <- unlist(sapply(cut_off_seq, function(x) sum(abs(obs_vect) >= abs(x))))

# TRY VARIABLE CUT-OFFS: plot FDR and nbrObservSignif ~ cut_off
outFile <- paste0(curr_outFold, "/", "FDR_var_cut_off_", curr_variable, ".", plotType)
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_FDR_with_observedSignif(yaxis_empFDR_vect= empFDR_seq,xaxis_cutoff= cut_off_seq, y2_obsSignif= nbrObservedSignif, variableName=curr_variable_plotName, feature_name="TADs")
cat(paste0("... written: ", outFile, "\n"))
foo <- dev.off()

# PLOT ALL VALUES AND AREA PERMUT WITH FDR for the given cut-off
k_cut_off <- 2
outFile <- paste0(curr_outFold, "/", "FDR_", k_cut_off, "_cut_off_", curr_variable, ".", plotType)
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
get_SAM_FDR(obs_vect, permutDT, cut_off = k_cut_off, variableName = curr_variable, symDir = "symmetric", withPlot = TRUE)
textTAD <- names(obs_vect)[which(obs_vect >= k_cut_off | obs_vect <= -k_cut_off ) ] # symmetric !
if(length(textTAD) > 0)
  text(y=obs_vect[textTAD], x = which(names(obs_vect) %in% textTAD), labels = textTAD, pos=2, col="gray")
# get_SAM_FDR(obs_vect, permutDT, cut_off = k_cut_off, variableName = curr_variable, symDir = "symmetric", withPlot = TRUE, minQuant=0, maxQuant = 1)
cat(paste0("... written: ", outFile, "\n"))
foo <- dev.off()

empFDR_list[[paste0("empFDR_", curr_variable)]] <- setNames(empFDR_seq, paste0(cut_off_seq))
empFDR_list[[paste0("nbrSignif_", curr_variable)]] <- setNames(nbrObservedSignif, paste0(cut_off_seq))
empFDR_list[[paste0("slopeEmpFDR_", curr_variable)]] <- slopeFDR

rm(obs_vect)
rm(permutDT)




##############################
outFile <- file.path(curr_outFold,"empFDR_list.Rdata")
save(empFDR_list, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

tmp <- unlist(empFDR_list[grep("empFDR", names(empFDR_list))])
tmp <- tmp[!is.na(tmp) & !is.infinite(tmp)]
if(any(tmp > 1)) {
  warning("!!! emp. FDR > 1 found !!!\n")
}

##############################
cat("***** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
