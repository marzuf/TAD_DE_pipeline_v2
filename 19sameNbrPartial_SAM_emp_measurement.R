#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))


# for each permutation, compute how many solutions greater than k
# R = average number of random solutions with logFC > k
# => represents the number of false positives
# => FDR = R/N

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")


pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4partial_runMeanTADCorr"
#script6_name <- "610000_runPermutationsMeanLogFC" # => output only SAM values for the correlation using permutation sameNbr
script7sameNbr_name <- "7sameNbrPartial_runPermutationsMeanTADCorr"


script_name <- "19sameNbrPartial_SAM_emp_measurement"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

# cat(paste0("setDir = ", setDir, "\n"))
source("main_settings.R") # setDir is the main_settings not in run_settings
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
source(paste0(pipScriptDir, "/", "TAD_DE_utils_meanCorr.R"))

registerDoMC(nCpu)


# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))


plotType <- "svg"
myHeight <- ifelse(plotType == "png", 480, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

fixCutOffSeq <- TRUE


k_cut_off <- 0.5

cut_off_seq_intraCorr <- seq(0,1,0.01)

sort_plotQuantileCutOff <- TRUE


corr_type <- "meanCorr"
txt <- paste0("taking sample correlation for corr_type\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)


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
## obs_logFC_file <- file.path(curr_outFold, script3_name, "all_meanLogFC_TAD.Rdata")
#obs_logFC_file <- file.path(pipOutFold, script3_name, "all_meanLogFC_TAD.Rdata")
#stopifnot(file.exists(obs_logFC_file))
#obs_vect_logFC <- eval(parse(text = load(obs_logFC_file)))

## shuff_logFC_file <- file.path(curr_outFold, script6_name, "meanLogFC_permDT.Rdata")
#shuff_logFC_file <- file.path(pipOutFold, script6_name, "meanLogFC_permDT.Rdata")
#stopifnot(file.exists(shuff_logFC_file))
#permutDT_logFC <- eval(parse(text = load(shuff_logFC_file)))

## because this is absolute logFC
## cut_off_seq_logFC <- seq(0,5,0.5)
#cut_off_seq_logFC <- seq(0,5,0.05) # MZ: UPDATE 16.07.2019

############# intraCorr
# obs_intraCorr_file <- file.path(curr_outFold, script4_name, "all_meanCorr_TAD.Rdata")
obs_intraCorr_file <- file.path(pipOutFold, script4_name, "all_meanCorr_TAD.Rdata")
stopifnot(file.exists(obs_intraCorr_file))
obs_vect_intraCorr <- eval(parse(text = load(obs_intraCorr_file)))

# shuff_intraCorr_file <- file.path(curr_outFold, script7_name, "meanCorr_permDT.Rdata")
# shuff_intraCorr_file <- file.path(pipOutFold, script7sameNbr_name, "meanCorr_sample_around_TADs_sameNbr.Rdata")
# stopifnot(file.exists(shuff_intraCorr_file))
# permutDT_intraCorr <- eval(parse(text = load(shuff_intraCorr_file)))


### PREPARE THE SAMPLE CORR VALUES FROM ALL DATASETS (COPY FROM SCRIPT 10)

### RETRIEVE ALL THE FILES IN THE FOLDER !!!
mainPipFold <- dirname(dirname(pipOutFold))
txt <- paste0("!!! take all the files matching \"meanCorr_sample_around_TADs_sameNbr.Rdata\" in ", mainPipFold, "\n")
printAndLog(txt, pipLogFile)

all_sampleCorr_files <- list.files(mainPipFold, pattern="meanCorr_sample_around_TADs_sameNbr.Rdata", full.names = TRUE, recursive = TRUE)
all_sampleCorr_files <- all_sampleCorr_files[grepl(script7sameNbr_name, all_sampleCorr_files)]
all_hicds <- list.files(mainPipFold)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(mainPipFold, x)))
cat(paste0(all_sampleCorr_files, collapse="\n"))
cat(paste0(all_exprds, collapse="\n"))
cat("\n")
cat(length(all_sampleCorr_files))
cat("\n")
cat(length(unlist(all_exprds)))
cat("\n")
stopifnot(length(all_sampleCorr_files) == length(unlist(all_exprds)))

### PREPARE THE SAMPLE CORR VALUES FROM ALL DATASETS
corr_file="/media/electron/mnt/etemp/marie/Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/7sameNbr_runPermutationsMeanTADCorr/meanCorr_sample_around_TADs_sameNbr.Rdata"
all_sample_corrValues <- foreach(corr_file = all_sampleCorr_files, .combine='c') %dopar% {
  stopifnot(file.exists(corr_file))
  corr_data <- eval(parse(text = load(corr_file)))
  all_samp_corrs <- as.numeric(sapply(corr_data, function(x) x[[paste0(corr_type)]]))
  stopifnot(!is.null(all_samp_corrs))
  all_samp_corrs <- na.omit(all_samp_corrs)  
  all_samp_corrs
}
nbrPermut <- length(all_sampleCorr_files) # each dataset corresponds to 1 permut



# RETRIEVE THE OBSERVED CORR DATA
obs_corr_file <- file.path(pipOutFold, script4_name, "all_meanCorr_TAD.Rdata")
stopifnot(file.exists(obs_corr_file))




obs_corr_values <- eval(parse(text = load(obs_corr_file)))


    

    
    
        
  cat("...... ", "get_SAM_FDR_aroundTADs", "\n")
  # higher: sum(obs_vect >= cut_off)
  empFDR_seq <- sapply(cut_off_seq_intraCorr, function(x) 
    get_SAM_FDR_aroundTADs(obs_vect = obs_corr_values, 
                           permut_values = all_sample_corrValues, 
                           cut_off = x, symDir = "higher", 
                           nPermut = nbrPermut,
                           withPlot = F))
  names(empFDR_seq) <- as.character(cut_off_seq_intraCorr)
  
  
  curr_variable_plotName <- "meanTADCorr"
  
  # toKeep <- !is.infinite(empFDR_seq) & !is.na(empFDR_seq)
  # slopeFDR <- as.numeric(coef(lm(empFDR_seq[toKeep] ~ cut_off_seq_intraCorr[toKeep]))["cut_off_seq_intraCorr[toKeep]"])
  
  # ! higher
  nbrObservedSignif <- sapply(cut_off_seq_intraCorr, function(x) sum(obs_corr_values >= x))
  names(nbrObservedSignif) <- cut_off_seq_intraCorr
  
  # TRY VARIABLE CUT-OFFS: plot FDR and nbrObservSignif ~ cut_off
  cat("...... ", "plot_FDR_with_observedSignif", "\n")
  outFile <- file.path(curr_outFold, paste0("FDR_var_cut_off_", curr_variable_plotName, "_", "sameNbr", ".", plotType))
  dir.create(dirname(outFile), recursive = T)
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot_FDR_with_observedSignif(yaxis_empFDR_vect= empFDR_seq,
                               xaxis_cutoff= cut_off_seq_intraCorr, 
                               y2_obsSignif= nbrObservedSignif, 
                               variableName=curr_variable_plotName, 
                               feature_name="TADs")
  cat(paste0("... written: ", outFile, "\n"))
  foo <- dev.off()
  
  # PLOT ALL VALUES AND AREA PERMUT WITH FDR for the given cut-off
  cat("...... ", "get_SAM_FDR_aroundTADs", "\n")
  outFile <- file.path(curr_outFold, paste0("FDR_", k_cut_off, "_cut_off_", curr_variable_plotName, "_", "sameNbr", ".", plotType))
  dir.create(dirname(outFile), recursive = T)
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  get_SAM_FDR_aroundTADs(obs_vect = obs_corr_values, 
                         permut_values = all_sample_corrValues, 
                         cut_off = k_cut_off, 
                         variableName = curr_variable_plotName, 
                         symDir = "higher", 
                         withPlot = TRUE,
                         sortToPlot = sort_plotQuantileCutOff)
  # get_SAM_FDR(obs_vect, permutDT, cut_off = k_cut_off, variableName = curr_variable, symDir = "higher", withPlot = TRUE, minQuant=0, maxQuant = 1)
  textTAD <- names(obs_corr_values)[which(obs_corr_values >= k_cut_off  ) ] # higher !
  if(length(textTAD) > 0)
    text(y=obs_corr_values[textTAD], x = which(names(obs_corr_values) %in% textTAD), labels = textTAD, pos=2, col="gray")
  
  cat(paste0("... written: ", outFile, "\n"))
  foo <- dev.off()
  
meanCorr_empFDR <- list(
    empFDR = empFDR_seq,
    nbrSignif = nbrObservedSignif
  )
  



outFile <- file.path(curr_outFold, "meanCorr_empFDR.Rdata")
save(meanCorr_empFDR, file= outFile)
cat(paste0("... written: ", outFile, "\n"))

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################


txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)
cat(paste0("*** DONE: ", script_name, "\n"))




