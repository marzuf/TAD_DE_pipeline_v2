#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script3: all_meanLogFC_TAD.Rdata
# - script6sameNbr: meanLogFC_sample_around_TADs_sameNbr.Rdata
################################################################################

# 8.4.20: rescaled version -> divide for each dataset the FC by the max value

################  OUTPUT
# - emp_pval_meanLogFC.Rdata
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
script3_name <- "3_runMeanTADLogFC"
script6sameNbr_name <- "6sameNbr_runPermutationsMeanLogFC"
script_name <- "9sameNbrRescaled_runEmpPvalMeanTADLogFC"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
#source("run_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

registerDoMC(ifelse(SSHFS,2, nCpu)) # loaded from main_settings.R

# if microarray was not set in the settings file -> by default set to  FALSE
if(!exists("microarray")) microarray <- FALSE


# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

fc_type <- "meanLogFC" #  "meanLogFC", "meanLogFC_right", "meanLogFC_left"

# ADDED 16.11.2018 to check using other files
txt <- paste0("inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", fc_type, "\n")
printAndLog(txt, pipLogFile)




### RETRIEVE ALL THE FILES IN THE FOLDER !!!
mainPipFold <- dirname(dirname(pipOutFold))
txt <- paste0("!!! take all the files matching \"meanLogFC_sample_around_TADs_sameNbr.Rdata\" in ", mainPipFold, "\n")
printAndLog(txt, pipLogFile)


#/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/6sameNbr_runPermutationsMeanLogFC/meanFC_sample_around_TADs_sameNbr.Rdata
all_sampleFC_files <- list.files(mainPipFold, pattern="meanFC_sample_around_TADs_sameNbr.Rdata", full.names = TRUE, recursive = TRUE)
stopifnot(length(all_sampleFC_files ) > 0)
all_sampleFC_files <- all_sampleFC_files[grepl(script6sameNbr_name, all_sampleFC_files)]  ### added 26.11.2019 otherwise match also the files for partial corr.
all_sampleFC_files <- all_sampleFC_files[!grepl("RANDOM", all_sampleFC_files) & !grepl("PERMUT", all_sampleFC_files)]

all_hicds <- list.files(mainPipFold)
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
stopifnot(length(all_hicds) == 30)

all_exprds <- sapply(all_hicds, function(x) list.files(file.path(mainPipFold, x)))


txt <- paste0("sampleFC_files used (n=", length(all_sampleFC_files), "):\n")
printAndLog(txt, pipLogFile)
txt <- paste0(all_sampleFC_files, collapse="\n")
printAndLog(txt, pipLogFile)



stopifnot(length(all_sampleFC_files) == length(unlist(all_exprds)))


### PREPARE THE SAMPLE CORR VALUES FROM ALL DATASETS
fc_file="/media/electron/mnt/etemp/marie/Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/7sameNbr_runPermutationsMeanTADLogFC/meanLogFC_sample_around_TADs_sameNbr.Rdata"
all_sample_fcValues <- foreach(fc_file = all_sampleFC_files, .combine='c') %dopar% {
  stopifnot(file.exists(fc_file))
  fc_data <- eval(parse(text = load(fc_file)))
  all_samp_fcs <- as.numeric(sapply(fc_data, function(x) x[[paste0(fc_type)]]))
  stopifnot(!is.null(all_samp_fcs))
  all_samp_fcs <- na.omit(all_samp_fcs)  
# added 8.4.20
  all_samp_fcs <- all_samp_fcs/max(abs(all_samp_fcs))
  stopifnot(!is.null(all_samp_fcs))
  stopifnot(!is.na(all_samp_fcs))
  all_samp_fcs
}


# RETRIEVE THE OBSERVED CORR DATA
obs_fc_file <- file.path(pipOutFold, script3_name, "all_meanLogFC_TAD.Rdata")
stopifnot(file.exists(obs_fc_file))
all_obs_fc <- eval(parse(text = load(obs_fc_file)))
 # added 8.4.20
all_obs_fc <- all_obs_fc/max(abs(all_obs_fc))

emp_pval_meanLogFC <- sapply(all_obs_fc, function(x) {



  if(x > 0){
    emp_pval <- sum(all_sample_fcValues >= x)   
  } else if (x < 0) {
    emp_pval <- sum(all_sample_fcValues <= x)   
  } else {
    emp_pval <- mean(c(sum(all_sample_fcValues >= x), sum(all_sample_fcValues <= x)), na.rm=T)
    warning("obs log FC == 0, will take the mean \n")
  }
  # emp_pval/length(shuff_logFC_TAD_region)
  ### ADD THE +1 => THIS IS FOR THE OBSERVED VALUE -> AVOID THE 0 VALUES !!!
  (emp_pval+1)/(length(all_sample_fcValues)+1)

})
names(emp_pval_meanLogFC) <- names(all_obs_fc)
stopifnot(all(emp_pval_meanLogFC > 0 & emp_pval_meanLogFC <= 1 ))


### CHECK RETRIEVE THE RIGHT VALUES
tadListFile <- file.path(pipOutFold, script0_name, "pipeline_regionList.Rdata")
stopifnot(file.exists(tadListFile))
pipeline_tadList <- eval(parse(text = load(tadListFile))) # not adjusted
stopifnot(setequal(names(emp_pval_meanLogFC), pipeline_tadList))


outFile <- file.path(curr_outFold, "emp_pval_meanLogFC.Rdata")
save(emp_pval_meanLogFC, file= outFile)
cat(paste0("... written: ", outFile, "\n"))

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################


txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)
cat(paste0("*** DONE: ", script_name, "\n"))
       
          
        
        
