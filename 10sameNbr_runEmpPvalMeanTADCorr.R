#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script4: all_meanCorr_TAD.Rdata
# - script7sameNbr: meanCorr_sample_around_TADs_sameNbr.Rdata
################################################################################

################  OUTPUT
# - emp_pval_meanCorr.Rdata
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
script4_name <- "4_runMeanTADCorr"
script7sameNbr_name <- "7sameNbr_runPermutationsMeanTADCorr"
script_name <- "10sameNbr_runEmpPvalMeanTADCorr"
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

# ADDED 16.11.2018 to check using other files
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)


corr_type <- "meanCorr"
txt <- paste0("taking sample correlation for corr_type\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)



### RETRIEVE ALL THE FILES IN THE FOLDER !!!
mainPipFold <- dirname(dirname(pipOutFold))
txt <- paste0("!!! take all the files matching \"meanCorr_sample_around_TADs_sameNbr.Rdata\" in ", mainPipFold, "\n")
printAndLog(txt, pipLogFile)

all_sampleCorr_files <- list.files(mainPipFold, pattern="meanCorr_sample_around_TADs_sameNbr.Rdata", full.names = TRUE, recursive = TRUE)
all_hicds <- list.files(mainPipFold)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(mainPipFold, x)))
stopifnot(length(all_sampleCorr_files) == length(unlist(all_exprds)))


### PREPARE THE SAMPLE CORR VALUES FROM ALL DATASETS
all_sample_corrValues <- foreach(corr_file = all_sampleCorr_files, .combine='c') %dopar% {
  stopifnot(file.exists(corr_file))
  corr_data <- load(corr_file)
  all_samp_corrs <-  unlist(lapply(corr_data, function(sub_data){
    lapply(sub_data, function(x) x[[paste0(corr_type)]])
  }))
  stopifnot(!is.null(all_samp_corrs))
  all_samp_corrs <- na.omit(all_samp_corrs)  
  all_samp_corrs
}


# RETRIEVE THE OBSERVED CORR DATA
obs_corr_file <- file.path(pipOutFold, script4_name, "all_meanCorr_TAD.Rdata")
stopifnot(file.exists(obs_corr_file))
all_obs_corr <- eval(parse(text = load(obs_corr_file)))
 
emp_pval_meanCorr <- sapply(all_obs_corr, function(x) {
  (sum(all_sample_corrValues >= x) + 1)/(length(all_sample_corrValues) + 1)
})
names(emp_pval_meanCorr) <- names(all_obs_corr)
stopifnot(all(emp_pval_meanCorr > 0 & emp_pval_meanCorr <= 1 ))


### CHECK RETRIEVE THE RIGHT VALUES
tadListFile <- file.path(pipOutFold, script0_name, "pipeline_regionList.Rdata")
stopifnot(file.exists(tadListFile))
pipeline_tadList <- eval(parse(text = load(tadListFile))) # not adjusted
stopifnot(setequal(names(emp_pval_meanCorr), pipeline_tadList))


outFile <- file.path(curr_outFold, "emp_pval_meanCorr.Rdata")
save(emp_pval_meanCorr, file= outFile)
cat(paste0("... written: ", outFile, "\n"))

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################


txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)
cat(paste0("*** DONE: ", script_name, "\n"))
       
          
        
        
