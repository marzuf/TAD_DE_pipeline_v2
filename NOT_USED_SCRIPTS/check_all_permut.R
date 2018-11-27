startTime <- Sys.time()

### PROBLEM WITH
# "GSE71119_dediffSM_MFSM",
# GSE79362_prog_nonprog

comp_folders <- c(
   "EZH2",
   "EZH2_CL",
  "EZH2_CL_withDupGene2TAD",
  "EZH2_withDupGene2TAD",
  "GSE101521_control_mdd",
  # "GSE102073_stic_nostic",
  "GSE51799_control_carrier",
  "GSE52166_prePf_postPf",
  "GSE64810_control_carrier",
  "GSE64813_prePTSD_postPTSD",
  "GSE65540_before_after",
  "GSE66306_before_after",
  "GSE68719_norm_park",
   # "GSE71119_dediffSM_MFSM",
   # "GSE71119_undiffSM_LMSM",
  "GSE73765_noninf_list",
  "GSE73765_noninf_salm",
  "GSE73765_salm_list",
  "GSE74927_neg_pos",
  "GSE77314_normal_tumor",
  "GSE77509_normal_ptt",
  "GSE77509_normal_tumor",
  "GSE77509_ptt_tumor",
  "GSE79209_dysp_nodysp",
  "GSE79362_prog_nonprog",
  "GSE81046_noninf_list",
  "GSE81046_noninf_salm",
  "GSE81046_salm_list",
  "GSE81089_normal_nsclc",
  "GSE84231_lhb_lsc",
  "GSE84231_lhb_rhb",
  "GSE84231_lsc_rsc",
  "GSE84231_rhb_rsc",
  "GSE86356_quadMD1_quadNorm",
  "GSE86356_tibMD1_quadMD1",
  "GSE86356_tibMD1_tibNorm",
  "GSE86422_myoPre_myoPost",
  "GSE87194_control_schi",
  "GSE87340_ad_nl",
  "GSE92592_control_ipf",
  "GSE94631_adMono_adDC",
  "GSE94736_old_young",
  "TCGA_brca_lum_bas",
  "TCGA_crc_msi_mss",
  "TCGA_stad_msi_gs",
  "TCGA_ucec_msi_cnl"
)


# comp_folders <- c(
#   "TCGA_brca_lum_bas",
#   "TCGA_crc_msi_mss",
#   "TCGA_stad_msi_gs",
#   "TCGA_ucec_msi_cnl"
# )

# comp_folders <- c( "GSE71119_dediffSM_MFSM")

comp_folders <- file.path("OUTPUT_FOLDER", comp_folders)

# args <- commandArgs(trailingOnly = TRUE)
# stopifnot(length(args) == 1)
# settingF <- args[1]
# stopifnot(file.exists(settingF))

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script5_name <- "5_runPermutationsMedian"
script8_name <- "8c_runAllDown"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script10_name <- "10_runEmpPvalMeanTADCorr"
script11_name <- "11_runEmpPvalCombined"
script_name <- "check_all_permut"
stopifnot(file.exists(paste0(script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

# settingF <- "BUILDDT_SETTING_FILES/run_settings_buildDT.R"

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
# source("main_settings.R")
# #source("run_settings.R")
# source(settingF)
# source("TAD_DE_utils.R")
# source("calc_smile_score.R")
# source("calc_fit_ks.R")
# source("plot_smile_fct.R")
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
registerDoMC(ifelse(SSHFS,2,20))

# "permThresh" the quantile of permutations to take is loaded from main_settings.R
# retrieved from main_settings.R toPlotRanked
#*********************************************************

###############################################################################################################################
###############################################################################################################################   CHECK THE GENE LIST
###############################################################################################################################

for(curr_ds in comp_folders) {
  
  pipList <- eval(parse(text = load(paste0(curr_ds, "/0_prepGeneData/pipeline_geneList.Rdata"))))
  deList <- eval(parse(text = load(paste0(curr_ds, "/1_runGeneDE/DE_geneList.Rdata"))))
                         
  if(length(pipList) != length(deList)) {
    txt <- paste0("> For ", curr_ds, "found:\t", 
    "length(pipList)", " = ", length(pipList), "\t",
    "length(deList)", " = ", length(deList), "\n")
    cat(txt)
  }
  stopifnot(length(pipList) == length(deList))
}


###############################################################################################################################
###############################################################################################################################   CHECK THE PERMUTATIONS
###############################################################################################################################

toPlotRanked <- c("ratioDown", "FCdown", "meanRatioFC" , "prodConcord", "prodMeanConcord", "prodLogRatioNbr", "prodRatioSum")

all_scores_DT <- foreach(curr_ds = comp_folders, .combine = 'rbind') %dopar% {
  
  cat(paste0("**** START computing scores for: ", basename(curr_ds), "\n"))
  
  inFile <- paste0(curr_ds, "/", script5_name, "/permutationsDT.Rdata")
  cat(inFile, "\n")
  stopifnot(file.exists(inFile))
  permutationsDT <- eval(parse(text = load(inFile)))
  stopifnot(length(unique(apply(permutationsDT, 2,  function(x) length(unique(x))))) == 1)
    
  for(curr_ratio_type in toPlotRanked){
    file1 <- paste0(curr_ds, "/", script8_name, "/all_obs_", curr_ratio_type, ".Rdata")
    cat(paste0(file1, "\n"))
    if(!file.exists(file1)) {
      stop(paste0("ERROR: ",  file1 , " does not exist !\n"))
    }
    file2 <- paste0(curr_ds, "/", script8_name, "/", curr_ratio_type, "_permDT.Rdata")
    if(!file.exists(file2)) {
      stop(paste0("ERROR: ",  file2, " does not exist !\n"))
    }
    cat(paste0(file2, "\n"))
    
    obs_curr_down <- eval(parse(text = load(file1)))
    permut_currDown <- eval(parse(text = load(file2)))
    # ensure I used the same set of TADs for the permutation and for the calculation
    # (NB: would also be possible to filter the obs_curr_down, but not the permut_currDown)
    stopifnot(all(names(obs_curr_down) %in% rownames(permut_currDown)))
    stopifnot(all(rownames(permut_currDown) %in% names(obs_curr_down)))
  }
}


cat("***** OK - Checked - DONE ***** \n ")