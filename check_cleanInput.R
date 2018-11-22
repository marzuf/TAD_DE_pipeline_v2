SSHFS <- F

startTime <- Sys.time()

require(edgeR)

args <- commandArgs(trailingOnly = T)
caller <- args[1]

if(is.na(caller)) caller <- "DI"

setDir <- ifelse(SSHFS, "/media/electron", "")
pipFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller))
stopifnot(file.exists(pipFold))

settingFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput")
stopifnot(file.exists(settingFold))

outFold <- "CHECK_cleanInput"
system(paste0("mkdir -p ", outFold))

logFile <- paste0(outFold, "/", "check_cleanInput_logFile_", caller, ".txt")
system(paste0("rm -f ", logFile))

printAndLog <- function(txt, myfile){
  cat(txt)
  if(!is.null(myfile)) cat(txt, file = myfile, append=TRUE)
}

minCpmRatio <- 20/888

########################################################################################
######################################################################################## CHECK BY ITERATING OVER SETTING FILES
########################################################################################

all_setting_files <- list.files(settingFold, full.names = T, pattern = ".R$")

# run_settings_GSE84231_rhb_rsc.R => raw
# run_settings_GSE58135_ERpos_adjERpos.R => FPKM
# run_settings_GSE71119_dediffSM_MFSM.R => DESeq2
# run_settings_GSE81046_noninf_list.R => RSEM
settingFile <- file.path(settingFold, "run_settings_GSE84231_rhb_rsc.R")
settingFile <- file.path(settingFold, "run_settings_GSE58135_ERpos_adjERpos.R")
settingFile <- file.path(settingFold, "run_settings_GSE71119_dediffSM_MFSM.R")
settingFile <- file.path(settingFold, "run_settings_GSE81046_noninf_list.R")
settingFile <- file.path(settingFold, "run_settings_GSE40419_normal_cancer.R")
settingFile <- file.path(settingFold, "run_settings_GSE52166_prePf_postPf.R")

settingFile <- file.path(settingFold, "run_settings_GSE81046_noninf_list.R")


# all_setting_files <- file.path(settingFold, "run_settings_GSE52166_prePf_postPf.R")


for(settingFile in all_setting_files) {
  
  source(settingFile)
  stopifnot(exists("inputDataType"))
  
  txt <- paste0("\n*** ", caller, " - START ", gsub("run_settings_", "", basename(settingFile)),"\tinput: ", inputDataType, "\n")
  printAndLog(txt, myfile = logFile)
  
  # if not finished, skip
  lastFile <- file.path(pipFold, pipOutFold, "14f2_cumulAllDown_limited_AUC", "allRatios_cumsum_obs_permut.svg")
  if(!file.exists(lastFile)) {
    txt <- paste0("...... - not finished yet, skip \n")
    printAndLog(txt, myfile = logFile)
    next  
  }
  
  ######################################################################################## CHECK1: "raw" should be integer [raw only]
  if(inputDataType == "raw") {
    # stopifnot(inRdata)
    
    
    # load the expression data
    if(inRdata) {
      exprDT <- eval(parse(text = load(file.path(setDir, rnaseqDT_file))))
    } else {
      if(geneID_loc == "rn"){
        exprDT <- read.delim(paste0(setDir, rnaseqDT_file), sep=my_sep, header=T, row.names = 1)
      } else {
        exprDT <- read.delim(paste0(setDir, rnaseqDT_file), sep=my_sep, header=T)
        rownames(exprDT) <- exprDT[, geneID_loc]
        exprDT <- exprDT[,-geneID_loc]
      }
    }
    
    exprDT_init <- exprDT
    
    # is.integer(exprDT[1,1])
    #stopifnot(all(exprDT%%1 == 0))
    if(! all(na.omit(exprDT)%%1 == 0)) {
      txt <- paste0("...... !!! WARNING: \"raw\" inputDataType but non integers found !!!\n")
      printAndLog(txt, myfile = logFile)
    }
    
    samp1_id <- eval(parse(text = load(file.path(setDir, sample1_file))))
    samp2_id <- eval(parse(text = load(file.path(setDir, sample2_file))))
    
  ######################################################################################## CHECK2: fpkmDT should be exprDT [FPKM and RSEM]
  } else if(inputDataType == "RSEM" |inputDataType == "FPKM") {
    
    
    # load the expression data
    if(inRdata) {
      exprDT <- eval(parse(text = load(file.path(setDir, rnaseqDT_file))))
    } else {
      if(geneID_loc == "rn"){
        exprDT <- read.delim(paste0(setDir, rnaseqDT_file), sep=my_sep, header=T, row.names = 1)
      } else {
        exprDT <- read.delim(paste0(setDir, rnaseqDT_file), sep=my_sep, header=T)
        rownames(exprDT) <- exprDT[, geneID_loc]
        exprDT <- exprDT[,-geneID_loc]
      }
    }
    exprDT_init <- exprDT
    
    fpkmDT <- eval(parse(text = load(file.path(pipFold, pipOutFold, "0_prepGeneData", "rna_fpkmDT.Rdata"))))
    
    samp1_id <- eval(parse(text = load(file.path(setDir, sample1_file))))
    samp2_id <- eval(parse(text = load(file.path(setDir, sample2_file))))
    
    stopifnot(setequal(colnames(fpkmDT), c(samp1_id,samp2_id)))
    stopifnot(c(samp1_id,samp2_id) %in% colnames(exprDT))
    exprDT <- exprDT[,c(samp1_id, samp2_id)]
    fpkmDT <- fpkmDT[,c(samp1_id, samp2_id)]
    
    stopifnot(all(rownames(fpkmDT) %in% rownames(exprDT)))
    exprDT <- exprDT[rownames(fpkmDT), ]
    
    exprDT[is.na(exprDT)] <- 0
    fpkmDT[is.na(fpkmDT)] <- 0
    # is_na_expr <- is.na(exprDT)
    # exprDT <- exprDT[!is_na_expr]
    # fpkmDT <- fpkmDT[!is_na_expr]
    
    if(as.character(all.equal(exprDT, fpkmDT, check.attributes=FALSE)) != "TRUE") {
      txt <- paste0("...... !!! WARNING: RSEM or FPKM but exprDT != fpkmDT \n")
      # cat(txt)
      printAndLog(txt, myfile = logFile)
    }
  }
  ######################################################################################## CHECK3: any colsums with 0 [all datasets]
  exprDT <- eval(parse(text = load(file.path(pipFold, pipOutFold, "0_prepGeneData", "rna_rnaseqDT.Rdata"))))
  
  if(any(colSums(exprDT) == 0)) {
    txt <- paste0("...... !!! WARNING: colSums with zero detected: ", sum(colSums(exprDT) == 0),"/", ncol(exprDT), "\n")
    printAndLog(txt, myfile = logFile)
  }
  
  ######################################################################################## CHECK4: when do calcNormFactors crashes [RSEM and raw only]  
  
  if(inputDataType == "raw" | inputDataType == "RSEM") {
    
    cpm_exprDT <- cpm(exprDT)
    # txt <- paste0(toupper(script_name), "> NA in cpm_exprDT: ", 
    #               sum(is.na(cpm_exprDT)), "/", dim(cpm_exprDT)[1]*dim(cpm_exprDT)[2], " (",
    #               round((sum(is.na(cpm_exprDT))/(dim(cpm_exprDT)[1]*dim(cpm_exprDT)[2]) * 100),2), "%)\n")
    # cat(txt)
    rowsToKeep <- rowSums(cpm_exprDT, na.rm=T) >= (minCpmRatio * ncol(exprDT))
    
    exprDT <- exprDT[rowsToKeep, ]
    
    
    my_group <- unlist(sapply(colnames(exprDT), function(x) {
      ifelse(x %in% samp1_id, "cond1", ifelse(x %in% samp2_id, "cond2", NA))
    }))
    
    seqDataTmp <- DGEList(exprDT, group=my_group, genes = rownames(exprDT))
    #calcNormFactors:
    #it is usual to apply scale normalization to RNA-seq read counts, e.g. TMM normalization
    seqDataTmp <- try(calcNormFactors(seqDataTmp))
    if(class(seqDataTmp) == "try-error") {
      seqData <- DGEList(exprDT, group=my_group, genes = rownames(exprDT))
      seqData <- calcNormFactors(seqData, method = "none")
      txt <- paste0("...... !!! WARNING: could not compute calcNormFactors with default method, used method = \"none\" \n")
      printAndLog(txt, myfile = logFile)
    }
  } 
} # end-iterating over setting files

cat(paste0("... written: ", logFile, "\n"))
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



  

