#!/usr/bin/Rscript

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: pipeline_regionList.Rdata
# - script0: pipeline_geneList.Rdata
# - script0: rna_fpkmDT.Rdata         # added 03.03.2019
# - script5: permutationsDT.Rdata
################################################################################

################  OUTPUT
# - wilcoxStat_permDT.Rdata
################################################################################

nCpu=70

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

source(file.path(pipScriptDir, "my_wilcox_onlyW.R"))


script0_name <- "0_prepGeneData"
script5_name <- "5_runPermutationsMedian"
script_name <- "6v2_runPermutationsWilcoxStat"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

wilcoxPaired <- TRUE
wilcoxAlternative <- "two.sided"

source("main_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(ifelse(SSHFS, 2, nCpu)) # from main_settings.R

# if microarray was not set in the settings file -> by default set to  FALSE
if(!exists("microarray")) microarray <- FALSE

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
dir.create(curr_outFold, recursive = TRUE)

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
file.remove(pipLogFile)

# ADDED 27.11.2018 to check using other files
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)

################################****************************************************************************************
####################################################### PREPARE INPUT
################################****************************************************************************************

# INPUT DATA
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

#*******************************************************************************
if(microarray) {
  ### CHANGED HERE FOR THE V2 VERSION !!! 02.03.2019
  # norm_rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_madnorm_rnaseqDT.Rdata"))))
  fpkm_rnaseqDT <- eval(parse(text = load(file.path(pipOutFold, script0_name, "rna_fpkmDT.Rdata"))))
}else {
  ### CHANGED HERE FOR THE V2 VERSION !!! 02.03.2019
  # norm_rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_qqnorm_rnaseqDT.Rdata")))) 
  fpkm_rnaseqDT <- eval(parse(text = load(file.path(pipOutFold, script0_name, "rna_fpkmDT.Rdata"))))
}


########## LOAD THE DATA
samp1 <- eval(parse(text=load(paste0(setDir, "/", sample1_file))))
samp2 <- eval(parse(text=load(paste0(setDir, "/", sample2_file))))

stopifnot(all(samp1 %in% colnames(fpkm_rnaseqDT)))
stopifnot(all(samp2 %in% colnames(fpkm_rnaseqDT)))

cat("... load permutation data ...\n")
permutationsDT <- eval(parse(text = load(paste0(pipOutFold, "/", script5_name, "/permutationsDT.Rdata"))))
if(ncol(permutationsDT) != nRandomPermut)
  stop("! NEED TO CHECK: different settings were used for running the permutations !\n")

pipeline_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))
if(!setequal(pipeline_geneList, rownames(permutationsDT))) {
  txtWarningGene <- paste0(toupper(script_name), "> Not the same set of genes in permutDT and pipeline_geneList\n")
  printAndLog(txtWarningGene, pipLogFile)    
  stop(txtWarningGene)
} else{
  txtWarningGene <- ""
}

gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% pipeline_geneList,]
stopifnot(nrow(gene2tadDT)  > 0)

pipeline_regionList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_regionList.Rdata"))))
if(useTADonly) {
  if(any(grepl("_BOUND", pipeline_regionList))) {
    stop("! data were not prepared for \"useTADonly\" !")
  }
}
if(!setequal(pipeline_regionList, permutationsDT[,2])) {
  txtWarningRegion <- paste0(toupper(script_name), "> Not the same set of regions in permutDT and pipeline_regionList\n")
  printAndLog(txtWarningRegion, pipLogFile)    
  stop(txtWarningRegion)
} else {
  txtWarningRegion <- ""
}

### TAKE ONLY THE GENES FOR WHICH A LOG FC VALUE IS AVAILABLE (I.E. THE ONES USED FOR DE ANALYSIS)
initNrow <- nrow(permutationsDT)

all_regions <- sort(unique(as.character(permutationsDT[,1])))

stopifnot(setequal(rownames(permutationsDT), pipeline_geneList))

if(useTADonly) {
  if(! all (regexpr("_TAD",  gene2tadDT$region[gene2tadDT$entrezID %in% rownames(permutationsDT)]) > 0 )) {
    stop("make not sense to filter TAD genes after permutations if permutations were run with genes belonging to non-TAD regions\n")
  }
  initLen <- length(all_regions)
  all_regions <- all_regions[grep("_TAD", all_regions)]
  # if want to use only TADs, make more sens if permutations were run without the TADs
  txt <- paste0(toupper(script_name), "> Take only TAD regions: ", length(all_regions), "/", initLen, "\n")
  printAndLog(txt, pipLogFile)    
  if(length(all_regions) < initLen){
    stop("make not sense to filter TAD regions after permutations if permutations were run with genes belonging to non-TAD regions\n")
  }
}

################################****************************************************************************************
######################################################## COMPUTE MEAN LOG FC BY TAD FOR THE PERMUTATIONS
################################****************************************************************************************

### REGIONS ARE STORED IN ROWNAMES OF PERMUTdt

cat("... start Wilcoxon permutDT \n")

wilcoxStat_permDT <- foreach(i_col = 1:ncol(permutationsDT), .combine='cbind') %dopar% {

  # wilcoxStat_permDT <- foreach(i_col = 1:2, .combine='cbind') %do% {
    
  g2t_permDT <- data.frame(entrezID = rownames(permutationsDT), 
                           region = permutationsDT[,i_col], stringsAsFactors = F)
  
  curr_regions <- unique(as.character(g2t_permDT$region))
  stopifnot(setequal(curr_regions, all_regions))
  stopifnot(length(curr_regions) == length(all_regions))
  
  cat("... start Wilcoxon tests permut ", i_col, "\t", ncol(permutationsDT) , "\n")
  # timed:   
  # sappply version: 1:20
  # utilisateur     système      écoulé 
  # 41.178       0.403      46.775 
  # %do% version: 1:20
  # utilisateur     système      écoulé 
  # 38.342       0.240      40.391 
  
  wilcox_pairedTAD_meanExpr_fpkm <- foreach(i_reg = c(1:length(curr_regions)), .combine='cbind') %do% {
    cat(paste0("...... wilcox tests, permut i_col: ", i_col, " - region ", i_reg, "/", length(curr_regions), "\n"))
    reg <- curr_regions[i_reg]
    perm_genes <- as.character(g2t_permDT$entrezID[g2t_permDT$region == reg])
    
    subData <- as.data.frame(fpkm_rnaseqDT[which(pipeline_geneList[rownames(fpkm_rnaseqDT)] %in% perm_genes),,drop=F])
    cond1_DT <-  subData[, samp1,drop=F]
    cond2_DT <- subData[, samp2,drop=F]
    
    reg_genes_avgExpr_cond1 <- rowMeans(cond1_DT)
    reg_genes_avgExpr_cond2 <- rowMeans(cond2_DT)
    stopifnot(!is.na(reg_genes_avgExpr_cond1))
    stopifnot(!is.na(reg_genes_avgExpr_cond2))
    stopifnot(names(reg_genes_avgExpr_cond1) == names(reg_genes_avgExpr_cond2) )
    
    wTest <- my_wilcox_onlyW.test(reg_genes_avgExpr_cond1, reg_genes_avgExpr_cond2, 
                         paired=wilcoxPaired, alternative=wilcoxAlternative,
                         exact=FALSE) # don't need for the pval
    as.numeric(wTest$statistic)
  } # end-wilcox stat for current permut
  
  stopifnot(length(wilcox_pairedTAD_meanExpr_fpkm) == length(all_regions))
  names(wilcox_pairedTAD_meanExpr_fpkm) <- curr_regions
  stopifnot(setequal(curr_regions, all_regions ))
  cat(paste0("... end Wilcoxon tests\n"))
  wilcox_pairedTAD_meanExpr_fpkm[all_regions]
} #end foreach all permut
                                                  # stopifnot(ncol(wilcoxStat_permDT) == ncol(permutationsDT)) # <<<<< WILL NEED TO  UPDATE
rownames(wilcoxStat_permDT) <- all_regions
colnames(wilcoxStat_permDT) <- paste0("permutation", 1:ncol(wilcoxStat_permDT))

cat("... end Wilcoxon tests permutDT \n")
outFile <- file.path(curr_outFold, "wilcoxStat_permDT.Rdata")
save(wilcoxStat_permDT, file= outFile)
cat(paste0("... written: ", outFile, "\n"))

################################****************************************************************************************
####################################################### WRITE OUTPUT
################################****************************************************************************************

txt <- paste0(toupper(script_name), "> Number of permutations for which wilcoxStat computed: ", ncol(wilcoxStat_permDT), "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> Number of regions for which wilcoxStat computed: ", nrow(wilcoxStat_permDT), "\n")
printAndLog(txt, pipLogFile)

cat(paste0(txtWarningGene))
cat(paste0(txtWarningRegion))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))




