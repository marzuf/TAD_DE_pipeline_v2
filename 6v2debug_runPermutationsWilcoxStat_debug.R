options(scipen=100)


SSHFS <- T
setDir <- ifelse(SSHFS, "/media/electron", "")
settingF <- file.path(setDir, "/mnt/etemp/marie/Cancer_HiC_data_TAD_DA/PIPELINE/INPUT_FILES/ENCSR401TBQ_Caki2_40kb/run_settings_TCGAkich_norm_kich.R")
pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

script0_name <- "0_prepGeneData"
script5_name <- "5_runPermutationsMedian"
script_name <- "6v2_runPermutationsWilcoxStat"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

wilcoxPaired <- TRUE
wilcoxAlternative <- "two.sided"

source(file.path(paste0(pipScriptDir, "_TopDom"), "main_settings.R"))
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(ifelse(SSHFS, 2, nCpu)) # from main_settings.R

# if microarray was not set in the settings file -> by default set to  FALSE
if(!exists("microarray")) microarray <- FALSE

pipLogFile <- ""

pipOutFold <- file.path(setDir, pipOutFold)

################################****************************************************************************************
####################################################### PREPARE INPUT
################################****************************************************************************************

# INPUT DATA
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

#*******************************************************************************

fpkm_rnaseqDT <- eval(parse(text = load(file.path( pipOutFold, script0_name, "rna_fpkmDT.Rdata"))))

########## LOAD THE DATA
samp1 <- eval(parse(text=load(paste0(setDir, "/", sample1_file))))
samp2 <- eval(parse(text=load(paste0(setDir, "/", sample2_file))))

stopifnot(all(samp1 %in% colnames(fpkm_rnaseqDT)))
stopifnot(all(samp2 %in% colnames(fpkm_rnaseqDT)))

cat("... load permutation data ...\n")

load(file.path(setDir, "/mnt/etemp/marie/Cancer_HiC_data_TAD_DA", "caki2_norm_kich_foo_5_permutDT.Rdata"))
permutationsDT <- caki2_norm_kich_foo_5_permutDT


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

gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% rownames(permutationsDT),]
stopifnot(nrow(gene2tadDT) == nrow(permutationsDT))

fpkm_rnaseqDT <- fpkm_rnaseqDT[rownames(fpkm_rnaseqDT) %in% as.character(gene2tadDT$entrezID), ]
rowMeans_samp1 <- rowMeans(fpkm_rnaseqDT[,samp1 ])
rowMeans_samp2 <- rowMeans(fpkm_rnaseqDT[,samp2])
stopifnot(names(rowMeans_samp1) == names(rowMeans_samp2))

################################****************************************************************************************
######################################################## COMPUTE MEAN LOG FC BY TAD FOR THE PERMUTATIONS
################################****************************************************************************************

### REGIONS ARE STORED IN ROWNAMES OF PERMUTdt

cat("... start Wilcoxon permutDT \n")

i_col=1

wilcoxStat_permDT <- foreach(i_col = 1:ncol(permutationsDT), .combine='cbind') %dopar% {

  curr_regions <- unique(as.character(permutationsDT[,i_col]))
  stopifnot(setequal(curr_regions, all_regions))
  
  # cat("... start Wilcoxon tests permut ", i_col, "\t", ncol(permutationsDT) , "\n")
  # wilcox_pairedTAD_meanExpr_fpkm <- foreach(i_reg = 1:length(curr_regions), .combine='c') %do% {
  #   cat(paste0("...... wilcox tests, region ", i_reg, "/", length(curr_regions), "\n"))
  #   reg <- curr_regions[i_reg]
  #   reg_genes <- gene2tadDT$entrezID[gene2tadDT$region == reg]
  #   subData <- as.data.frame(fpkm_rnaseqDT[which(pipeline_geneList[rownames(fpkm_rnaseqDT)] %in% reg_genes),,drop=F])
  #   
  #   cond1_DT <-  as.data.frame(fpkm_rnaseqDT[which(pipeline_geneList[rownames(fpkm_rnaseqDT)] %in% reg_genes), samp1,drop=F])
  #   cond2_DT <- as.data.frame(fpkm_rnaseqDT[which(pipeline_geneList[rownames(fpkm_rnaseqDT)] %in% reg_genes), samp2,drop=F])
  #   
  #   reg_genes_avgExpr_cond1 <- rowMeans(cond1_DT)
  #   reg_genes_avgExpr_cond2 <- rowMeans(cond2_DT)
  #   stopifnot(!is.na(reg_genes_avgExpr_cond1))
  #   stopifnot(!is.na(reg_genes_avgExpr_cond2))
  #   stopifnot(names(reg_genes_avgExpr_cond1) == names(reg_genes_avgExpr_cond2) )
  #   
  #   wTest <- wilcox.test(reg_genes_avgExpr_cond1, reg_genes_avgExpr_cond2, paired=TRUE, alternative="two.sided") 
  #   as.numeric(wTest$statistic)
  # } # end-wilcox stat for current permut
  
  cat("... start Wilcoxon tests permut ", i_col, "\t", ncol(permutationsDT) , "\n")
  
  i_reg=1
  
  
  system.time(wilcox_pairedTAD_meanExpr_fpkm <- sapply(c(1:length(curr_regions))[1:20], function(i_reg) {
    cat(paste0("...... wilcox tests, region ", i_reg, "/", length(curr_regions), "\n"))
    reg <- curr_regions[i_reg]
    reg_genes <- gene2tadDT$entrezID[gene2tadDT$region == reg]
    subData <- as.data.frame(fpkm_rnaseqDT[which(pipeline_geneList[rownames(fpkm_rnaseqDT)] %in% reg_genes),,drop=F])
    
    cond1_DT <-  as.data.frame(fpkm_rnaseqDT[which(pipeline_geneList[rownames(fpkm_rnaseqDT)] %in% reg_genes), samp1,drop=F])
    cond2_DT <- as.data.frame(fpkm_rnaseqDT[which(pipeline_geneList[rownames(fpkm_rnaseqDT)] %in% reg_genes), samp2,drop=F])
    
    reg_genes_avgExpr_cond1 <- rowMeans(cond1_DT)
    reg_genes_avgExpr_cond2 <- rowMeans(cond2_DT)
    stopifnot(!is.na(reg_genes_avgExpr_cond1))
    stopifnot(!is.na(reg_genes_avgExpr_cond2))
    stopifnot(names(reg_genes_avgExpr_cond1) == names(reg_genes_avgExpr_cond2) )
    
    wTest <- wilcox.test(reg_genes_avgExpr_cond1, reg_genes_avgExpr_cond2, paired=wilcoxPaired, alternative=wilcoxAlternative) 
    as.numeric(wTest$statistic)
  })) # end-wilcox stat for current permut
  
  # sappply version: 1:20
  # utilisateur     système      écoulé 
  # 41.178       0.403      46.775 
  
  # %do% version: 1:20
  # utilisateur     système      écoulé 
  # 38.342       0.240      40.391 
  
  
  # sappply version: 1:20 - wTest
  utilisateur     système      écoulé 
  52.632       0.434      60.116 
  
  # sappply version: 1:20 - wTest2
  utilisateur     système      écoulé 
  41.246       0.275      43.531 
  
  # sappply version: 1:20 - wTest2 - changed subset
  utilisateur     système      écoulé 
  15.141       0.109      15.392 
  
  user  system elapsed 
  0.203   0.000   0.203 

  system.time(wilcox_pairedTAD_meanExpr_fpkm <- foreach(i_reg=c(1:length(curr_regions))[1:20], .combine='c') %do% {
    cat(paste0("...... wilcox tests, region ", i_reg, "/", length(curr_regions), "\n"))
    reg <- curr_regions[i_reg]
    reg_genes <- gene2tadDT$entrezID[gene2tadDT$region == reg]
    subData <- as.data.frame(fpkm_rnaseqDT[which(pipeline_geneList[rownames(fpkm_rnaseqDT)] %in% reg_genes),,drop=F])
    
    cond1_DT <-  as.data.frame(fpkm_rnaseqDT[which(pipeline_geneList[rownames(fpkm_rnaseqDT)] %in% reg_genes), samp1,drop=F])
    cond2_DT <- as.data.frame(fpkm_rnaseqDT[which(pipeline_geneList[rownames(fpkm_rnaseqDT)] %in% reg_genes), samp2,drop=F])
    
    reg_genes_avgExpr_cond1 <- rowMeans(cond1_DT)
    reg_genes_avgExpr_cond2 <- rowMeans(cond2_DT)
    stopifnot(!is.na(reg_genes_avgExpr_cond1))
    stopifnot(!is.na(reg_genes_avgExpr_cond2))
    stopifnot(names(reg_genes_avgExpr_cond1) == names(reg_genes_avgExpr_cond2) )
    
    wTest <- wilcox.test(reg_genes_avgExpr_cond1, reg_genes_avgExpr_cond2, paired=wilcoxPaired, alternative=wilcoxAlternative) 
    as.numeric(wTest$statistic)
  }) # end-wilcox stat for current permut
  
  
  
  system.time(wilcox_pairedTAD_meanExpr_fpkm <- foreach(i_reg=c(1:length(curr_regions))[1:20], .combine='c') %do% {
    cat(paste0("...... wilcox tests, region ", i_reg, "/", length(curr_regions), "\n"))
    reg <- curr_regions[i_reg]
    reg_genes <- gene2tadDT$entrezID[gene2tadDT$region == reg]
    subData <- as.data.frame(fpkm_rnaseqDT[which(pipeline_geneList[rownames(fpkm_rnaseqDT)] %in% reg_genes),,drop=F])
    
    cond1_DT <- subData[, samp1,drop=F]
    cond2_DT <- subData[, samp2,drop=F]
    
    reg_genes_avgExpr_cond1 <- rowMeans(cond1_DT)
    reg_genes_avgExpr_cond2 <- rowMeans(cond2_DT)
    stopifnot(!is.na(reg_genes_avgExpr_cond1))
    stopifnot(!is.na(reg_genes_avgExpr_cond2))
    stopifnot(names(reg_genes_avgExpr_cond1) == names(reg_genes_avgExpr_cond2) )
    
    wTest <- wilcox.test(reg_genes_avgExpr_cond1, reg_genes_avgExpr_cond2, paired=wilcoxPaired, alternative=wilcoxAlternative) 
    as.numeric(wTest$statistic)
  }) # end-wilcox stat for current permut
  
  
  
  source(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/my_wilcox_onlyW.R"))
  system.time(wilcox_pairedTAD_meanExpr_fpkm2 <- foreach(i_reg=c(1:length(curr_regions))[1:20], .combine='c') %do% {
    cat(paste0("...... wilcox tests, region ", i_reg, "/", length(curr_regions), "\n"))
    reg <- curr_regions[i_reg]
    reg_genes <- gene2tadDT$entrezID[gene2tadDT$region == reg]
    subData <- as.data.frame(fpkm_rnaseqDT[which(pipeline_geneList[rownames(fpkm_rnaseqDT)] %in% reg_genes),,drop=F])
    
    cond1_DT <- subData[, samp1,drop=F]
    cond2_DT <- subData[, samp2,drop=F]
    
    reg_genes_avgExpr_cond1 <- rowMeans(cond1_DT)
    reg_genes_avgExpr_cond2 <- rowMeans(cond2_DT)
    stopifnot(!is.na(reg_genes_avgExpr_cond1))
    stopifnot(!is.na(reg_genes_avgExpr_cond2))
    stopifnot(names(reg_genes_avgExpr_cond1) == names(reg_genes_avgExpr_cond2) )
    
    wTest2 <- my_wilcox_onlyW.test(reg_genes_avgExpr_cond1, reg_genes_avgExpr_cond2, paired=wilcoxPaired, alternative=wilcoxAlternative) 
    
    as.numeric(wTest2$statistic)
    
    
  }) # end-wilcox stat for current permut
  
  
  
  source(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/my_wilcox_onlyW.R"))
  system.time(wilcox_pairedTAD_meanExpr_fpkm3 <- foreach(i_reg=c(1:length(curr_regions))[1:1000], .combine='c') %do% {
    cat(paste0("...... wilcox tests, region ", i_reg, "/", length(curr_regions), "\n"))
    reg <- curr_regions[i_reg]
    reg_genes <- gene2tadDT$entrezID[gene2tadDT$region == reg]
    stopifnot(reg_genes %in% names(rowMeans_samp1))
    reg_genes_avgExpr_cond1 <- rowMeans_samp1[reg_genes]
    reg_genes_avgExpr_cond2 <- rowMeans_samp2[reg_genes]
    stopifnot(!is.na(reg_genes_avgExpr_cond1))
    stopifnot(!is.na(reg_genes_avgExpr_cond2))
    stopifnot(names(reg_genes_avgExpr_cond1) == names(reg_genes_avgExpr_cond2) )
    wTest2 <- my_wilcox_onlyW.test(reg_genes_avgExpr_cond1, reg_genes_avgExpr_cond2, paired=wilcoxPaired, alternative=wilcoxAlternative) 
    as.numeric(wTest2$statistic)
  }) # end-wilcox stat for current permut
  user  system elapsed 
  0.493   0.047   0.534 
  user  system elapsed 
  0.220   0.002   0.219 
  
  source(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/my_wilcox_onlyW.R"))
  system.time(wilcox_pairedTAD_meanExpr_fpkm3 <- sapply(c(1:length(curr_regions))[1:1000],function(i_reg){
    cat(paste0("...... wilcox tests, region ", i_reg, "/", length(curr_regions), "\n"))
    reg <- curr_regions[i_reg]
    reg_genes <- gene2tadDT$entrezID[gene2tadDT$region == reg]
    stopifnot(reg_genes %in% names(rowMeans_samp1))
    reg_genes_avgExpr_cond1 <- rowMeans_samp1[reg_genes]
    reg_genes_avgExpr_cond2 <- rowMeans_samp2[reg_genes]
    stopifnot(!is.na(reg_genes_avgExpr_cond1))
    stopifnot(!is.na(reg_genes_avgExpr_cond2))
    stopifnot(names(reg_genes_avgExpr_cond1) == names(reg_genes_avgExpr_cond2) )
    wTest2 <- my_wilcox_onlyW.test(reg_genes_avgExpr_cond1, reg_genes_avgExpr_cond2, paired=wilcoxPaired, alternative=wilcoxAlternative) 
    as.numeric(wTest2$statistic)
    wTest2
  })) # end-wilcox stat for current permut
  user  system elapsed 
  0.419   0.038   0.453 
  user  system elapsed 
  0.200   0.000   0.198 
  
  
  
  stopifnot(length(wilcox_pairedTAD_meanExpr_fpkm) == length(all_regions))
  names(wilcox_pairedTAD_meanExpr_fpkm) <- curr_regions
  stopifnot(setequal(curr_regions, all_regions ))
  cat(paste0("... end Wilcoxon tests\n"))
  wilcox_pairedTAD_meanExpr_fpkm[all_regions]
} #end foreach all permut
stopifnot(ncol(wilcoxStat_permDT) == ncol(permutationsDT))
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




