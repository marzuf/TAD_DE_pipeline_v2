#!/usr/bin/Rscript

startTime <- Sys.time()

### 03.03.2019
# => v2 version: combine meanCorr v2 and wilcoxStat v2 !!!

# -> wilcoxStat replaces wilcoxStat
# -> meanCorrv2 uses another permutation scheme !

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script2v2: wilcox_pairedTAD_meanExpr_wilcoxStat.Rdata.Rdata
# - script8: all_obs_ratioDown.Rdata
# - script9v2: emp_pval_wilcoxStat.Rdata
# - script10v2: emp_pval_meanCorr.Rdata
################################################################################

################  OUTPUT
# - emp_pval_combined.Rdata + tables
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
script3_name <- "3_runMeanTADWilcoxStat"
script8_name <- "8c_runAllDown"
script9_name <- "9v2_runEmpPvalWilcoxStat"
script10_name <- "10v2_runEmpPvalMeanTADCorr"
script_name <- "11v2_runEmpPvalCombined"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar

# create the directories
curr_outFold <- file.path(pipOutFold,  script_name)
dir.create(curr_outFold, recursive = TRUE)

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
file.remove(pipLogFile)

# ADDED 16.11.2018 to check using other files
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)

################################****************************************************************************************
####################################################### PREPARE INPUT
################################****************************************************************************************
# load emp. p-val wilcoxStat 
emp_pval_wilcoxStat <- eval(parse(text = load(paste0(pipOutFold, "/", script9_name, "/emp_pval_meanWilcoxStat.Rdata"))))

# load emp. p-val intraTAD corr
emp_pval_intraCorr <- eval(parse(text = load(paste0(pipOutFold, "/", script10_name, "/emp_pval_meanCorr.Rdata"))))

intersectRegions <- sort(intersect(names(emp_pval_wilcoxStat), names(emp_pval_intraCorr)))
txt <- paste0(toupper(script_name), "> Take regions in common between permutations and observed data \n")
printAndLog(txt, pipLogFile)

### filter TAD only regions
if(useTADonly) {
    if( length(grep("_TAD", names(emp_pval_wilcoxStat))) > 0 ) {
        txt <- paste0(toupper(script_name), "> !!! WARNING: empirical p-val wilcoxStat data contain non-TAD regions as well !!!\n")
        printAndLog(txt, pipLogFile)    
    }
    if( length(grep("_TAD", names(emp_pval_intraCorr))) > 0 ) {
        txt <- paste0(toupper(script_name), "> !!! WARNING: empirical p-val intraCorr data contain non-TAD regions as well !!!\n")
        printAndLog(txt, pipLogFile)    
    }
    initLen <- length(intersectRegions)
    intersectRegions <- intersectRegions[grep("_TAD", intersectRegions)]
    txt <- paste0(toupper(script_name), "> Take only TAD regions: ", length(intersectRegions), "/", initLen, "\n")
    printAndLog(txt, pipLogFile)    
}

initLen <- length(emp_pval_wilcoxStat)
emp_pval_wilcoxStat <- emp_pval_wilcoxStat[intersectRegions]
txt <- paste0(toupper(script_name), "> ... -> emp. p-val wilcoxStat: ", length(emp_pval_wilcoxStat), "/", initLen, "\n")
printAndLog(txt, pipLogFile)

initLen <- length(emp_pval_intraCorr)
emp_pval_intraCorr <- emp_pval_intraCorr[intersectRegions]
txt <- paste0(toupper(script_name), "> ... -> emp. p-val intraCorr: ", length(emp_pval_intraCorr), "/", initLen, "\n")
printAndLog(txt, pipLogFile)

stopifnot(!any(is.na(emp_pval_wilcoxStat)))
stopifnot(!any(is.na(emp_pval_intraCorr)))

stopifnot(all(names(emp_pval_wilcoxStat) == names(emp_pval_intraCorr)))

################################****************************************************************************************
####################################################### CALCULATE EMP. PVAL INTRA-TAD CORR & WRITE OUTPUT
################################****************************************************************************************

# COMBINE EMPIRICAL P-VALUES
emp_pval_combined <- unlist(sapply(seq_along(intersectRegions), function(x) 
                  stouffer(c(emp_pval_intraCorr[x], emp_pval_wilcoxStat[x]), two.tails = TRUE)))
names(emp_pval_combined) <- intersectRegions

stopifnot(length(emp_pval_combined) == length(intersectRegions))

save(emp_pval_combined, file= paste0(curr_outFold, "/emp_pval_combined.Rdata"))
cat(paste0("... written: ", curr_outFold, "/emp_pval_combined.Rdata", "\n"))

#***** build and save table
# TAD | meanFC | ratioDown | emp. p-val combined | genes list comma separated
# load emp. p-val wilcoxStat 
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
pipeline_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))
DE_table <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_topTable.Rdata"))))

entrezDT <- read.delim(entrezDT_file, header=T, stringsAsFactors = F)
entrezDT$entrezID <- as.character(entrezDT$entrezID)

DE_table <- DE_table[DE_table$genes %in% names(pipeline_geneList),]
stopifnot(nrow(DE_table) > 0)
stopifnot(all(DE_table$genes %in% names(pipeline_geneList)))

DE_entrez <- as.character(unlist(sapply(DE_table$genes, function(x) pipeline_geneList[x])))
gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% DE_entrez,]

stopifnot(all(DE_entrez %in% entrezDT$entrezID))

################################****************************************************************************************
####################################################### BUILD TABLES & WRITE OUTPUT
################################****************************************************************************************

obs_wilcoxStat <- eval(parse(text = load(paste0(pipOutFold, "/", script3_name, "/wilcox_pairedTAD_meanExpr_wilcoxStat.Rdata.Rdata"))))
obs_ratioDown <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/all_obs_ratioDown.Rdata"))))

interReg <- sort(Reduce(intersect, list(names(obs_wilcoxStat), names(obs_ratioDown), names(emp_pval_combined))))

txt <- paste0(toupper(script_name), "> Number of TADs in wilcoxStat: ", length(obs_wilcoxStat), "\n")
printAndLog(txt, pipLogFile)    
txt <- paste0(toupper(script_name), "> Number of TADs in ratioDown: ", length(obs_ratioDown), "\n")
printAndLog(txt, pipLogFile)    
txt <- paste0(toupper(script_name), "> Number of TADs in emp. p-val combined: ", length(emp_pval_combined), "\n")
printAndLog(txt, pipLogFile)    
txt <- paste0(toupper(script_name), "> Number of TADs in the intersect: ", length(interReg), "\n")
printAndLog(txt, pipLogFile)    

interReg_empPval <- emp_pval_combined[names(emp_pval_combined) %in% interReg]
interReg_empPval_sort <- sort(interReg_empPval)

pvalDT <- foreach(i_reg = 1:length(interReg_empPval_sort), .combine = 'rbind') %do% {
  reg <- names(interReg_empPval_sort)[i_reg]
  reg_genes_entrez <- gene2tadDT$entrezID[gene2tadDT$region == reg]
  reg_genes_symbol <- unlist(sapply(reg_genes_entrez, function(x) entrezDT$symbol[entrezDT$entrezID == x]))
  reg_genes_symbol_list <- paste0(reg_genes_symbol, collapse = ",")
  data.frame(rank_pval = i_reg,
             TAD = reg,
             meanWilcoxStat = obs_wilcoxStat[reg],
             ratioDown = obs_ratioDown[reg],
             emp_pval_comb = emp_pval_combined[reg],
             TAD_genes = reg_genes_symbol_list
             )
}

outFile <- file.path(curr_outFold, "TAD_ratioDown_wilcoxStat_empPvalComb.txt")
write.table(pvalDT, file = outFile, col.names =T , row.names = F, quote=F, sep="\t")
cat(paste0("... written: ", outFile, "\n"))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)
cat(paste0("*** DONE: ", script_name, "\n"))
