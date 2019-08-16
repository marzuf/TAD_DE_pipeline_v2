#!/usr/bin/Rscript

startTime <- Sys.time()


options(scipen=100)
### !!! UPDATE 24.07.19: STOUFFER ONE-SIDED

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script3: all_meanLogFC_TAD.Rdata
# - script8: all_obs_ratioDown.Rdata
# - script9rank: emp_pval_meanLogFC_rank.Rdata
# - script10rank: emp_pval_meanCorr_rank_rank.Rdata
################################################################################

################  OUTPUT
# - emp_pval_combined_rank.Rdata + tables
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
script8_name <- "8c_runAllDown"
script9rank_name <- "9rank_runEmpPvalMeanTADLogFC"
script10rank_name <- "10rank_runEmpPvalMeanTADCorr"
script_name <- "11_runEmpPvalCombined"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

# ADDED 16.11.2018 to check using other files
txt <- paste0("inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)

################################****************************************************************************************
####################################################### PREPARE INPUT
################################****************************************************************************************
# load emp. p-val logFC 
emp_pval_logFC_rank <- eval(parse(text = load(paste0(pipOutFold, "/", script9rank_name, "/emp_pval_meanLogFC_rank.Rdata"))))

# load emp. p-val intraTAD corr
emp_pval_intraCorr_rank <- eval(parse(text = load(paste0(pipOutFold, "/", script10rank_name, "/emp_pval_meanCorr_rank.Rdata"))))

intersectRegions <- sort(intersect(names(emp_pval_logFC_rank), names(emp_pval_intraCorr_rank)))
txt <- paste0(toupper(script_name), "> Take regions in common between permutations and observed data \n")
printAndLog(txt, pipLogFile)

### filter TAD only regions
if(useTADonly) {
    if( length(grep("_TAD", names(emp_pval_logFC_rank))) > 0 ) {
        txt <- paste0(toupper(script_name), "> !!! WARNING: empirical p-val logFC data contain non-TAD regions as well !!!\n")
        printAndLog(txt, pipLogFile)    
    }
    if( length(grep("_TAD", names(emp_pval_intraCorr_rank))) > 0 ) {
        txt <- paste0(toupper(script_name), "> !!! WARNING: empirical p-val intraCorr data contain non-TAD regions as well !!!\n")
        printAndLog(txt, pipLogFile)    
    }
    initLen <- length(intersectRegions)
    intersectRegions <- intersectRegions[grep("_TAD", intersectRegions)]
    txt <- paste0(toupper(script_name), "> Take only TAD regions: ", length(intersectRegions), "/", initLen, "\n")
    printAndLog(txt, pipLogFile)    
}

initLen <- length(emp_pval_logFC_rank)
emp_pval_logFC_rank <- emp_pval_logFC_rank[intersectRegions]
txt <- paste0(toupper(script_name), "> ... -> emp. p-val logFC: ", length(emp_pval_logFC_rank), "/", initLen, "\n")
printAndLog(txt, pipLogFile)

initLen <- length(emp_pval_intraCorr_rank)
emp_pval_intraCorr_rank <- emp_pval_intraCorr_rank[intersectRegions]
txt <- paste0(toupper(script_name), "> ... -> emp. p-val intraCorr: ", length(emp_pval_intraCorr_rank), "/", initLen, "\n")
printAndLog(txt, pipLogFile)

stopifnot(!any(is.na(emp_pval_logFC_rank)))
stopifnot(!any(is.na(emp_pval_intraCorr_rank)))

stopifnot(all(names(emp_pval_logFC_rank) == names(emp_pval_intraCorr_rank)))

################################****************************************************************************************
####################################################### CALCULATE EMP. PVAL INTRA-TAD CORR & WRITE OUTPUT
################################****************************************************************************************

# COMBINE EMPIRICAL P-VALUES
emp_pval_combined_rank <- unlist(sapply(seq_along(intersectRegions), function(x) 
                  stouffer(c(emp_pval_intraCorr_rank[x], emp_pval_logFC_rank[x]), two.tails = FALSE)))
names(emp_pval_combined_rank) <- intersectRegions

stopifnot(length(emp_pval_combined_rank) == length(intersectRegions))

save(emp_pval_combined_rank, file= paste0(curr_outFold, "/emp_pval_combined_rank.Rdata"))
cat(paste0("... written: ", curr_outFold, "/emp_pval_combined_rank.Rdata", "\n"))

#***** build and save table
# TAD | meanFC | ratioDown | emp. p-val combined | genes list comma separated
# load emp. p-val logFC 
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

obs_logFC <- eval(parse(text = load(paste0(pipOutFold, "/", script3_name, "/all_meanLogFC_TAD.Rdata"))))
obs_ratioDown <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/all_obs_ratioDown.Rdata"))))

interReg <- sort(Reduce(intersect, list(names(obs_logFC), names(obs_ratioDown), names(emp_pval_combined_rank))))

txt <- paste0(toupper(script_name), "> Number of TADs in logFC: ", length(obs_logFC), "\n")
printAndLog(txt, pipLogFile)    
txt <- paste0(toupper(script_name), "> Number of TADs in ratioDown: ", length(obs_ratioDown), "\n")
printAndLog(txt, pipLogFile)    
txt <- paste0(toupper(script_name), "> Number of TADs in emp. p-val combined: ", length(emp_pval_combined_rank), "\n")
printAndLog(txt, pipLogFile)    
txt <- paste0(toupper(script_name), "> Number of TADs in the intersect: ", length(interReg), "\n")
printAndLog(txt, pipLogFile)    

interReg_empPval <- emp_pval_combined_rank[names(emp_pval_combined_rank) %in% interReg]
interReg_empPval_sort <- sort(interReg_empPval)

pvalDT <- foreach(i_reg = 1:length(interReg_empPval_sort), .combine = 'rbind') %do% {
  reg <- names(interReg_empPval_sort)[i_reg]
  reg_genes_entrez <- gene2tadDT$entrezID[gene2tadDT$region == reg]
  reg_genes_symbol <- unlist(sapply(reg_genes_entrez, function(x) entrezDT$symbol[entrezDT$entrezID == x]))
  reg_genes_symbol_list <- paste0(reg_genes_symbol, collapse = ",")
  data.frame(rank_pval = i_reg,
             TAD = reg,
             meanLogFC = obs_logFC[reg],
             ratioDown = obs_ratioDown[reg],
             emp_pval_comb = emp_pval_combined_rank[reg],
             TAD_genes = reg_genes_symbol_list
             )
}

write.table(pvalDT, file = paste0(curr_outFold, "/TAD_ratioDown_logFC_empPvalComb.txt"), col.names =T , row.names = F, quote=F, sep="\t")
cat(paste0("... written: ", paste0(curr_outFold, "/TAD_ratioDown_logFC_empPvalComb.txt"), "\n"))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)
cat(paste0("*** DONE: ", script_name, "\n"))
