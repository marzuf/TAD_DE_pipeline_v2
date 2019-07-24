#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()


#### 03.03.2019
# -> v2 version, use of wilcox stat instead of logFC
# http://users.sussex.ac.uk/~grahamh/RM1web/WilcoxonHandoout2011.pdf
# The  Wilcoxon  test  statistic  "W"  is  simply  the  smaller  of  the  rank  totals.  
# The SMALLER  it  is  (taking  into  account  how  many  participants  you  have)  then  the  less likely it is to have occurred by chance. 
# A table of critical values of W shows you how likely  it  is  to  obtain  your  particular  value  of  W  purely  by  chance.  
# Note  that  the Wilcoxon  test  is  unusual  in  this  respect:  
#   normally,  the  BIGGER  the test statistic,  the less likely it is to have occurred by chance). 
# 
# For a two tailed test the test statistic is the smaller of W+ and W-
#   If the null hypothesis were true, you would expect W+  and W-  to have roughly the same value.  

# Freudenberg 2010: respective W statistic assuming that datasets with the smallest W statistic are the most similar.

# => emp. pval: the SMALLER the W, the LESS likely occured by chance
# => so the emp pval is # obs W <= permut W

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script2v2: wilcox
# - script6v2:wilcoxStat_permDT.Rdata
################################################################################

################  OUTPUT
# - emp_pval_wilcoxStat.Rdata + plots
################################################################################

nCpu <- 70

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script3_name <- "2v2_runWilcoxonTAD" # wilcox_pairedTAD_meanExpr_fpkm.Rdata
script6_name <- "6v2onlyW_runPermutationsWilcoxStat"
script_name <- "9v2_runEmpPvalWilcoxStat"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(ifelse(SSHFS, 2, nCpu)) # from main_settings.R

# create the directories
curr_outFold <- file.path(pipOutFold, script_name)
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

# ADDED 16.11.2018 to check using other files
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)

#****************************************************************************************************************************
plotType  <- "svg"
# annotate the plot for the nTop top ranking regions
nTop <- 10
#****************************************************************************************************************************

################################****************************************************************************************
####################################################### PREPARE INPUT
################################****************************************************************************************

# LOAD DATA AND DISCARD THE REGIONS WITH NA 
obs_TADwilcoxStat <- eval(parse(text = load(file.path(pipOutFold, script3_name, "wilcox_pairedTAD_meanExpr_wilcoxStat.Rdata"))))
initLen <- length(obs_TADwilcoxStat)
stopifnot(!is.na(obs_TADwilcoxStat))

permut_TADwilcoxStat_DT <- eval(parse(text = load(file.path(pipOutFold, script6_name, "wilcoxStat_permDT.Rdata"))))
initNrow <- nrow(permut_TADwilcoxStat_DT)
permut_TADwilcoxStat_DT <- na.omit(permut_TADwilcoxStat_DT)
txt <- paste0(toupper(script_name), "> Discard rows with NA in permutation Wilcox. Stat. Retain: ", nrow(permut_TADwilcoxStat_DT), "/", initNrow, "\n")
printAndLog(txt, pipLogFile)

intersectRegions <- sort(intersect(names(obs_TADwilcoxStat), rownames(permut_TADwilcoxStat_DT)))
txt <- paste0(toupper(script_name), "> Take regions in common between permutations and observed data \n")
printAndLog(txt, pipLogFile)

initLen <- length(obs_TADwilcoxStat)
obs_TADwilcoxStat <- obs_TADwilcoxStat[intersectRegions]
txt <- paste0(toupper(script_name), "> ... -> observed: ", length(obs_TADwilcoxStat), "/", initLen, "\n")
printAndLog(txt, pipLogFile)

initNrow <- nrow(permut_TADwilcoxStat_DT)
permut_TADwilcoxStat_DT <- permut_TADwilcoxStat_DT[intersectRegions,]
txt <- paste0(toupper(script_name), "> ... -> permutations: ", nrow(permut_TADwilcoxStat_DT), "/", initNrow, "\n")
printAndLog(txt, pipLogFile)

stopifnot(is.numeric(permut_TADwilcoxStat_DT[1,1]))

################################****************************************************************************************
####################################################### CALCULATE EMP. PVAL wilcoxStat
################################****************************************************************************************

### filter TAD only regions
if(useTADonly) {
    if( length(grep("_TAD", rownames(permut_TADwilcoxStat_DT))) > 0 ) {
        txt <- paste0(toupper(script_name), "> !!! WARNING: wilcoxStat permutations data generated using non-TAD regions as well !!!\n")
        printAndLog(txt, pipLogFile)    
    }
    initLen <- length(intersectRegions)
    intersectRegions <- intersectRegions[grep("_TAD", intersectRegions)]
    txt <- paste0(toupper(script_name), "> Take only TAD regions: ", length(intersectRegions), "/", initLen, "\n")
    printAndLog(txt, pipLogFile)    
}

emp_pval_wilcoxStat <- foreach(reg = intersectRegions, .combine='c') %dopar% {
  # select the random values for this region
  shuff_wilcoxStat_TAD_region <- permut_TADwilcoxStat_DT[reg,]
  stopifnot(length(shuff_wilcoxStat_TAD_region) == ncol(permut_TADwilcoxStat_DT))
  # select the observed value for this region
  obs_wilcoxStat <- obs_TADwilcoxStat[reg]
  stopifnot(length(obs_wilcoxStat) == 1 )
  
  # SMALLER W stat are less likely to occur
  # so the pval is the # the obs. W <= permut. W
  emp_pval <- sum(obs_wilcoxStat <= shuff_wilcoxStat_TAD_region )   

  # emp_pval/length(shuff_wilcoxStat_TAD_region)
  ### ADD THE +1 => THIS IS FOR THE OBSERVED VALUE -> AVOID THE 0 VALUES !!!
  (emp_pval+1)/(length(shuff_wilcoxStat_TAD_region)+1)
}
names(emp_pval_wilcoxStat) <- intersectRegions

stopifnot(all(emp_pval_wilcoxStat > 0 & emp_pval_wilcoxStat <= 1 ))

################################****************************************************************************************
####################################################### WRITE OUTPUT
################################****************************************************************************************

save(emp_pval_wilcoxStat, file= paste0(curr_outFold, "/emp_pval_wilcoxStat.Rdata"))
cat(paste0("... written: ", curr_outFold, "/emp_pval_wilcoxStat.Rdata", "\n"))

outFile <- paste0(curr_outFold, "/", "volcano_plot_empPval_wilcoxStat.", plotType)
do.call(plotType, list(outFile, height=ifelse(plotType=="png", 480, 7), width=ifelse(plotType=="png", 600, 10)))

plot( x = obs_TADwilcoxStat[intersectRegions], y = log10(emp_pval_wilcoxStat[intersectRegions]),
      pch=16, cex  = 0.7,
      main = paste0("volcano plot - emp. p-val Wilcox. Stat."),
      xlab="observed wilcoxStat",
      ylab = "log10 (# observed Wilcox. Stat. <= permuted Wilcox. Stat.)")

emp_pval_wilcoxStat_rank <- rank(emp_pval_wilcoxStat, ties="min")
emp_pval_wilcoxStat_rank <- emp_pval_wilcoxStat_rank[order(emp_pval_wilcoxStat_rank)]
top10_x <- names(sort(obs_TADwilcoxStat[intersectRegions], decreasing = T))[1:nTop]
top10_y <- names(sort(log10(emp_pval_wilcoxStat[intersectRegions])))[1:nTop]
tads_to_plot <- union(top10_x, top10_y)
text(x = obs_TADwilcoxStat[tads_to_plot],
     y = log10(emp_pval_wilcoxStat[tads_to_plot]),
     col ="red",
     pos=4,
     srt = 45,
     labels = tads_to_plot)
legend("bottomright", bty="n", legend=c(paste0("(Top", nTop, " - union x and y, unsorted)"), tads_to_plot), cex=0.8)  
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)
cat(paste0("*** DONE: ", script_name, "\n"))

