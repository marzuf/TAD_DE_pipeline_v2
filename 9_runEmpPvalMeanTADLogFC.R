#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script3: all_meanLogFC_TAD.Rdata
# - script6: meanLogFC_permDT.Rdata
################################################################################

################  OUTPUT
# - emp_pval_meanLogFC.Rdata + plots
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
script6_name <- "6_runPermutationsMeanLogFC"
script_name <- "9_runEmpPvalMeanTADLogFC"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(ifelse(SSHFS, 2, nCpu)) # from main_settings.R


# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

# ADDED 27.11.2018 to check using other files
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)


# ADDED 16.11.2018 to check using other files
txt <- paste0("inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
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
obs_TADLogFC <- eval(parse(text = load(paste0(pipOutFold, "/", script3_name, "/all_meanLogFC_TAD.Rdata"))))
initLen <- length(obs_TADLogFC)
obs_TADLogFC <- na.omit(obs_TADLogFC)
txt <- paste0(toupper(script_name), "> Discard rows with NA in observed logFC. Retain: ", length(obs_TADLogFC), "/", initLen, "\n")
printAndLog(txt, pipLogFile)

permut_TADLogFC_DT <- eval(parse(text = load(paste0(pipOutFold, "/", script6_name, "/meanLogFC_permDT.Rdata"))))
initNrow <- nrow(permut_TADLogFC_DT)
permut_TADLogFC_DT <- na.omit(permut_TADLogFC_DT)
txt <- paste0(toupper(script_name), "> Discard rows with NA in permutation logFC. Retain: ", nrow(permut_TADLogFC_DT), "/", initNrow, "\n")
printAndLog(txt, pipLogFile)

intersectRegions <- sort(intersect(names(obs_TADLogFC), rownames(permut_TADLogFC_DT)))
txt <- paste0(toupper(script_name), "> Take regions in common between permutations and observed data \n")
printAndLog(txt, pipLogFile)

initLen <- length(obs_TADLogFC)
obs_TADLogFC <- obs_TADLogFC[intersectRegions]
txt <- paste0(toupper(script_name), "> ... -> observed: ", length(obs_TADLogFC), "/", initLen, "\n")
printAndLog(txt, pipLogFile)

initNrow <- nrow(permut_TADLogFC_DT)
permut_TADLogFC_DT <- permut_TADLogFC_DT[intersectRegions,]
txt <- paste0(toupper(script_name), "> ... -> permutations: ", nrow(permut_TADLogFC_DT), "/", initNrow, "\n")
printAndLog(txt, pipLogFile)

stopifnot(is.numeric(permut_TADLogFC_DT[1,1]))

################################****************************************************************************************
####################################################### CALCULATE EMP. PVAL LOGFC
################################****************************************************************************************

### filter TAD only regions
if(useTADonly) {
    if( length(grep("_TAD", rownames(permut_TADLogFC_DT))) > 0 ) {
        txt <- paste0(toupper(script_name), "> !!! WARNING: logFC permutations data generated using non-TAD regions as well !!!\n")
        printAndLog(txt, pipLogFile)    
    }
    initLen <- length(intersectRegions)
    intersectRegions <- intersectRegions[grep("_TAD", intersectRegions)]
    txt <- paste0(toupper(script_name), "> Take only TAD regions: ", length(intersectRegions), "/", initLen, "\n")
    printAndLog(txt, pipLogFile)    
}


emp_pval_meanLogFC <- foreach(reg = intersectRegions, .combine='c') %dopar% {
  # select the random values for this region
  shuff_logFC_TAD_region <- permut_TADLogFC_DT[reg,]
  stopifnot(length(shuff_logFC_TAD_region) == ncol(permut_TADLogFC_DT))
  # select the observed value for this region
  obs_meanLogFC <- obs_TADLogFC[reg]
  stopifnot(length(obs_meanLogFC) == 1 )
  if(obs_meanLogFC > 0){
    emp_pval <- sum(shuff_logFC_TAD_region >= obs_meanLogFC)   
  } else if (obs_meanLogFC < 0) {
    emp_pval <- sum(shuff_logFC_TAD_region <= obs_meanLogFC)   
  } else {
    emp_pval <- mean(c(sum(shuff_logFC_TAD_region <= obs_meanLogFC), sum(shuff_logFC_TAD_region >= obs_meanLogFC)), na.rm=T)
    warning("obs log FC == 0, will take the mean \n")
  }
  # emp_pval/length(shuff_logFC_TAD_region)
  ### ADD THE +1 => THIS IS FOR THE OBSERVED VALUE -> AVOID THE 0 VALUES !!!
  (emp_pval+1)/(length(shuff_logFC_TAD_region)+1)
}
names(emp_pval_meanLogFC) <- intersectRegions

stopifnot(all(emp_pval_meanLogFC > 0 & emp_pval_meanLogFC <= 1 ))

################################****************************************************************************************
####################################################### WRITE OUTPUT
################################****************************************************************************************

save(emp_pval_meanLogFC, file= paste0(curr_outFold, "/emp_pval_meanLogFC.Rdata"))
cat(paste0("... written: ", curr_outFold, "/emp_pval_meanLogFC.Rdata", "\n"))


outFile <- paste0(curr_outFold, "/", "volcano_plot_empPval_meanLogFC.", plotType)
do.call(plotType, list(outFile, height=ifelse(plotType=="png", 480, 7), width=ifelse(plotType=="png", 600, 10)))

plot( x = obs_TADLogFC[intersectRegions], y = log10(emp_pval_meanLogFC[intersectRegions]),
      pch=16, cex  = 0.7,
      main = paste0("volcano plot - emp. p-val mean logFC"),
      xlab="mean observed logFC",
      ylab = "log10 (# permuted mean logFC <> observed mean log FC)")

emp_pval_meanLogFC_rank <- rank(emp_pval_meanLogFC, ties="min")
emp_pval_meanLogFC_rank <- emp_pval_meanLogFC_rank[order(emp_pval_meanLogFC_rank)]
top10_x <- names(sort(obs_TADLogFC[intersectRegions], decreasing = T))[1:nTop]
top10_y <- names(sort(log10(emp_pval_meanLogFC[intersectRegions])))[1:nTop]
tads_to_plot <- union(top10_x, top10_y)
text(x = obs_TADLogFC[tads_to_plot],
     y = log10(emp_pval_meanLogFC[tads_to_plot]),
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

