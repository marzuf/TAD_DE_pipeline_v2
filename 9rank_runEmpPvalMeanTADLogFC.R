#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script3: all_meanLogFC_TAD.Rdata
# - script6: meanLogFC_permDT.Rdata
################################################################################

################  OUTPUT
# - emp_pval_meanLogFC_rank.Rdata + plots
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
script_name <- "9rank_runEmpPvalMeanTADLogFC"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(ifelse(SSHFS, 2, nCpu)) # from main_settings.R

tiesMet <- "min"

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
txt <- paste0("inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("tiesMet\t=\t", tiesMet, "\n")
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
obs_TADLogFC_rank <- rank(-abs(obs_TADLogFC), ties = tiesMet)
txt <- paste0(toupper(script_name), "> ... -> observed: ", length(obs_TADLogFC), "/", initLen, "\n")
printAndLog(txt, pipLogFile)

initNrow <- nrow(permut_TADLogFC_DT)
permut_TADLogFC_DT <- permut_TADLogFC_DT[intersectRegions,]
txt <- paste0(toupper(script_name), "> ... -> permutations: ", nrow(permut_TADLogFC_DT), "/", initNrow, "\n")
printAndLog(txt, pipLogFile)

stopifnot(is.numeric(permut_TADLogFC_DT[1,1]))

# FC -> consider up- and down- equally -> take the abs value
permut_TADLogFC_DT <- abs(permut_TADLogFC_DT)

# take -x to have decreasing ranking
rank_meanFC_permDT <- apply(permut_TADLogFC_DT, 2 , function(x) rank(-x, ties=tiesMet))
max1_idx <- which(rownames(permut_TADLogFC_DT) == names(which(rank_meanFC_permDT[,1] == 1)))
stopifnot(permut_TADLogFC_DT[-max1_idx,1] <  permut_TADLogFC_DT[max1_idx,1])

stopifnot( setequal(intersectRegions, rownames(rank_meanFC_permDT) ))





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


emp_pval_meanLogFC_rank <- unlist(foreach(reg = intersectRegions, .combine='c') %dopar% {
  # select the random values for this region
  all_shuff_meanFC_rank <- rank_meanFC_permDT[paste0(reg),]
  stopifnot(length(all_shuff_meanFC_rank) == ncol(rank_meanFC_permDT))
  # select the observed value for this region
  obs_rank <- obs_TADLogFC_rank[reg]
  stopifnot(length(obs_rank) == 1 )
  # the number of times the permut rank is smaller than the observed one
  # => the pvalue will be high if the permutation often a smaller rank
  emp_pval <- sum(all_shuff_meanFC_rank <= obs_rank) 
  # emp_pval/length(shuff_intraTADcorr_region)
  ### ADD THE +1 => THIS IS FOR THE OBSERVED VALUE -> AVOID THE 0 VALUES !!!
  (emp_pval+1)/(length(all_shuff_meanFC_rank)+1)
})
names(emp_pval_meanLogFC_rank) <- intersectRegions

stopifnot(all(emp_pval_meanLogFC_rank > 0 & emp_pval_meanLogFC_rank <= 1 ))

################################****************************************************************************************
####################################################### WRITE OUTPUT
################################****************************************************************************************

save(emp_pval_meanLogFC_rank, file= paste0(curr_outFold, "/emp_pval_meanLogFC_rank.Rdata"))
cat(paste0("... written: ", curr_outFold, "/emp_pval_meanLogFC_rank.Rdata", "\n"))


outFile <- paste0(curr_outFold, "/", "volcano_plot_empPval_meanLogFC.", plotType)
do.call(plotType, list(outFile, height=ifelse(plotType=="png", 480, 7), width=ifelse(plotType=="png", 600, 10)))

plot( x = obs_TADLogFC[intersectRegions], y = log10(emp_pval_meanLogFC_rank[intersectRegions]),
      pch=16, cex  = 0.7,
      main = paste0("volcano plot - emp. p-val mean logFC"),
      xlab="mean observed logFC",
      ylab = "log10 (# permuted mean logFC <> observed mean log FC)")

emp_pval_meanLogFC_rank <- emp_pval_meanLogFC_rank[order(emp_pval_meanLogFC_rank)]
top10_x <- names(sort(obs_TADLogFC[intersectRegions], decreasing = T))[1:nTop]
top10_y <- names(sort(log10(emp_pval_meanLogFC_rank[intersectRegions])))[1:nTop]
tads_to_plot <- union(top10_x, top10_y)
text(x = obs_TADLogFC[tads_to_plot],
     y = log10(emp_pval_meanLogFC_rank[tads_to_plot]),
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

