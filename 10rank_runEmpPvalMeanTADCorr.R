#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script4: all_meanCorr_TAD.Rdata
# - script7: meanCorr_permDT.Rdata
################################################################################

################  OUTPUT
# - emp_pval_meanCorr_rank.Rdata + plots
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
script7_name <- "7_runPermutationsMeanTADCorr"
script_name <- "10rank_runEmpPvalMeanTADCorr"
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

plotType  <- "svg"

tiesMet <- "min"

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
txt <- paste0("tiesMet\t=\t", tiesMet, "\n")
printAndLog(txt, pipLogFile)

#****************************************************************************************************************************
withDiago <- FALSE
txt <- paste0(toupper(script_name), "> Take the diagonal when computing the mean of lower triangle: ", as.character(withDiago), "\n")
printAndLog(txt, pipLogFile)

# annotate the plot for the nTop top ranking regions
nTop <- 10
#*******************************************************************************

################################****************************************************************************************
####################################################### PREPARE INPUT
################################****************************************************************************************
# LOAD DATA AND DISCARD THE REGIONS WITH NA 
obs_intraTADcorr <- eval(parse(text = load(paste0(pipOutFold, "/", script4_name, "/all_meanCorr_TAD.Rdata"))))
initLen <- length(obs_intraTADcorr)
obs_intraTADcorr <- na.omit(obs_intraTADcorr)
txt <- paste0(toupper(script_name), "> Discard rows with NA in observed intraTAD correlation. Retain: ", length(obs_intraTADcorr), "/", initLen, "\n")
printAndLog(txt, pipLogFile)

permut_intraTADcorr_DT <- eval(parse(text = load(paste0(pipOutFold, "/", script7_name, "/meanCorr_permDT.Rdata"))))
initNrow <- nrow(permut_intraTADcorr_DT)
permut_intraTADcorr_DT <- na.omit(permut_intraTADcorr_DT)
txt <- paste0(toupper(script_name), "> Discard rows with NA in permutation intraTAD correlation. Retain: ", nrow(permut_intraTADcorr_DT), "/", initNrow, "\n")
printAndLog(txt, pipLogFile)


intersectRegions <- sort(intersect(names(obs_intraTADcorr), rownames(permut_intraTADcorr_DT)))
txt <- paste0(toupper(script_name), "> Take regions in common between permutations and observed data \n")
printAndLog(txt, pipLogFile)

initLen <- length(obs_intraTADcorr)
obs_intraTADcorr <- obs_intraTADcorr[intersectRegions]
obs_intraTADcorr_rank <- rank(-obs_intraTADcorr, ties=tiesMet)
txt <- paste0(toupper(script_name), "> ... -> observed: ", length(obs_intraTADcorr), "/", initLen, "\n")
printAndLog(txt, pipLogFile)

initNrow <- nrow(permut_intraTADcorr_DT)
permut_intraTADcorr_DT <- permut_intraTADcorr_DT[intersectRegions,]
txt <- paste0(toupper(script_name), "> ... -> permutations: ", nrow(permut_intraTADcorr_DT), "/", initNrow, "\n")
printAndLog(txt, pipLogFile)

stopifnot( setequal(intersectRegions, rownames(permut_intraTADcorr_DT) ))

# take -x to have decreasing ranking
rank_meanCorr_permDT <- apply(permut_intraTADcorr_DT, 2 , function(x) rank(-x, ties=tiesMet))
max1_idx <- which(rownames(permut_intraTADcorr_DT) == names(which(rank_meanCorr_permDT[,1] == 1)))
stopifnot(permut_intraTADcorr_DT[-max1_idx,1] <  permut_intraTADcorr_DT[max1_idx,1])

stopifnot( setequal(intersectRegions, rownames(rank_meanCorr_permDT) ))



### filter TAD only regions
if(useTADonly) {
    if( length(grep("_TAD", rownames(permut_intraTADcorr_DT))) > 0 ) {
        txt <- paste0(toupper(script_name), "> !!! WARNING: intraTADcorr permutations data generated using non-TAD regions as well !!!\n")
        printAndLog(txt, pipLogFile)    
    }
    initLen <- length(intersectRegions)
    intersectRegions <- intersectRegions[grep("_TAD", intersectRegions)]
    txt <- paste0(toupper(script_name), "> Take only TAD regions: ", length(intersectRegions), "/", initLen, "\n")
    printAndLog(txt, pipLogFile)    
}

stopifnot(is.numeric(permut_intraTADcorr_DT[1,1]))

################################****************************************************************************************
####################################################### CALCULATE EMP. PVAL INTRA-TAD CORR
################################****************************************************************************************
  

stopifnot( setequal(intersectRegions, rownames(rank_meanCorr_permDT) ))

emp_pval_meanCorr_rank <- unlist(foreach(reg = intersectRegions, .combine='c') %dopar% {
  # select the random values for this region
  all_shuff_meanCorr_rank <- rank_meanCorr_permDT[paste0(reg),]
  stopifnot(length(all_shuff_meanCorr_rank) == ncol(rank_meanCorr_permDT))
  # select the observed value for this region
  obs_rank <- obs_intraTADcorr_rank[reg]
  stopifnot(length(obs_rank) == 1 )
  # the number of times the permut rank is smaller than the observed one
  # => the pvalue will be high if the permutation often a smaller rank
  emp_pval <- sum(all_shuff_meanCorr_rank <= obs_rank) 
  # emp_pval/length(shuff_intraTADcorr_region)
  ### ADD THE +1 => THIS IS FOR THE OBSERVED VALUE -> AVOID THE 0 VALUES !!!
  (emp_pval+1)/(length(all_shuff_meanCorr_rank)+1)
})
names(emp_pval_meanCorr_rank) <- intersectRegions

stopifnot(all(emp_pval_meanCorr_rank > 0 & emp_pval_meanCorr_rank <= 1 ))

################################****************************************************************************************
####################################################### WRITE OUTPUT
################################****************************************************************************************


save(emp_pval_meanCorr_rank, file= paste0(curr_outFold, "/emp_pval_meanCorr_rank.Rdata"))
cat(paste0("... written: ", curr_outFold, "/emp_pval_meanCorr_rank.Rdata", "\n"))


#### PLOT: X = MEAN OBSERVED CORRELATION; Y=log10 EMPIRICAL PVALUES FOR THE CORRELATION
outFile <- paste0(curr_outFold, "/", "volcano_plot_empPval_meanCorr.", plotType)
do.call(plotType, list(outFile, height=ifelse(plotType=="png", 480, 7), width=ifelse(plotType=="png", 600, 10)))
plot( x = obs_intraTADcorr[intersectRegions], y = log10(emp_pval_meanCorr_rank[intersectRegions]),
      pch=16, cex  = 0.7,
      main = paste0("volcano plot - emp. p-val intraTAD correlation"),
      xlab="mean observed intraTAD expression correlation",
      ylab = "log10 (# permuted mean expression correlation > observed mean correlation)")


emp_pval_meanCorr_rank <- emp_pval_meanCorr_rank[order(emp_pval_meanCorr_rank)]
top10_x <- names(sort(obs_intraTADcorr[intersectRegions], decreasing = T))[1:nTop]
top10_y <- names(sort(log10(emp_pval_meanCorr_rank[intersectRegions])))[1:nTop]
tads_to_plot <- union(top10_x, top10_y)
text(x = obs_intraTADcorr[tads_to_plot],
     y = log10(emp_pval_meanCorr_rank[tads_to_plot]),
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
