library(foreach)
library(doMC)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

registerDoMC(ifelse(SSHFS, 2, 20))

# script0_name <- "0_prepGeneData"
# script3_name <- "3_runMeanTADLogFC"
# script4_name <- "4_runMeanTADCorr"
# script6_name <- "6_runPermutationsMeanLogFC"
# script7_name <- "7_runPermutationsMeanTADCorr"
# script8_name <- "8c_runAllDown"

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

script19_name <- "19_SAM_emp_measurement"

cat(paste0("> Rscript cmp_19_SAM.R\n"))

# cat(paste0("setDir = ", setDir, "\n"))
# source("main_settings.R") # setDir is the main_settings not in run_settings
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

caller <- "TopDom"

mainOutFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER")

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "cmp_19_SAM")
system(paste0("mkdir -p ", outFold))

myHeight <- 10
myWidth <- 12

all_datasets <- list.files(mainOutFold, full.names = TRUE)


all_datasets_list <- foreach(curr_dataset = all_datasets) %dopar% {
  eval(parse(text = load(file.path(curr_dataset, script19_name, "empFDR_list.Rdata"))))
  # ratioLists <- eval(parse(text = load(file.path(curr_dataset, script19_name, "empFDR_list.Rdata"))))
  # all_datasets_list[[paste0(basename(curr_dataset))]] <- ratioLists
}
names(all_datasets_list) <- basename(all_datasets)


all_ratios <- c("logFC", "intraTADcorr", "prodSignedRatio", "ratioDown")

# all_ratios <- all_ratios[1]

for(curr_ratio in all_ratios) {
  
  cat("> START", curr_ratio, "\n")
  
  # retrieve the x-axis seq:
  xseq_FDR <- unique(lapply(all_datasets_list, function(x) names(x[[paste0("empFDR_", curr_ratio)]])))
  stopifnot(length(xseq_FDR) == 1)
  xseq_NbrSignif <- unique(lapply(all_datasets_list, function(x) names(x[[paste0("nbrSignif_", curr_ratio)]])))
  stopifnot(length(xseq_NbrSignif) == 1)
  stopifnot(xseq_FDR[[1]] == xseq_NbrSignif[[1]])
  xseq <- as.numeric(as.character(xseq_FDR[[1]]))
  stopifnot(!is.na(xseq))
  
  # retrieve all y-values for FDR
  y_FDR <- lapply(all_datasets_list, function(x) x[[paste0("empFDR_", curr_ratio)]])
  
  # retrieve all y-values for nbr signif
  y_NbrSignif <- lapply(all_datasets_list, function(x) x[[paste0("nbrSignif_", curr_ratio)]])
  
  y_FDR_range <- range(unlist(y_FDR)[!is.na(unlist(y_FDR)) & !is.infinite(unlist(y_FDR))]) * 100
  stopifnot(!is.na(y_FDR_range) & ! is.infinite(y_FDR_range))
  y_NbrSignif_range <- c(0, log10(max(unlist(y_NbrSignif)[!is.na(unlist(y_NbrSignif)) & !is.infinite(unlist(y_NbrSignif))])))
  stopifnot(!is.na(y_NbrSignif_range) & ! is.infinite(y_NbrSignif_range))
  
  rightAxisCol <- "steelblue"
  feature_name <- "TADs"
  variableName <- curr_ratio
  
  # plot_FDR_with_observedSignif <- function(yaxis_empFDR_vect, xaxis_cutoff, y2_obsSignif, variableName,  feature_name = "genes",) {
   # draw first empty plot
  
     
  outFile <- paste0(outFold, "/", curr_ratio, "_all_datasets_empFDR_nbrSignif.svg")   
  svg(outFile, height = myHeight, width = myWidth)   
  initmar <- par()$mar
  par(mar = c(5,5,2,5))
  
  plot(NULL,
       xlim = range(xseq),
       ylim = y_FDR_range,
       type="o", pch=16, 
       main=paste0("Empirical FDR with variable cut-off of ", variableName),
       xlab=paste0(variableName, " cut-off"),
       ylab="% FDR", axes=F)
  axis(1, at=xseq, labels = xseq)
  axis(2)
  
  for(i in 1:length(y_FDR)) {
    lines(x = xseq, y = 100*y_FDR[[i]], pch=16, type="o", cex = 0.6)
  }
                                            # plot((100*yaxis_empFDR_vect) ~ xaxis_cutoff, type="o", pch=16, 
                                            #      main=paste0("Empirical FDR with variable cut-off of ", variableName),
                                            #      xlab=paste0(variableName, " cut-off"),
                                            #      ylab="% FDR", axes=F)
  par(new = T)
  plot(NULL, 
       xlim = range(xseq),
       ylim = y_NbrSignif_range,
       type="o", 
       col=rightAxisCol, 
       pch=16, 
       axes=F, xlab=NA, ylab=NA, bty="l")
                                            # plot(log10(y2_obsSignif) ~ xaxis_cutoff, type="o", 
                                            #      col=rightAxisCol, pch=16, ylim=c(0, max(log10(y2_obsSignif))),
                                            #      axes=F, xlab=NA, ylab=NA, bty="l")
  box(bty="u")
  axis(side = 4, col=rightAxisCol)
  mtext(side = 4, line = 3, paste0('Observed # of significant ', feature_name, ' (log10)'), col = rightAxisCol)
  mtext(paste0(caller, " - ", "all datasets"), font = 3, line = -1)
  
  for(i in 1:length(y_FDR)) {
    lines(x = xseq, y = log10(y_NbrSignif[[i]]), pch=16, type="o", col = rightAxisCol, cex = 0.6)
  }
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}

