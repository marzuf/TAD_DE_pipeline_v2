startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline")

script9_name <- "9_runEmpPvalMeanTADLogFC"
script10_name <- "10_runEmpPvalMeanTADCorr"
script11_name <- "11_runEmpPvalCombined"
script_name <- "16_cross_comp"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
#source("run_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar

registerDoMC(ifelse(SSHFS,2, nCpu)) # loaded from main_settings.R

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))


#*********************************************************

txt <- paste0(toupper(script_name), "> Cross-comparison for following datasets: ", paste0(comp_folders, collapse = ", "), "\n")
printAndLog(txt, pipLogFile)

#*********************************************************

# set the name of the Rdata associated with script name
script_Rdata <- setNames(c("emp_pval_meanLogFC.Rdata","emp_pval_meanCorr.Rdata", "emp_pval_combined.Rdata" ), c(script9_name, script10_name, script11_name))

#*********************************************************

gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

symbDT <- read.delim(paste0(symbolDT_file), header=TRUE, stringsAsFactors = F)
symbDT$entrezID <- as.character(symbDT$entrezID)

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 480 , 7)
myWidth <- ifelse(plotType == "png", 600, 10)


########################################### INTERSECT FOR THE MEAN TAD CORR

for(i_script in 1:length(script_Rdata)) {
  all_DS_rank <- foreach(ds = comp_folders) %dopar% {
    pval_file <- paste0(ds, "/", names(script_Rdata)[i_script], "/", script_Rdata[i_script])
    cat(paste0("... pval_file: ", pval_file, "\n"))
    pval_data <- eval(parse(text = load(pval_file)))
    rank(pval_data, ties="min")
  }
  names(all_DS_rank) <- comp_folders
  interSize <- foreach(curr_rank = 1:upToRank, .combine='c') %dopar% {
    # for all callers, take all the regions that are up to current rank
  curr_rank_DS <- lapply(all_DS_rank, function(x) names(x[x <= curr_rank]))
  length(Reduce(intersect, curr_rank_DS))
  }
  # draw intersect size ~ rank
  outFile <- paste0(curr_outFold, "/", sub(".Rdata", "", as.character(script_Rdata[i_script])), "_datasets_intersect_vs_rank_upto", upToRank, ".", plotType)
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(interSize ~ c(1:upToRank), 
       bty="l", xlab = "up to rank", ylab ="intersect size", 
       pch=16, cex=0.7, type='o',
       main=paste0(gsub("_", " " , sub(".Rdata", "", as.character(script_Rdata[i_script]))), " - datasets intersect"))
  legend("topleft", legend = comp_folders, cex=0.7, bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  if(all(interSize < 1))  next
  # take max 10 at the intersect and draw
  rankToPlot <- max(which(interSize <= 10))
  # retrieve the TADs at the intersect at this rank
  plot_rank_DS <- lapply(all_DS_rank, function(x) names(x[x <= rankToPlot]))
  TADs_to_plot <- Reduce(intersect, plot_rank_DS)
  stopifnot(length(TADs_to_plot) <= 10)
  
  rank_tad_DT <- foreach(i = TADs_to_plot, .combine = 'rbind') %dopar% {
    # retrieve the genes
    entrez_genes <- gene2tadDT$entrezID[gene2tadDT$region == i]
    symbol_genes <- unlist(sapply(entrez_genes, function(x) symbDT$symbol[symbDT$entrezID == x][1]))
    data.frame(TAD = i, entrezID = entrez_genes, symbol = symbol_genes)
  }
  outFile <- paste0(curr_outFold, "/", sub(".Rdata", "", as.character(script_Rdata[i_script])), "_datasets_intersect_upto_10_intersect_TAD.txt")
  write.table(rank_tad_DT, file = outFile, sep="\t", quote=F, col.names = T, row.names = F)
  cat(paste0("... written: ", outFile, "\n"))
  
} # end iterating over the different p-values


txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)
cat(paste0("*** DONE: ", script_name, "\n"))



