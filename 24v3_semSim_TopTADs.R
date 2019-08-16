#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()


# v2: randomize only once for each TAD size

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: pipeline_geneList.Rdata
# - script11: emp_pval_combined.Rdata
################################################################################

################  OUTPUT
# - topTADs_semSim.Rdata
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
script11_name <- "11_runEmpPvalCombined"

script_name <- "24v3_semSim_TopTADs"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
#source("run_settings.R")
source(settingF)
# settingF = "SETTING_FILES_NOVOOM/run_settings_GSE71119_dediffSM_MFSM.R"
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(GOSemSim, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 

registerDoMC(ifelse(SSHFS, 2, nCpu)) # from main_settings.R
# registerDoMC(2)

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

# ADDED 27.11.2018 to check using other files
txt <- paste0("inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)


plotType <- "svg"
myHeight <- ifelse(plotType == "png", 480, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

### RETRIEVE FROM MAIN_SETTINGS.R
stopifnot(exists("nTopTADs_semSim"))
stopifnot(exists("nRandomSemSim"))
stopifnot(exists("semSimMetric"))
stopifnot(exists("combineSemSimMethod"))
stopifnot(semSimMetric %in%  c("Resnik", "Lin", "Rel", "Jiang", "Wang" ))
if(semSimMetric== "Wang") stop("-- for Wang: use 24v2Wang -> will be faster\n")
stopifnot(combineSemSimMethod %in%  c("max", "avg", "rcmax", "BMA"))
stopifnot(is.numeric(nTopTADs_semSim))
stopifnot(is.numeric(nRandomSemSim))

cat(paste0("... found setting: semSimMetric = ", semSimMetric, "\n"))
cat(paste0("... found setting: combineSemSimMethod = ", combineSemSimMethod, "\n"))
cat(paste0("... found setting: nTopTADs_semSim = ", nTopTADs_semSim, "\n"))
cat(paste0("... found setting: nRandomSemSim = ", nRandomSemSim, "\n"))

# #****************************** FOR DEBUG
# semSimMetric <- "Wang"
# combineSemSimMethod <- "BMA"
# nTopTADs_semSim <- 10
# nRandomSemSim <- 100
# tad_genes <- c("835", "5261","241")
#******************************

################################****************************************************************************************
####################################################### PREPARE INPUT
################################****************************************************************************************
# if(semSimMetric == "Wang") {
#   hsGO <- godata('org.Hs.eg.db', ont="BP", computeIC = semSimMetric!="Wang")
# } else {
#   hsGO <- godata('org.Hs.eg.db', ont="BP")
# }
# hsGO <- godata('org.Hs.eg.db', ont="BP", computeIC = semSimMetric!="Wang")
hsGO <- godata('org.Hs.eg.db', ont="BP", computeIC = semSimMetric!="Wang")

# INPUT DATA
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))


curr_g2t_DT <- gene2tadDT[gene2tadDT$entrezID %in% geneList,]

pval_combined <- eval(parse(text = load(paste0(pipOutFold, "/", script11_name, "/", "emp_pval_combined.Rdata"))))
adj_pval_combined <- sort(p.adjust(pval_combined, method="BH"))

# retain the top nTopTADs TADs based on the combined pval
# use rank because of the ties !

topTADs <- names(rank(adj_pval_combined, ties="min")[rank(adj_pval_combined, ties="min") <= nTopTADs_semSim])

topTADs_genes <- lapply(topTADs, function(x) as.character(curr_g2t_DT$entrezID[curr_g2t_DT$region == x]))
names(topTADs_genes) <- topTADs

### GET ALL LENGTH OF THE TOPTADS SO THAT I RUN PERMUT SEMSIM ONLY ONCE FOR EACH SIZE

topTADs_size <- unique(c(sapply(topTADs_genes, length)))

cat(paste0("... found unique sizes of topTADs: ", length(topTADs_size), "\n"))


all_sizes_random_semSim <- lapply(topTADs_size, function(x) {
  cat("... start random semSim for size:", x, "genes \n")
  # need to use while because of potential NA returned by mgeneSim
  all_random_semSim <- c()
  while(length(all_random_semSim) < nRandomSemSim) {
    cat("...... random semSim:", length(all_random_semSim), "/", nRandomSemSim, "\n")
    stillMissing <- nRandomSemSim - length(all_random_semSim)
    
    missing_random_semSim <- foreach(curr_perm = seq_len(stillMissing), .combine='c') %do% { 
      
      cat(paste0("... for topTADs_size: ", x, " - curr_perm = ", curr_perm, "/", stillMissing, "\n"))
      
      # select random genes
      # randomGenes <- sample(x=curr_g2t_DT$entrezID[! curr_g2t_DT$entrezID %in% tad_genes], size = length(tad_genes))
      randomGenes <- sample(x=as.character(curr_g2t_DT$entrezID), size = x)
      
      random_semSim <- mgeneSim(genes=randomGenes,
                                semData=hsGO,
                                combine=combineSemSimMethod,
                                measure=semSimMetric,
                                verbose=FALSE)
      cat("...... mgeneSim -> done\n")
      
      if(length(random_semSim) == 1) {
        random_semSim <- random_semSim
      } else{
        random_semSim <- mean(random_semSim[lower.tri(random_semSim, diag=F)], na.rm=T)
      }
      random_semSim
    } 
    # but random_semSim can return NaN if entrez not in the database
    missing_random_semSim <- missing_random_semSim[!is.na(missing_random_semSim)]
    all_random_semSim <- c(all_random_semSim, missing_random_semSim)
  }
  stopifnot(length(all_random_semSim) == nRandomSemSim)
  all_random_semSim
})

names(all_sizes_random_semSim) <- topTADs_size

topTADs_semSim <- lapply(names(topTADs_genes), function(x) {
  # topTADs_semSim <- foreach(x=names(topTADs_genes)) %do% {

    
  cat("... compute TRUE semantic similarity for TAD:", x, "\n")
  tad_genes <- topTADs_genes[[x]]
  tad_semSim <-  mgeneSim(genes=tad_genes,
             semData=hsGO,
             combine=combineSemSimMethod,
             measure=semSimMetric,
             verbose=FALSE)
  cat("...... mgeneSim -> done\n")
  
  if(length(tad_semSim) == 1) {
    tad_semSim <- tad_semSim
  } else{
    tad_semSim <- mean(tad_semSim[lower.tri(tad_semSim, diag=F)], na.rm=T)
  }
  
  stopifnot(length(tad_genes) %in% names(all_sizes_random_semSim))
  curr_size_random_semSim <- all_sizes_random_semSim[[paste0(length(tad_genes))]]
  
  emp_p_val_semSimTAD <- (sum(tad_semSim <=  curr_size_random_semSim) + 1) / (length(curr_size_random_semSim) + 1)
  
  list(
    TAD_semSim = tad_semSim,
    TAD_semSim_enpPval = emp_p_val_semSimTAD,
    settings = c(semSimMetric=semSimMetric, 
                 combineSemSimMethod=combineSemSimMethod, 
                 nTopTADs_semSim=nTopTADs_semSim, 
                 nRandomSemSim=nRandomSemSim)
    )
# }
})
names(topTADs_semSim) <- names(topTADs_genes)

outFile <- paste0(curr_outFold, "/", "topTADs_semSim", semSimMetric, ".Rdata")
save(topTADs_semSim, file = outFile)
cat(paste0("... written: ", outFile, "\n"))



############################################################################################################################
############################################################################################################################
############################################################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))



