#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

#### UPDATE: do not take raw counts but fpkm data !!!

#set.seed(20180202) # this row was added 08.03.18, the files in OUTPUTFOLDER so far without set.seed

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: pipeline_regionList.Rdata
# - script0: rna_geneList.Rdata
# - script0: pipeline_geneList.Rdata
# - script0: rna_madnorm_rnaseqDT.Rdata
# - script0: rna_fpkmDT.Rdata # UPDATE 
################################################################################

################  OUTPUT
# - permutationsDT.Rdata
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
script_name <- "5v2_runPermutations"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(ifelse(SSHFS, 2, nCpu))

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

nRandom <- ifelse(SSHFS, 5, nRandomPermut)


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

initList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_geneList.Rdata"))))
geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))

txt <- paste0(toupper(script_name), "> Start with # genes: ", length(geneList), "/", length(initList), "\n")
printAndLog(txt, pipLogFile)

stopifnot(!duplicated(names(geneList)))
#*******************************************************************************

# INPUT DATA
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% as.character(geneList),]

### take only the filtered data according to initial settings
pipeline_regionList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_regionList.Rdata"))))
if(useTADonly) {
  if(any(grepl("_BOUND", pipeline_regionList))) {
    stop("! data were not prepared for \"useTADonly\" !")
  }
}
initLen <- length(unique(gene2tadDT$region))
gene2tadDT <- gene2tadDT[gene2tadDT$region %in% pipeline_regionList,]
txt <- paste0(toupper(script_name), "> Take only filtered regions: ", length(unique(gene2tadDT$region)), "/", initLen, "\n")
printAndLog(txt, pipLogFile)

if(useTADonly) {
  initLen <- length(geneList)
  TAD_gene2tadDT <- gene2tadDT[grep("_TAD",gene2tadDT$region),]
  rowsToKeep <- which(geneList %in% as.character(TAD_gene2tadDT$entrezID))
  geneList <- geneList[rowsToKeep]
  gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% as.character(geneList),]
  txt <- paste0(toupper(script_name), "> Take only genes that are within TADs: ", length(geneList), "/", initLen, "\n")
  printAndLog(txt, pipLogFile)
}

################################****************************************************************************************
####################################################### RUN PERMUTATIONS
################################****************************************************************************************

cat("... Start permutations\n")

nGenesByTAD <- setNames(as.numeric(table(gene2tadDT$region)), as.character(names(table(gene2tadDT$region))))

stopifnot(pipeline_regionList %in% names(nGenesByTAD))
stopifnot(names(nGenesByTAD) %in% pipeline_regionList )

tad_set <- names(nGenesByTAD)



shuffleData <- foreach(i_simu = 1:nRandom, .combine='cbind') %dopar% {
  
  cat("> start ", i_simu, "/", nRandom, "\n")
  
  i=0
  # while(nrow(random_g2t) < nrow(gene2tadDT)){
  while(TRUE){  # => in the while block, I should generate random_g2t
    first_gene_idx <- 1
    last_gene_idx <- NA
    cat("... ", i, "\n")
    random_tad_set <- sample(x=tad_set, size=length(tad_set), replace=FALSE)
    random_g2t <- data.frame(entrezID = character(0),
                             chromo = character(0),
                             start = numeric(0),
                             end = numeric(0),
                             region = character(0),
                             stringsAsFactors = FALSE)
    for(tad in random_tad_set)  {
      cat("tad: ", tad, "\n")
      stopifnot(tad %in% names(nGenesByTAD))
      tad_nGenes <- as.numeric(nGenesByTAD[tad])
      stopifnot(!is.na(tad_nGenes))
      last_gene_idx <- first_gene_idx + tad_nGenes - 1
      
      if(last_gene_idx > nrow(gene2tadDT)){
        txt <- paste0("i=\t", i, ":  last_gene_idx > nrow(gene2tadDT) => break\n")
        printAndLog(txt, pipLogFile)
        break
      }
      if(nrow(random_g2t) == nrow(gene2tadDT))
        break
      
      tad_geneDT <- gene2tadDT[first_gene_idx:last_gene_idx,,drop=FALSE]
      stopifnot(nrow(tad_geneDT) == tad_nGenes)
      tad_geneDT$region <- tad
      
      true_tad <- unique(gene2tadDT$region[gene2tadDT$entrezID %in% tad_geneDT$entrezID]) # in which TADs where the genes I selected ?
      
      ### OCCURS TOO OFTEN ??!
      if(length(true_tad) == 1){
        gene2tadDT[gene2tadDT$region == true_tad,]
        
        if(nGenesByTAD[true_tad] == nrow(tad_geneDT) & nGenesByTAD[true_tad] > 500){
          cat("break TRUE TAD\n")
          # stop("erreur")
          break
        }
      }
      
      random_g2t <- rbind(random_g2t, tad_geneDT)
      first_gene_idx <- last_gene_idx + 1      
      
      cat(paste0("nrow(random_g2t) =" , nrow(random_g2t), "\n"))
      cat(paste0("nrow(gene2tadDT) =" , nrow(gene2tadDT), "\n"))
      
      
    }    
    if(nrow(random_g2t) == nrow(gene2tadDT))
      break
    cat(paste0("...... nrow(random_g2t) = ", nrow(random_g2t), "\n"))
    cat(paste0("...... nrow(gene2tadDT) = ", nrow(gene2tadDT), "\n"))
    i=i+1
  } # end-while all genes reassigned !
  random_nGenesByTAD <- setNames(as.numeric(table(random_g2t$region)), as.character(names(table(random_g2t$region))))
  stopifnot(!duplicated(random_g2t$entrezID))
  stopifnot(random_nGenesByTAD == nGenesByTAD)
  stopifnot(random_g2t$entrezID %in% geneList)
  stopifnot(geneList %in% random_g2t$entrezID)
  rownames(random_g2t) <- as.character(random_g2t$entrezID)
  stopifnot(pipeline_regionList %in% random_g2t$region)
  stopifnot(random_g2t$region %in% pipeline_regionList )
  new_assign_DT <- random_g2t[geneList, "region", drop=FALSE]
  new_assign_DT
} # end-foreach nRandom permut


################################****************************************************************************************
####################################################### WRITE OUTPUT
################################****************************************************************************************
permutationsDT <- shuffleData
colnames(permutationsDT) <- paste0("region", c(1:ncol(permutationsDT)))

stopifnot(ncol(permutationsDT) == (nRandom))

txt <- paste0(toupper(script_name), "> Number of permutations: ", nRandom, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> WARNING: genes that map to duplicated entrezID are removed ! \n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> -> number of genes retained for permutations: ", nrow(permutationsDT), "/", length(geneList), "\n")
printAndLog(txt, pipLogFile)

outFile <- file.path(curr_outFold, "permutationsDT.Rdata")
save(permutationsDT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

obs_dt <- gene2tadDT

nbrSameTAD_obs_permut <- foreach(i=1:ncol(permutationsDT)) %dopar% {
  cat(i, "***\n")
  stopifnot(geneList == rownames(permutationsDT))
  
  all_permut_tads <- unique(as.character(permutationsDT[,i]))
  stopifnot(setequal(all_permut_tads, obs_dt$region))
  
  permut_g2t <- setNames(permutationsDT[,i], rownames(permutationsDT) )
  
  exact_match <- sapply(all_permut_tads, function(x) {
    
    permut_genes <- names(permut_g2t)[permut_g2t == x]
    stopifnot(length(permut_genes) > 0)
    
    obs_tads <- obs_dt$region[obs_dt$entrezID %in% permut_genes ]
    stopifnot(length(obs_tads) > 0)
    
    obs_genes <- obs_dt$entrezID[obs_dt$region %in% obs_tads]
    stopifnot(length(obs_genes) > 0)
    stopifnot(permut_genes %in% obs_genes)
    
    as.numeric(setequal(permut_genes, obs_genes))
  })
  sum(exact_match)
  exact_match
}
names(nbrSameTAD_obs_permut) <- colnames(permutationsDT)

outFile <- file.path(curr_outFold, "nbrSameTAD_obs_permut.Rdata")
save(nbrSameTAD_obs_permut, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))

