#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

# purityFilter version -> keep only the TADs that are TADs_to_keep from prep_data_purityFilter.R (cf. lines with 19.08.20 additions)


#### UPDATE: do not take raw counts but fpkm data !!!

set.seed(20180202) # this row was added 08.03.18, the files in OUTPUTFOLDER so far without set.seed => but not reproducible on multiple cores ?? cf. trial of the 16.08.2019

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
script_name <- "5pF_runPermutationsMedian"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

source(paste0(pipScriptDir, "/", "my_save_pigz.R")) # UPDATE 16.08.2019 -> to use customed fastSave save.pigz() => faster save
source(paste0(pipScriptDir, "/", "TAD_DE_utils_fasterPermut.R")) # UPDATE 16.08.2019 -> modified function for tad shuffling => faster permuts

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

nRandom <- ifelse(SSHFS, 5, nRandomPermut)

if(withExprClass)
  nClass <- permutExprClass  # number of class of expression



#####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<>>>>>>>><<< 19.08.20: add purity filter !!!
curr_exprds <- basename(pipOutFold)
curr_hicds <- basename(dirname(pipOutFold))

purityFilterFile <- file.path("//mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA_PURITYFILTER/PREP_DATA_PURITYFILTER/aran/CPE/log10/TADs_to_keep.Rdata")
stopifnot(file.exists(purityFilterFile))
all_TADs_to_keep <- get(load(purityFilterFile))

cat(head(dirname(all_TADs_to_keep)), "\n")

cat(file.path(curr_hicds, curr_exprds), "\n")

ds_TADs_to_keep <- all_TADs_to_keep[dirname(all_TADs_to_keep) == file.path(curr_hicds, curr_exprds)]
if(length(ds_TADs_to_keep) == 0) stop("--- found no TADs to keep after purity filter ! \n")
#####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<>>>>>>>><<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



# ADDED 16.11.2018 to check using other files
txt <- paste0("withExprClass\t=\t", withExprClass, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("nClass\t=\t", nClass, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)
  
#******************************************************************************* !! HARD CODED
aggregFction <- "median"
# withExprClass <- TRUE # loaded from main_settings.R !
#*******************************************************************************

################################****************************************************************************************
####################################################### PREPARE INPUT
################################****************************************************************************************

# UPDATE SELECT THE GENES ACCORDING TO THE SETTINGS PREPARED IN 0_PREPGENEDATA
if(withExprClass) {
    # rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_rnaseqDT.Rdata"))))
    rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_fpkmDT.Rdata"))))
}
initList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_geneList.Rdata"))))
geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))

txt <- paste0(toupper(script_name), "> Start with # genes: ", length(geneList), "/", length(initList), "\n")
printAndLog(txt, pipLogFile)

rnaseqDT <- rnaseqDT[names(geneList),]    
stopifnot(all(rownames(rnaseqDT) == names(geneList)))
stopifnot(!any(duplicated(names(geneList))))
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
    stopifnot(all(rownames(rnaseqDT) == names(geneList)))
    rowsToKeep <- which(geneList %in% as.character(TAD_gene2tadDT$entrezID))
    geneList <- geneList[rowsToKeep]
    rnaseqDT <- rnaseqDT[rowsToKeep,]
    gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% as.character(geneList),]
    txt <- paste0(toupper(script_name), "> Take only genes that are within TADs: ", length(geneList), "/", initLen, "\n")
    printAndLog(txt, pipLogFile)
}

#####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<>>>>>>>><<< 19.08.20: add purity filter !!!


init_nTads <- length(unique(gene2tadDT$region))
init_nGenes <- length(unique(gene2tadDT$entrezID))
stopifnot(init_nGenes == length(geneList))

pF_gene2tadDT <- gene2tadDT[gene2tadDT$region %in% basename(ds_TADs_to_keep),]
stopifnot(nrow(pF_gene2tadDT) > 0 )

stopifnot(all(rownames(rnaseqDT) == names(geneList)))
rowsToKeep <- which(geneList %in% as.character(pF_gene2tadDT$entrezID))
geneList <- geneList[rowsToKeep]
rnaseqDT <- rnaseqDT[rowsToKeep,]
stopifnot(as.character(geneList) %in% gene2tadDT$entrezID)
gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% as.character(geneList),]

stopifnot(gene2tadDT$region %in% basename(ds_TADs_to_keep))
stopifnot(length(gene2tadDT$entrezID) == length(geneList))

pF_nTads <- length(unique(gene2tadDT$region))
pF_nGenes <- length(unique(gene2tadDT$entrezID))

txt <- paste0(toupper(script_name), "> Take only TADs that pass purity filter: ", pF_nTads, "/", init_nTads, "\n")
printAndLog(txt, pipLogFile)

txt <- paste0(toupper(script_name), "> Take only genes that pass purity filter: ", pF_nGenes, "/", init_nGenes, "\n")
printAndLog(txt, pipLogFile)



#####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<>>>>>>>><<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


################################****************************************************************************************



################################****************************************************************************************
####################################################### RUN PERMUTATIONS
################################****************************************************************************************

cat("... Start permutations\n")
if(withExprClass) {
    shuffleData <- get_multiShuffledPositions_vFunct(g2TADdt=gene2tadDT, RNAdt=rnaseqDT, 
                                              geneIDlist=geneList, nClass = nClass, TADonly=F, nSimu=nRandom, 
                                              withExprClass=withExprClass, nCpu=nCpu, aggregFun = aggregFction)
} else {
    shuffleData <- get_multiShuffledPositions_vFunct(g2TADdt=gene2tadDT, RNAdt=NULL, 
                                              geneIDlist=geneList, nClass = NULL, TADonly=F, nSimu=nRandom,
                                              withExprClass=withExprClass, nCpu=nCpu, aggregFun = aggregFction)
}

################################****************************************************************************************
####################################################### WRITE OUTPUT
################################****************************************************************************************
permutationsDT <- shuffleData
rownames(permutationsDT) <- permutationsDT[,1]
permutationsDT <- permutationsDT[,-1]

txt <- paste0(toupper(script_name), "> Number of permutations: ", nRandom, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> With expression classes: ", as.character(withExprClass), "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> WARNING: genes that map to duplicated entrezID are removed ! \n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> -> number of genes retained for permutations: ", nrow(permutationsDT), "/", length(geneList), "\n")
printAndLog(txt, pipLogFile)

#save(permutationsDT, file = paste0(curr_outFold, "/permutationsDT.Rdata"))
# update 16.08.2019 => faster save version
my_save.pigz(permutationsDT, pigz_exec_path = pigz_exec_path, file = paste0(curr_outFold, "/permutationsDT.Rdata") )
cat(paste0("... written: ", paste0(curr_outFold, "/permutationsDT.Rdata"), "\n"))



txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))

cat("stopifnot ", ncol(permutationsDT), " == ", (nRandom), "\n")
stopifnot(ncol(permutationsDT) == (nRandom))
#####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<>>>>>>>><<< 19.08.20: add purity filter !!!
cat("stopifnot ", nrow(permutationsDT), " == ", (pF_nGenes), "\n")
stopifnot(nrow(permutationsDT) == (pF_nGenes))


cat("dim(permutationsDT)\n")
cat(dim(permutationsDT))
cat("\n")

cat("... using faster save and permut\n")

