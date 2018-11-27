startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script5_name <- "5_runPermutationsMedian"
script_name <- "8c_runAllDown"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
#source("run_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
source(paste0(pipScriptDir, "/", "all_down_ratio_fct.R"))

registerDoMC(ifelse(SSHFS,2, 10))

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

#allDown <- c("ratioDown", "FCdown", "prodConcord", "prodMeanConcord", "prodLogRatioNbr", "prodRatioSum") # loaded from main_settings.R

#allDown <- "meanRatioFC"

allDown_fct <- paste0("get_", allDown, "ByRegion_v2")
names(allDown_fct) <- allDown

#****************************************************************************************************************************

# INPUT DATA
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

DE_topTable <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_topTable.Rdata"))))
DE_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_geneList.Rdata"))))

stopifnot(all(DE_topTable$genes %in% names(DE_geneList)))
stopifnot(!any(duplicated(names(DE_geneList))))
### TO USE THE FUNCTIONS FOR RATIO DOWN DE_TOPTABLE$GENES => SHOULD GIVE ENTREZ ID
entrezList <- unlist(sapply(DE_topTable$genes, function(x) DE_geneList[x]))
DE_topTable$genes <- entrezList

##########################################################################################
########################################################################################## take the filtered list
##########################################################################################
pipeline_regionList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_regionList.Rdata"))))
pipeline_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))

initNrow <- nrow(DE_topTable)
DE_topTable <- DE_topTable[DE_topTable$genes %in% pipeline_geneList,]
txt <- paste0(toupper(script_name), "> Take only filtered genes: ", nrow(DE_topTable), "/", initNrow, "\n")
printAndLog(txt, pipLogFile)

gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% entrezList,]
initLen <- length(unique(gene2tadDT$region))
gene2tadDT <- gene2tadDT[gene2tadDT$region %in% pipeline_regionList,]
txt <- paste0(toupper(script_name), "> Take only filtered regions: ", length(unique(gene2tadDT$region)), "/", initLen, "\n")
printAndLog(txt, pipLogFile)

cat("... load permutation data ...\n")
permutationsDT <- eval(parse(text = load(paste0(pipOutFold, "/", script5_name, "/permutationsDT.Rdata"))))
### TAKE ONLY THE GENES FOR WHICH A LOG FC VALUE IS AVAILABLE (I.E. THE ONES USED FOR DE ANALYSIS)
initNrow <- nrow(permutationsDT)
permutationsDT <- permutationsDT[which(rownames(permutationsDT) %in% DE_topTable$genes),]
txt <- paste0(toupper(script_name), "> Take only the genes for which logFC value is available (the ones used in DE analysis): ", nrow(permutationsDT), "/", initNrow, "\n")
printAndLog(txt, pipLogFile)

all_regions <- sort(unique(as.character(permutationsDT[,2])))

### ===>> change EZH2 data, gene can map to boundary and tads, so change from stop() to warning()
if(useTADonly) {
    if(! all (regexpr("_TAD",  gene2tadDT$region[gene2tadDT$entrezID %in% rownames(permutationsDT)]) > 0 )) {
        warning("make not sense to filter TAD genes after permutations if permutations were run with genes belonging to non-TAD regions\n")
    }
    initLen <- length(all_regions)
    all_regions <- all_regions[grep("_TAD", all_regions)]
    # if want to use only TADs, make more sens if permutations were run without the TADs
    txt <- paste0(toupper(script_name), "> Take only TAD regions: ", length(all_regions), "/", initLen, "\n")
    printAndLog(txt, pipLogFile)    
    if(length(all_regions) < initLen){
        warning("make not sense to filter TAD regions after permutations if permutations were run with genes belonging to non-TAD regions\n")
    }
}

##########################################################################################
########################################################################################## CALCULATE OBSERVED RATIO DOWN
##########################################################################################

for(curr_ratio_type in allDown) {
  cat(paste0("... start calculate ", curr_ratio_type, "/TAD\n"))
  all_obs_ratio <- do.call(allDown_fct[curr_ratio_type], list(g2TADdt=gene2tadDT, DEdt=DE_topTable))
  txt <- paste0(toupper(script_name), "> Number of observed ", curr_ratio_type, " computed: ", length(all_obs_ratio), "\n")
  printAndLog(txt, pipLogFile)
  # assign(curr_ratio_type,get("b"))
  save(all_obs_ratio, file = paste0(curr_outFold, "/all_obs_", curr_ratio_type, ".Rdata"))
  cat(paste0("... written: ",  paste0(curr_outFold, "/all_obs_", curr_ratio_type, ".Rdata"), "\n"))

  cat("... start calculate permutations ", curr_ratio_type, "/TAD\n")
  ratio_permDT <- get_statFromShuffle_para(DEdt=DE_topTable, shuffData=permutationsDT, stat_fct=allDown_fct[curr_ratio_type], ncpu=nCpu, TADonly=F)
  txt <- paste0(toupper(script_name), "> Number of permutations for which ", curr_ratio_type," computed: ", ncol(ratio_permDT), "\n")
  printAndLog(txt, pipLogFile)
  txt <- paste0(toupper(script_name), "> Number of regions for which ", curr_ratio_type, " computed: ", nrow(ratio_permDT), "\n")
  printAndLog(txt, pipLogFile)
  save(ratio_permDT, file = paste0(curr_outFold, "/", curr_ratio_type, "_permDT.Rdata"))
  cat(paste0("... written: ",  paste0(curr_outFold, "/", curr_ratio_type, "_permDT.Rdata"), "\n"))
}


txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))







