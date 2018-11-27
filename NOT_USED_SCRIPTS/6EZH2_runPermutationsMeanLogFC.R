startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script5_name <- "5_runPermutationsMedian"
script_name <- "6_runPermutationsMeanLogFC"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
#source("run_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(ifelse(SSHFS, 2, nCpu)) # from main_settings.R

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

#****************************************************************************************************************************

# INPUT DATA
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

DE_topTable <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_topTable.Rdata"))))
DE_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_geneList.Rdata"))))

stopifnot(all(DE_topTable$genes %in% names(DE_geneList)))
stopifnot(!any(duplicated(names(DE_geneList))))

entrezList <- unlist(sapply(DE_topTable$genes, function(x) DE_geneList[x]))
stopifnot(length(entrezList) == length(DE_topTable$genes))

# replace the gene symbol rownames by ensemblID rownames
logFC_DT <- data.frame(entrezID =  entrezList,
                       logFC = DE_topTable$logFC, stringsAsFactors = F)
rownames(logFC_DT) <- NULL

gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% entrezList,]

cat("... load permutation data ...\n")
permutationsDT <- eval(parse(text = load(paste0(pipOutFold, "/", script5_name, "/permutationsDT.Rdata"))))
# stopifnot(all(entrezList %in% rownames(permutationsDT) ))
# stopifnot(all(rownames(permutationsDT) %in%  entrezList))
# permutationsDT[1:5,1:5]


pipeline_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))
if(!setequal(pipeline_geneList, rownames(permutationsDT))) {
  txtWarningGene <- paste0(toupper(script_name), "> Not the same set of genes in permutDT and pipeline_geneList\n")
  printAndLog(txtWarningGene, pipLogFile)    
} else{
  txtWarningGene <- ""
}

pipeline_regionList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_regionList.Rdata"))))
if(!setequal(pipeline_geneList, permutationsDT[,2])) {
  txtWarningRegion <- paste0(toupper(script_name), "> Not the same set of regions in permutDT and pipeline_regionList\n")
  printAndLog(txtWarningRegion, pipLogFile)    
} else {
  txtWarningRegion <- ""
}


### TAKE ONLY THE GENES FOR WHICH A LOG FC VALUE IS AVAILABLE (I.E. THE ONES USED FOR DE ANALYSIS)
initNrow <- nrow(permutationsDT)
permutationsDT <- permutationsDT[which(rownames(permutationsDT) %in% logFC_DT$entrezID),]
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

#****************************************************************************************************************************
#****************************************************************************************************************************
# COMPUTE MEAN LOG FC BY TAD FOR THE PERMUTATIONS
#****************************************************************************************************************************
#****************************************************************************************************************************

### REGIONS ARE STORED IN ROWNAMES OF PERMUTdt

cat("... start logFC permutDT \n")

logFC_permDT <- foreach(i_col = 1:ncol(permutationsDT), .combine='cbind') %dopar% {
  cat(paste0("...TAD logFC for permutation: ", i_col, "/", ncol(permutationsDT), "\n"))  
  
  g2t_permDT <- data.frame(entrezID = rownames(permutationsDT), 
                           region = permutationsDT[,i_col], stringsAsFactors = F)
  g2t_permDT$entrezID <- as.character(g2t_permDT$entrezID)
  g2t_permDT$region <-  as.character(g2t_permDT$region)

  unlist(sapply(unique(all_regions), function(x) {
    reg_genes <- g2t_permDT$entrezID[which(g2t_permDT$region == x)]
    # head(reg_genes)
    mean(logFC_DT$logFC[logFC_DT$entrezID %in% reg_genes], na.rm=T)
  }))
}

cat("... end logFC permutDT \n")

meanLogFC_permDT <- as.data.frame(logFC_permDT)
stopifnot(ncol(meanLogFC_permDT) == ncol(permutationsDT))  
colnames(meanLogFC_permDT) <- paste0("permutation",  c(1:ncol(permutationsDT)))
stopifnot(nrow(meanLogFC_permDT) == length(all_regions))
rownames(meanLogFC_permDT) <- all_regions

txt <- paste0(toupper(script_name), "> Number of permutations for which mean logFC computed: ", ncol(meanLogFC_permDT), "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> Number of regions for which mean logFC computed: ", nrow(meanLogFC_permDT), "\n")
printAndLog(txt, pipLogFile)

save(meanLogFC_permDT, file= paste0(curr_outFold, "/meanLogFC_permDT.Rdata"))
cat(paste0("... written: ", curr_outFold, "/meanLogFC_permDT.Rdata", "\n"))

cat(paste0(txtWarningGene))
cat(paste0(txtWarningRegion))


txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))




