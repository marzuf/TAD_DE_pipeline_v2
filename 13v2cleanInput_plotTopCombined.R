#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: rna_rnaseqDT.Rdata
# - script0: rna_geneList.Rdata
# - script0: pipeline_geneList.Rdata
# - script2v2: all_wilcoxStat_TAD.Rdata
# - script7v2: meanCorr_permDT.Rdata
# - script11v2: emp_pval_combined.Rdata
################################################################################

################  OUTPUT
# - emp_pval_meanCorr.Rdata + plots
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
script2v2_name <- "2v2_runWilcoxonTAD"
script11_name <- "11v2_runEmpPvalCombined"
script_name <- "13v2_plotTopCombined"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
#source("run_settings.R")
source(settingF)
# settingF = "SETTING_FILES_NOVOOM/run_settings_GSE71119_dediffSM_MFSM.R"
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(grid, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(tools, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 480, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

geneRankOnPlot <- TRUE


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

# INPUT DATA
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

entrez2symbDT <- read.delim(entrezDT_file, header=T, stringsAsFactors=F)
entrez2symbDT <- entrez2symbDT[,c("entrezID", "symbol")]
colnames(entrez2symbDT) <- c("entrezID", "geneName")
entrez2symbDT$entrezID <- as.character(entrez2symbDT$entrezID)

meanTADlogFC <- eval(parse(text = load(paste0(pipOutFold, "/", script3_name, "/", "all_meanLogFC_TAD.Rdata"))))

### wilcoxon p-values
all_wilcox_ttest <- eval(parse(text = load(paste0(pipOutFold, "/", script2_name, "/", "wilcox_ttest_meanTAD_qq.Rdata"))))
all_wilcox_pval <- unlist(sapply(all_wilcox_ttest, function(x) x[["wilcox_pval"]]))
all_wilcox_pval <- sort(all_wilcox_pval)
initLen <- length(all_wilcox_pval)
# plot only the TADs
all_wilcox_pval <- all_wilcox_pval[grep("_TAD", names(all_wilcox_pval))]
txt <- paste0(toupper(script_name), "> Plot only TAD regions, Wilcoxon p-values: ", length(all_wilcox_pval), "/", initLen, "\n")
printAndLog(txt, pipLogFile)    
if(useTADonly & length(all_wilcox_pval) < initLen) {
    txt <- paste0(toupper(script_name), "> !!! WARNING: useTADonly==TRUE, but found non-TAD regions in Wilcoxon p-values\n")
    printAndLog(txt, pipLogFile)    
}

### combined p-values
all_combined_pval <- eval(parse(text = load(paste0(pipOutFold, "/", script11_name, "/", "emp_pval_combined.Rdata"))))
all_combined_pval <- sort(all_combined_pval)
initLen <- length(all_combined_pval)
# plot only the TADs
all_combined_pval <- all_combined_pval[grep("_TAD", names(all_combined_pval))]
txt <- paste0(toupper(script_name), "> Plot only TAD regions, Combined p-values: ", length(all_combined_pval), "/", initLen, "\n")
printAndLog(txt, pipLogFile)    
if(useTADonly & length(all_combined_pval) < initLen) {
    txt <- paste0(toupper(script_name), "> !!! WARNING: useTADonly==TRUE, but found non-TAD regions in Combined p-values\n")
    printAndLog(txt, pipLogFile)    
}
samp1 <- eval(parse(text = load(paste0(setDir, "/", sample1_file))))
samp2 <- eval(parse(text = load(paste0(setDir, "/", sample2_file))))

# UPDATE SELECT THE GENES ACCORDING TO THE SETTINGS PREPARED IN 0_PREPGENEDATA
rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_rnaseqDT.Rdata"))))
initList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_geneList.Rdata"))))
geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))

txt <- paste0(toupper(script_name), "> Start with # genes: ", length(geneList), "/", length(initList), "\n")
printAndLog(txt, pipLogFile)

## for some datasets, the count matrix contains negative values => run without voom and no cpm
#if(any(unlist(sapply(runWithoutVoom, function(x) regexpr(x, settingF) > 0)))) {
#  applyVoomAndCPM <- FALSE
#} else{
#  applyVoomAndCPM <- TRUE
#}
#cat(paste0("> applyVoomAndCPM: ", as.character(applyVoomAndCPM), "\n"))

rnaseqDT <- rnaseqDT[names(geneList),]    
stopifnot(all(rownames(rnaseqDT) == names(geneList)))
stopifnot(!any(duplicated(names(geneList))))

DE_topTable <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_topTable.Rdata"))))
stopifnot(!any(duplicated(names(geneList))))
DE_topTable <- DE_topTable[DE_topTable$genes %in% names(geneList),]
stopifnot(nrow(DE_topTable) > 0)

gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% geneList,]

DE_topTable$genes <- unlist(sapply(DE_topTable$genes, function(x) geneList[x]))
rownames(DE_topTable) <- NULL
 
# ! duplicated row.names are not allowed !
dupEntrez <- unique(geneList[duplicated(geneList)])
geneList <- geneList[! geneList %in% dupEntrez]
stopifnot(!any(duplicated(geneList)))

DE_topTable <- DE_topTable[!DE_topTable$genes %in% dupEntrez,]
stopifnot(!any(duplicated(DE_topTable$genes)))

initNrow <- nrow(rnaseqDT)
rnaseqDT <- rnaseqDT[which(rownames(rnaseqDT) %in% names(geneList)),]
txt <- paste0(toupper(script_name), "> Discard duplicated symbol, retain: ", nrow(rnaseqDT), "/", initNrow , " genes\n")
printAndLog(txt, pipLogFile)

stopifnot(all(rownames(rnaseqDT) == names(geneList)))
stopifnot(is.numeric(rnaseqDT[1,1]))
#if(applyVoomAndCPM) {
if(inputDataType == "raw" | inputDataType == "RSEM") {
  log2_rnaseqDT <- log2(rnaseqDT + 0.0001) 
} else{
  log2_rnaseqDT <- rnaseqDT
}
meanExpr <- rowMeans(log2_rnaseqDT, na.rm=T)
names(meanExpr) <- unlist(sapply(names(meanExpr), function(x) geneList[x]))

################################****************************************************************************************
####################################################### DO THE PLOTS 
################################****************************************************************************************
###################################################################### PLOT LOLLI TAD

library(rlang) # used internally for the package used to draw the lolliplots
as_dictionary <- as_data_pronoun


# retrieve which condition is cond1, i.e. the one that is more expressed when logFC is positive
geneHighestLogFC <- names(geneList[geneList == DE_topTable$genes[which.max(DE_topTable$logFC)] ])

samp1_vect <- rnaseqDT[geneHighestLogFC,samp1, drop=F]
stopifnot(dim(samp1_vect) == c(1,length(samp1)))
samp2_vect <- rnaseqDT[geneHighestLogFC,samp2, drop=F]
stopifnot(dim(samp2_vect) == c(1,length(samp2)))

plot_cond1 <- ifelse(as.numeric(rowMeans(samp1_vect, na.rm = T)) > as.numeric(rowMeans(samp2_vect, na.rm = T)), cond1, cond2)
plot_cond2 <- ifelse(as.numeric(rowMeans(samp1_vect, na.rm = T)) > as.numeric(rowMeans(samp2_vect, na.rm = T)), cond2, cond1)

n_topTAD_toplot <- nTopLolliPlot

outWidth <- 20
outHeight <- min(c(7 * n_topTAD_toplot/2, 49))

all_tops <- c("all_combined_pval")

for(i_top in all_tops) {
  tads_to_plot <- names(eval(parse(text = i_top)))
  p_val_met <- gsub(".+_(.+)_pval", "\\1", i_top)
  vect_plot <- list()
  for(i_plot in 1:n_topTAD_toplot)  {
    tad_to_plot <- tads_to_plot[i_plot]
    vect_plot[[i_plot]] <- plot_lolliTAD(TAD_to_plot = tad_to_plot,
                                         meanExpr_vect = meanExpr, 
                                         DE_table = DE_topTable,
                                         g2t_table = gene2tadDT, 
                                         id2name_table=entrez2symbDT, 
                                         geneList = geneList,
                                         textLeft =  meanTADlogFC[tad_to_plot] > 0,
                                         # orderBy = "wilcoxStat", cond1=plot_cond1, cond2=plot_cond2)
                                         orderBy = "wilcoxStat", cond1=plot_cond1, cond2=plot_cond2,
                                         labelWithRank = geneRankOnPlot)
  }
  all_plots <- do.call(grid.arrange, c(vect_plot,  list(ncol=2, top=textGrob(paste("Top", nTopLolliPlot, toTitleCase(p_val_met), "p-val"),
                                                                             gp=gpar(fontsize=20,font=2)))))
  outFile <- paste0(curr_outFold, "/", p_val_met, "_pval_top", n_topTAD_toplot, ".", plotType)
  ggsave(filename = outFile, all_plots, width=outWidth, height = outHeight)
  cat("... written: ", outFile, "\n")
}




txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))



