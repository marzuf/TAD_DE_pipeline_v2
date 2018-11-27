#!/usr/bin/Rscript

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: rna_rnaseqDT.Rdata
# - script0: rna_geneList.Rdata
# - script0: pipeline_geneList.Rdata
# - script3: all_meanLogFC_TAD.Rdata
# - script4: all_meanCorr_TAD.Rdata
# - script7: meanCorr_permDT.Rdata
# - script11: emp_pval_combined.Rdata
# - script21: signifTADs.Rdata
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
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script11_name <- "11_runEmpPvalCombined"
script21_name <- "21_getSignifFDR"
script_name <- "22_plotTopCombinedFDR"
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

meanTADintraCorr <- eval(parse(text = load(paste0(pipOutFold, "/", script4_name, "/", "all_meanCorr_TAD.Rdata"))))

########################################################################################
########### PREPARE THE DATA FROM EMP FDR
########################################################################################
# retrieve the signif TADs 
signifTADs_empFDR <- eval(parse(text = load(paste0(pipOutFold, "/", script21_name, "/", "signifTADs.Rdata"))))

txt <- paste0(toupper(script_name), "> Retrieved signif. TADs for emp. FDR, found: ", length(signifTADs_empFDR), " TADs\n")
printAndLog(txt, pipLogFile)

if(any(grepl("_BOUND", signifTADs_empFDR))) {
  if(useTADonly) {
    txt <- paste0(toupper(script_name), "> !!! WARNING: useTADonly==TRUE, but found non-TAD regions in empFDR signif TADs\n")
    printAndLog(txt, pipLogFile)    
  }
}
# thresh_FDRfile <- paste0(myfold, script21_name, "all_empFDR_signifThresh.Rdata")
# empFDR_thresh <- eval(parse(text = load(paste0(thresh_FDRfile))))
# empFDR_thresh <- as.numeric(empFDR_thresh["desired_empFDR_thresh"])

stopifnot(signifTADs_empFDR %in% names(meanTADlogFC))
stopifnot(signifTADs_empFDR %in% names(meanTADintraCorr))

# to rank the TADs for the empFDR, combine abs(logFC) + intraCorr
meanTAD_logFC_intraCorr <- abs(meanTADlogFC[signifTADs_empFDR]) + meanTADintraCorr[signifTADs_empFDR] 

txt <- paste0(toupper(script_name), "> Signif. TADs for emp. FDR are ranked by decreasing score=abs(meanLogFC) + meanIntraCorr\n")
printAndLog(txt, pipLogFile)


# ADDED 27.11.2018 to check using other files
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)


all_empFDR_scores <- sort(meanTAD_logFC_intraCorr, decreasing=T)

########################################################################################
########### PREPARE THE COMBINED P-VALUES
########################################################################################
### combined p-values
all_combined_pval <- eval(parse(text = load(paste0(pipOutFold, "/", script11_name, "/", "emp_pval_combined.Rdata"))))
all_combined_pval <- p.adjust(all_combined_pval, method="BH")
# sort the TADs by increasing adj. combined p-val
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


########################################################################################
########### PREPARE THE GENES TO PLOT THE logFC 
########################################################################################
samp1 <- eval(parse(text = load(paste0(setDir, "/", sample1_file))))
samp2 <- eval(parse(text = load(paste0(setDir, "/", sample2_file))))

# UPDATE SELECT THE GENES ACCORDING TO THE SETTINGS PREPARED IN 0_PREPGENEDATA
rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_rnaseqDT.Rdata"))))
initList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_geneList.Rdata"))))
geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))

txt <- paste0(toupper(script_name), "> Start with # genes: ", length(geneList), "/", length(initList), "\n")
printAndLog(txt, pipLogFile)

rnaseqDT <- rnaseqDT[names(geneList),]    
stopifnot(all(rownames(rnaseqDT) == names(geneList)))
stopifnot(!any(duplicated(names(geneList))))
stopifnot(!any(duplicated(geneList)))

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

all_tops <- c("all_combined_pval", "all_empFDR_scores")

cat("> Start drawing top-scoring TADs\n")
for(i_top in all_tops) {
  cat("... start for", i_top, "\n")
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
                                         orderBy = "logFC", cond1=plot_cond1, cond2=plot_cond2)
  }
  all_plots <- do.call(grid.arrange, c(vect_plot,  list(ncol=2, top=textGrob(paste("Top", nTopLolliPlot, toTitleCase(p_val_met), "p-val"),
                                                                             gp=gpar(fontsize=20,font=2)))))
  outFile <- paste0(curr_outFold, "/", p_val_met, "_pval_top", n_topTAD_toplot, ".", plotType)
  ggsave(filename = outFile, all_plots, width=outWidth, height = outHeight)
  cat("... written: ", outFile, "\n")
}


###################################################################### INTERSECT empFDR pvalCOMBINED
cat("> Start drawing Venn diagram intersect TADs\n")

nTopVennPlot <- min(c(length(all_empFDR_scores), length(all_combined_pval), nTopVennPlot))

# all_empFDR_nTop <- all_empFDR_scores[1:nTopVennPlot]  
# all_combined_pval_nTop <- all_combined_pval[1:nTopVennPlot] 


nTopVennPlot_empFDR <- min(c(length(all_empFDR_scores),  nTopVennPlot))
nTopVennPlot_pval <- min(c(length(all_combined_pval), nTopVennPlot))

# intersectReg <- intersect(names(all_empFDR_nTop), names(all_combined_pval_nTop))

if(nTopVennPlot_pval > 0 & nTopVennPlot_empFDR > 0) {
  
  all_empFDR_nTop <- all_empFDR_scores[1:nTopVennPlot_empFDR]
  all_combined_pval_nTop <- all_combined_pval[1:nTopVennPlot_pval]
  
  vennF <- paste0(curr_outFold, "/", "venn_combined_empFDR_intersect_top", nTopVennPlot, ".", plotType)
  do.call(plotType, list(file=vennF, height=myHeight, width=myWidth))
  draw_venn_font(vd_tit = paste0("Intersect combinedPval/empFDR top ", nTopVennPlot, " p-values/scores"),
                 my_font ="time",  my_cex = 1.2,
                 top_empFDR = names(all_empFDR_nTop),
                 top_combined = names(all_combined_pval_nTop))
  foo <- dev.off()
  
  system(paste0("rm -f VennDiagram*-*-*.log"))  
  
} else {
  cat("... find to draw: # empFDR signif = ", nTopVennPlot_empFDR, "\t # combPval signif = ", nTopVennPlot_pval,"\n")
  cat("... -> nothing to draw for Venn diagram\n")
}

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))



