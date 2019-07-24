#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: pipeline_geneList.Rdata
# - script1: DE_topTable.Rdata
# - script11: emp_pval_combined.Rdata
################################################################################

################  OUTPUT
# - [...].Rdata + plots
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

script_name <- "23topTable_funcEnrich_TopTableGenes_TADsGenes"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
#source("run_settings.R")
source(settingF)
# settingF = "SETTING_FILES_NOVOOM/run_settings_GSE71119_dediffSM_MFSM.R"
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(clusterProfiler, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(org.Hs.eg.db, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 

# library(foreach)
# library(doMC)

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

# ADDED 27.11.2018 to check using other files
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

# number of topTADs for GO
stopifnot(exists("funcTopTable_nPerm_gseGO"))
stopifnot(exists("funcTopTable_p_adj_DE_genes"))
stopifnot(exists("funcTopTable_nShowPlot"))

cat(paste0("... found setting: funcTopTable_nPerm_gseGO = ", funcTopTable_nPerm_gseGO, "\n"))
cat(paste0("... found setting: funcTopTable_p_adj_DE_genes = ", funcTopTable_p_adj_DE_genes, "\n"))
cat(paste0("... found setting: funcTopTable_nShowPlot = ", funcTopTable_nShowPlot, "\n"))


current_settings <- c(
  funcTopTable_nPerm_gseGO=funcTopTable_nPerm_gseGO,
  funcTopTable_p_adj_DE_genes=funcTopTable_p_adj_DE_genes,
  funcTopTable_nShowPlot=funcTopTable_nShowPlot
)

outFile <- paste0(curr_outFold, "/", "current_settings.Rdata")
save(current_settings, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

################################****************************************************************************************
####################################################### PREPARE INPUT
################################****************************************************************************************

# INPUT DATA
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))

DE_topTable <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_topTable.Rdata"))))
stopifnot(DE_topTable$genes == rownames(DE_topTable))
stopifnot(names(geneList) %in% DE_topTable$genes)
curr_topTable_DT <- DE_topTable[DE_topTable$genes %in% names(geneList),]
stopifnot(curr_topTable_DT$genes %in% names(geneList)) # not true
stopifnot(nrow(curr_topTable_DT) > 0)
rm(DE_topTable)
curr_topTable_DT$entrezID <- unlist(sapply(curr_topTable_DT$genes, function(x) geneList[names(geneList) == x]))
stopifnot(!any(is.na(curr_topTable_DT$entrezID)))

curr_g2t_DT <- gene2tadDT[gene2tadDT$entrezID %in% geneList,]

pval_combined <- eval(parse(text = load(paste0(pipOutFold, "/", script11_name, "/", "emp_pval_combined.Rdata"))))
adj_pval_combined <- sort(p.adjust(pval_combined, method="BH"))

# retain the top nTopTADs TADs based on the combined pval
# use rank because of the ties !

# select the same number of genes, taking the same # based on limma pval
curr_topTable_DT <- curr_topTable_DT[order(curr_topTable_DT$adj.P.Val, curr_topTable_DT$P.Value),]

# main change here -> select the genes based on their adjusted p-val
topTopTable_genes_symb <- curr_topTable_DT$genes[curr_topTable_DT$adj.P.Val <= funcTopTable_p_adj_DE_genes]
topTopTable_genes <- unlist(sapply(topTopTable_genes_symb, function(x) as.character(geneList[x])))

outFile <- paste0(curr_outFold, "/", "topTopTable_genes.Rdata")
save(topTopTable_genes, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

customBarplotTitle <- paste0("Top ", funcTopTable_nShowPlot, " most signif. enriched GO")
customBarplotSubtit <- paste0("# DE genes tot.: ", length(topTopTable_genes))


################################################################################################################
######################################################## GO Enrichment Analysis of a gene set (enrichMSigDB)
################################################################################################################


################################################################################################################
######################################################## Using MSigDB gene set collections
################################################################################################################

# gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
# c5 <- read.gmt(gmtfile)
gmtfile <- file.path(pipScriptDir, "h.all.v6.1.entrez.gmt")
h_msigdb <- read.gmt(gmtfile)

topTopTable_genes_enrichMSigDB <- enricher(topTopTable_genes, 
                                       # TERM2GENE=c5,
                                       TERM2GENE=h_msigdb,
                                       pvalueCutoff = 1, 
                                       pAdjustMethod = "BH", 
                                       minGSSize = 1, 
                                       maxGSSize = 500, 
                                       qvalueCutoff = 1
)


if(!is.null(topTopTable_genes_enrichMSigDB)) {
  
  outFile <- paste0(curr_outFold, "/", "topTopTable_genes_enrichMSigDB_dotplot.", plotType)
  # do.call(plotType, list(file=outFile, height=myHeight, width=myWidth*2))
  # dotplot(topTopTable_genes_enrichMSigDB, showCategory=funcTADs_nShowPlot)
  # foo <- dev.off()
  p <- dotplot(topTopTable_genes_enrichMSigDB, showCategory=funcTADs_nShowPlot)
  ggsave(p, filename = outFile, height=myHeight, width=myWidth*2)
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- paste0(curr_outFold, "/", "topTopTable_genes_enrichMSigDB_barplot.", plotType)
  # do.call(plotType, list(file=outFile, height=myHeight, width=myWidth*2))
  # barplot(topTopTable_genes_enrichMSigDB, showCategory=funcTADs_nShowPlot, order=T)
  # foo <- dev.off()
  p <- barplot(topTopTable_genes_enrichMSigDB, showCategory=funcTADs_nShowPlot, order=T)
  ggsave(p, filename = outFile, height=myHeight, width=myWidth*2)
  cat(paste0("... written: ", outFile, "\n"))
  
  topTopTable_genes_enrichMSigDB_resultDT <- topTopTable_genes_enrichMSigDB@result
  topTopTable_genes_enrichMSigDB_resultDT <- topTopTable_genes_enrichMSigDB_resultDT[order(topTopTable_genes_enrichMSigDB_resultDT$p.adjust, topTopTable_genes_enrichMSigDB_resultDT$pvalue),]
  
  # topTopTable_genes_enrichMSigDB <- topTopTable_genes_enrichMSigDB_resultDT[1:min(nTopGO, nrow(topTopTable_genes_enrichMSigDB_resultDT)),]
  # outFile <- paste0(curr_outFold, "/", "topTopTable_genes_enrichMSigDB.Rdata")
  # save(topTopTable_genes_enrichMSigDB, file = outFile)
  # cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- paste0(curr_outFold, "/", "topTopTable_genes_enrichMSigDB_custombarplot_padj_withNbr.", plotType)
  try(dev.off())
  do.call(plotType, list(file=outFile, height=myHeight*1.5, width=myWidth*2))
  barplot_funcEnrich_results(resultDT = topTopTable_genes_enrichMSigDB_resultDT, 
                             nTop = funcTADs_nShowPlot, 
                             xVar = "p.adjust", 
                             g2t_DT = gene2tadDT, 
                             DE_genes = topTopTable_genes, 
                             plotHoriz=TRUE,
                             barcol="slategray")
  title(main = paste0(customBarplotTitle, " - enrichMSigDB"))
  mtext(text = customBarplotSubtit, side = 3, font = 3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}else{
  topTopTable_genes_enrichMSigDB_resultDT <- NULL
  txt <- paste0(toupper(script_name), "> ", "topTopTable_genes_enrichMSigDB is NULL\n")
  printAndLog(txt, pipLogFile)
}

outFile <- paste0(curr_outFold, "/", "topTopTable_genes_enrichMSigDB_resultDT.Rdata")
save(topTopTable_genes_enrichMSigDB_resultDT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))



################################################################################################################
######################################################## GO Enrichment Analysis of a gene set (enrichGO)
################################################################################################################

# use an over-representation test rather than an enrichment, because here I have already kind of thresholded list
# of genes as I take only genes from the topTADs (I don't care about their scores)

### FOR THE topTopTable_genes
topTopTable_genes_scoreDT <- curr_topTable_DT[curr_topTable_DT$genes %in% topTopTable_genes_symb, c("genes","adj.P.Val", "entrezID")]

if("keytype" %in% formalArgs(enrichGO)) {
topTopTable_genes_enrichGO <- enrichGO(gene = topTopTable_genes_scoreDT$entrezID,
                                       OrgDb = org.Hs.eg.db,
                                       # change default min gene set size (10)
                                       minGSSize = 1,
                                       keytype = "ENTREZID",
                                       pvalueCutoff = 1,
                                       qvalueCutoff  = 1,
                                       pAdjustMethod = "BH",
                                       ont="ALL")
} else {
  topTopTable_genes_enrichGO <- enrichGO(gene = topTopTable_genes_scoreDT$entrezID,
                                         OrgDb = org.Hs.eg.db,
                                         # change default min gene set size (10)
                                         minGSSize = 1,
                                         keyType = "ENTREZID",
                                         pvalueCutoff = 1,
                                         qvalueCutoff  = 1,
                                         pAdjustMethod = "BH",
                                         ont="ALL")
}

if(!is.null(topTopTable_genes_enrichGO)) {
  
  outFile <- paste0(curr_outFold, "/", "topTopTable_genes_enrichGO_dotplot.", plotType)
  # do.call(plotType, list(file=outFile, height=myHeight, width=myWidth*2))
  # # try(dotplot(topTopTable_genes_enrichGO, showCategory=funcTopTable_nShowPlot))
  # dotplot(topTopTable_genes_enrichGO, showCategory=funcTopTable_nShowPlot)
  # foo <- dev.off()
  p <- dotplot(topTopTable_genes_enrichGO, showCategory=funcTopTable_nShowPlot)
  ggsave(p, filename = outFile, height=myHeight, width=myWidth*2)
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- paste0(curr_outFold, "/", "topTopTable_genes_enrichGO_barplot.", plotType)
  # do.call(plotType, list(file=outFile, height=myHeight, width=myWidth*2))
  # # try(barplot(topTopTable_genes_enrichGO, showCategory=funcTopTable_nShowPlot, order=T))
  # barplot(topTopTable_genes_enrichGO, showCategory=funcTopTable_nShowPlot, order=T)
  # foo <- dev.off()
  p <- barplot(topTopTable_genes_enrichGO, showCategory=funcTopTable_nShowPlot, order=T)
  ggsave(p, filename = outFile, height=myHeight, width=myWidth*2)
  cat(paste0("... written: ", outFile, "\n"))
  
  topTopTable_genes_enrichGO_resultDT <- topTopTable_genes_enrichGO@result
  topTopTable_genes_enrichGO_resultDT <- topTopTable_genes_enrichGO_resultDT[order(topTopTable_genes_enrichGO_resultDT$p.adjust, topTopTable_genes_enrichGO_resultDT$pvalue),]
  
  
  outFile <- paste0(curr_outFold, "/", "topTopTable_genes_enrichGO_custombarplot_padj_withNbr.", plotType)
  try(dev.off())
  do.call(plotType, list(file=outFile, height=myHeight*1.5, width=myWidth*2))
  barplot_funcEnrich_results(resultDT = topTopTable_genes_enrichGO_resultDT, 
                             nTop = funcTADs_nShowPlot, 
                             xVar = "p.adjust", 
                             g2t_DT = gene2tadDT, 
                             DE_genes = topTopTable_genes, 
                             plotHoriz=TRUE,
                             barcol="slategray")
  title(main = paste0(customBarplotTitle, " - enrichGO"))
  mtext(text = customBarplotSubtit, side = 3, font = 3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  # topTopTable_genes_enrichGO_topResultDT <- topTopTable_genes_enrichGO_resultDT[1:min(nTopGO, nrow(topTopTable_genes_enrichGO_resultDT)),]
  # outFile <- paste0(curr_outFold, "/", "topTopTable_genes_enrichGO_topResultDT.Rdata")
  # save(topTopTable_genes_enrichGO_topResultDT, file = outFile)
  # cat(paste0("... written: ", outFile, "\n"))
  

} else{
  topTopTable_genes_enrichGO_resultDT <- NULL
  txt <- paste0(toupper(script_name), "> ", "topTopTable_genes_enrichGO is NULL\n")
  printAndLog(txt, pipLogFile)
}
outFile <- paste0(curr_outFold, "/", "topTopTable_genes_enrichGO_resultDT.Rdata")
save(topTopTable_genes_enrichGO_resultDT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

################################################################################################################
######################################################## Gene Set Enrichment Analysis of Gene Ontology (gseGO)
################################################################################################################
# not sure but I think this is not appropriate
# this is useful when using the full set of genes, so that we don't need thresholding the genes of interest
# but in my case I already have a kind of threshold
# UPDATE: do gseGO for the full list of genes and then compare to the enrichGO for the reduced list of genes from top-ranking TADs/TopTable

### PERFORM GENE SET ON ALL THE GENES WITH ADJ. PVAL AS SCORES
all_genes_scoreDT <- curr_topTable_DT[, c("genes","P.Value", "adj.P.Val", "entrezID")]
all_genes_scoreDT <- all_genes_scoreDT[order(all_genes_scoreDT$adj.P.Val, all_genes_scoreDT$P.Value),]
all_genes_scores <- setNames(-log10(all_genes_scoreDT$adj.P.Val), all_genes_scoreDT$entrezID)
stopifnot(!is.na(names(all_genes_scores)))


# outFile <- paste0(curr_outFold, "/", "tmp_all_genes_scores_topTable.Rdata")
# save(all_genes_scores, file = outFile)
# stop("--ok\n")

# HERE USE OR NOT MIN SIZE GENE SET ???
all_genes_gseGO <- gseGO(geneList = all_genes_scores,
                         OrgDb = org.Hs.eg.db,
                         nPerm = funcTopTable_nPerm_gseGO,
                         pvalueCutoff = 1,
                         verbose=F,
                         ont="ALL")

if(!is.null(all_genes_gseGO)) {
  outFile <- paste0(curr_outFold, "/", "all_genes_gseGO_dotplot.", plotType)
  # do.call(plotType, list(file=outFile, height=myHeight, width=myWidth*2))
  # dotplot(all_genes_gseGO, showCategory=funcTopTable_nShowPlot)
  # foo <- dev.off()
  p <- dotplot(all_genes_gseGO, showCategory=funcTopTable_nShowPlot)
  ggsave(p, filename = outFile, height=myHeight, width=myWidth*2)
  cat(paste0("... written: ", outFile, "\n"))
  
  all_genes_gseGO_resultDT <- all_genes_gseGO@result
  all_genes_gseGO_resultDT <- all_genes_gseGO_resultDT[order(all_genes_gseGO_resultDT$p.adjust, all_genes_gseGO_resultDT$pvalue),]
  
  all_genes_signifGO <- all_genes_gseGO_resultDT$ID[all_genes_gseGO_resultDT$p.adjust <= p_adj_topGOthresh]
  txt <- paste0(toupper(script_name), "> found signif GO for all genes (gseGO): ", length(all_genes_signifGO), "\n")
  printAndLog(txt, pipLogFile)
  
  all_genes_nTopGO <- all_genes_gseGO_resultDT$ID[1:min(nTopGO, nrow(all_genes_gseGO_resultDT))]
  txt <- paste0(toupper(script_name), "> found enough nTop gseGO for all genes: ", length(all_genes_nTopGO), "/", nTopGO, "\n")
  printAndLog(txt, pipLogFile)
  
  # all_genes_gseGO_topResultDT <- all_genes_gseGO_resultDT[1:min(nTopGO, nrow(all_genes_gseGO_resultDT)),]
  # outFile <- paste0(curr_outFold, "/", "all_genes_gseGO_topResultDT.Rdata")
  # save(all_genes_gseGO_topResultDT, file = outFile)
  # cat(paste0("... written: ", outFile, "\n"))
  outFile <- paste0(curr_outFold, "/", "all_genes_gseGO_custombarplot_padj_withNbr.", plotType)
  try(dev.off())
  do.call(plotType, list(file=outFile, height=myHeight*1.5, width=myWidth*2))
  barplot_funcEnrich_results(resultDT = all_genes_gseGO_resultDT, 
                             nTop = funcTADs_nShowPlot, 
                             xVar = "p.adjust", 
                             g2t_DT = gene2tadDT, 
                             DE_genes = topTopTable_genes, 
                             plotHoriz=TRUE,
                             barcol="slategray")
  title(main = paste0(customBarplotTitle, " - gseGO"))
  mtext(text = customBarplotSubtit, side = 3, font = 3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  

}else{
  all_genes_gseGO_resultDT <- NULL
  txt <- paste0(toupper(script_name), "> ", "all_genes_gseGO is NULL\n")
  printAndLog(txt, pipLogFile)
}

outFile <- paste0(curr_outFold, "/", "all_genes_gseGO_resultDT.Rdata")
save(all_genes_gseGO_resultDT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

                      


############################################################################################################################
############################################################################################################################
############################################################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))



