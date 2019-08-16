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

script_name <- "23_funcEnrich_TopTableGenes_TADsGenes"
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

stopifnot(exists("nTopTADs_Fisher"))
stopifnot(exists("nTopTADs"))
stopifnot(exists("nTopGSEA"))
stopifnot(exists("nTopGO"))
stopifnot(exists("nTopKEGG"))
stopifnot(exists("nPerm_gseGO"))
stopifnot(exists("nPerm_gseKEGG"))
stopifnot(exists("nPerm_GSEA"))
stopifnot(exists("p_adj_topKEGGthresh"))
stopifnot(exists("p_adj_topGOthresh"))
stopifnot(exists("p_adj_DE_thresh_Fisher"))

cat(paste0("... found setting: nTopTADs = ", nTopTADs, "\n"))
cat(paste0("... found setting: nTopTADs_Fisher = ", nTopTADs_Fisher, "\n"))
cat(paste0("... found setting: nTopGSEA = ", nTopGSEA, "\n"))
cat(paste0("... found setting: nTopGO = ", nTopGO, "\n"))
cat(paste0("... found setting: nTopKEGG = ", nTopKEGG, "\n"))
cat(paste0("... found setting: nPerm_gseGO = ", nPerm_gseGO, "\n"))
cat(paste0("... found setting: nPerm_gseKEGG = ", nPerm_gseKEGG, "\n"))
cat(paste0("... found setting: nPerm_GSEA = ", nPerm_GSEA, "\n"))
cat(paste0("... found setting: p_adj_DE_thresh_Fisher = ", p_adj_DE_thresh_Fisher, "\n"))
cat(paste0("... found setting: p_adj_topGOthresh = ", p_adj_topGOthresh, "\n"))
cat(paste0("... found setting: p_adj_topKEGGthresh = ", p_adj_topKEGGthresh, "\n"))

#****************************** FOR DEBUG
# nTopTADs <- 10
# nTopGSEA <- 10
# nPerm_gseGO <- 100
# nPerm_GSEA <- 100
# p_adj_topGOthresh <- 0.05
# nTopGO <- 20
#******************************

nShowPlot <- 10
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

topTADs <- names(rank(adj_pval_combined, ties="min")[rank(adj_pval_combined, ties="min") <= nTopTADs])
topTADs_genes <- as.character(c(unlist(sapply(topTADs, function(x) curr_g2t_DT$entrezID[curr_g2t_DT$region == x]))))
topTADs_genes_symb <- as.character(unlist(sapply(topTADs_genes, function(x) names(geneList[geneList == x]))))

# select the same number of genes, taking the same # based on limma pval
curr_topTable_DT <- curr_topTable_DT[order(curr_topTable_DT$adj.P.Val, curr_topTable_DT$P.Value),]

topTopTable_genes_symb <- curr_topTable_DT$genes[1:length(topTADs_genes)]
topTopTable_genes <- unlist(sapply(topTopTable_genes_symb, function(x) geneList[x]))

################################################################################################################
######################################################## HYPERGEO tests -> 1 single test of overrepresentation of DE genes in pooled topTADs
################################################################################################################

        topTADs_Fisher <- names(rank(adj_pval_combined, ties="min")[rank(adj_pval_combined, ties="min") <= nTopTADs_Fisher])
        topTADs_Fisher_genes <- as.character(c(unlist(sapply(topTADs_Fisher, function(x) curr_g2t_DT$entrezID[curr_g2t_DT$region == x]))))
        topTADs_Fisher_genes_symb <- as.character(unlist(sapply(topTADs_Fisher_genes, function(x) names(geneList[geneList == x]))))

        cat(paste0("... # of genes derived from topTADs_Fisher: ", length(topTADs_Fisher_genes_symb), "\n"))

        # select the DE genes
        curr_topTable_DT <- curr_topTable_DT[order(curr_topTable_DT$adj.P.Val, curr_topTable_DT$P.Value),]
        TopTable_DE_genes_symb <- curr_topTable_DT$genes[curr_topTable_DT$adj.P.Val <= p_adj_DE_thresh_Fisher]
        TopTable_DE_genes <- unlist(sapply(TopTable_DE_genes_symb, function(x) geneList[x]))

        TopTable_notDE_genes_symb <- curr_topTable_DT$genes[curr_topTable_DT$adj.P.Val >= p_adj_DE_thresh_Fisher]
        TopTable_notDE_genes <- unlist(sapply(TopTable_notDE_genes_symb, function(x) geneList[x]))


        ###
        #                    DE gene | not DE gene
        #  -----------------------------------
        #  intopTADs     |    #      |   #
        # -------------------------------------
        #  notIntopTADs  |    #      |   #
        #
        #
        # filled by column
        DE_inTop <- sum(TopTable_DE_genes %in% topTADs_Fisher_genes)
        DE_notInTop <- sum(!TopTable_DE_genes %in% topTADs_Fisher_genes)
        stopifnot(sum(topTADs_Fisher_genes %in% TopTable_DE_genes) == DE_inTop)
        notDE_inTop <- sum(!topTADs_Fisher_genes %in% TopTable_DE_genes)
        stopifnot(sum(TopTable_notDE_genes %in% topTADs_Fisher_genes) == notDE_inTop)
        notDE_notInTop <- sum(! TopTable_notDE_genes %in% topTADs_Fisher_genes)

        testMatrix <- data.frame(DE_genes = c(DE_inTop, DE_notInTop),
                                 notDE_genes = c(notDE_inTop, notDE_notInTop),
                                 row.names=c("inTopTADs", "notInTopTADs"))
        stopifnot(sum(testMatrix) == length(curr_topTable_DT$adj.P.Val))
        print(testMatrix)

        FT_result <- fisher.test(testMatrix, alternative = "greater")

        DEgenes_in_topTADs_fisherTest <- list(
          FT_matrix = testMatrix,
          FT_pval = FT_result$p.value,
          FT_odds = FT_result$estimate,
          FT_settings = c(nTopTADs = nTopTADs_Fisher, DE_p_adj_thresh = p_adj_DE_thresh_Fisher)
        )

        outFile <- paste0(curr_outFold, "/", "DEgenes_in_topTADs_fisherTest.Rdata")
        save(DEgenes_in_topTADs_fisherTest, file = outFile)
        cat(paste0("... written: ", outFile, "\n"))


################################################################################################################
######################################################## HYPERGEO tests -> 1 test of overrepresentation of DE genes separately for each topTADs
################################################################################################################

    DEgenes_in_topTAD_separate_fisherTest <- foreach(currTAD = topTADs_Fisher) %dopar% {

      currTADs_genes <- as.character(c(unlist(sapply(currTAD, function(x) curr_g2t_DT$entrezID[curr_g2t_DT$region == x]))))

      # filled by column
      DE_inCurrTAD <- sum(TopTable_DE_genes %in% currTADs_genes)
      DE_notInCurrTAD <- sum(!TopTable_DE_genes %in% currTADs_genes)
      stopifnot(sum(currTADs_genes %in% TopTable_DE_genes) == DE_inCurrTAD)
      notDE_inCurrTAD <- sum(!currTADs_genes %in% TopTable_DE_genes)
      stopifnot(sum(TopTable_notDE_genes %in% currTADs_genes) == notDE_inCurrTAD)
      notDE_notInCurrTAD <- sum(! TopTable_notDE_genes %in% currTADs_genes)

      testMatrix <- data.frame(DE_genes = c(DE_inCurrTAD, DE_notInCurrTAD),
                               notDE_genes = c(notDE_inCurrTAD, notDE_notInCurrTAD),
                               row.names=c("inCurrTAD", "notInCurrTAD"))
      print(testMatrix)

      FT_result <- fisher.test(testMatrix, alternative = "greater")

      c(FT_pval = FT_result$p.value,
        FT_odds = FT_result$estimate,
        set_size = length(currTADs_genes)
      )

    }
    names(DEgenes_in_topTAD_separate_fisherTest) <- topTADs_Fisher

    outFile <- paste0(curr_outFold, "/", "DEgenes_in_topTAD_separate_fisherTest.Rdata")
    save(DEgenes_in_topTAD_separate_fisherTest, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))


################################################################################################################
######################################################## GSEA FOR THE TOP TADs
################################################################################################################

                      stopifnot(length(topTADs) >= nTopGSEA)

                      tad_topTable_DT <- curr_topTable_DT
                      tad_topTable_DT$topTAD <- unlist(sapply(tad_topTable_DT$entrezID, function(x) {
                        x_tad <- curr_g2t_DT$region[curr_g2t_DT$entrezID == x]
                        i_top <- which(topTADs == x_tad)
                        if(length(i_top) == 0) {
                          return("NO")
                        } else {
                          # return(paste0("topTAD", i_top))
                          return(topTADs[i_top])
                        }
                      }))
                      tad_topTable_DT <- tad_topTable_DT[order(tad_topTable_DT$adj.P.Val, tad_topTable_DT$P.Value, decreasing = T),]
                      score_vect <- setNames(tad_topTable_DT$adj.P.Val, tad_topTable_DT$genes)
                      term2gene <- tad_topTable_DT[,c("topTAD", "genes")]
                      colnames(term2gene) <- c("term", "name")
                      topTADs_GSEA <- GSEA(geneList= score_vect, TERM2GENE=term2gene, nPerm=nPerm_GSEA, pvalueCutoff=1, pAdjustMethod= "BH", minGSSize = 3)
                      # gseaplot(topTADs_GSEA, topTADs[1])
                      # summary(topTADs_GSEA)

                      topTADs_adjPvalGSEA <- setNames(topTADs_GSEA@result$p.adjust, topTADs_GSEA@result$ID)
                      topTADs_adjPvalGSEA <- topTADs_adjPvalGSEA[match(topTADs,names(topTADs_adjPvalGSEA))]
                      stopifnot(names(topTADs_adjPvalGSEA) == topTADs)

                      topTADs_GSEA_resultDT <- topTADs_GSEA@result
                      topTADs_GSEA_topResultDT <- topTADs_GSEA_resultDT[1:min(nTopGO, nrow(topTADs_GSEA_resultDT)),]

                      outFile <- paste0(curr_outFold, "/", "topTADs_GSEA_topResultDT.Rdata")
                      save(topTADs_GSEA_topResultDT, file = outFile)
                      cat(paste0("... written: ", outFile, "\n"))


                      ################################################################################################################
                      ######################################################## GO Enrichment Analysis of a gene set (enrichGO)
                      ################################################################################################################

                      # use an over-representation test rather than an enrichment, because here I have already kind of thresholded list
                      # of genes as I take only genes from the topTADs (I don't care about their scores)


                      ### FOR THE topTADs_genes
                      topTADs_genes_scoreDT <- curr_topTable_DT[curr_topTable_DT$genes %in% topTADs_genes_symb, c("genes","adj.P.Val", "entrezID")]

                      # for the moment, background genes set to universe
                      topTADs_genes_enrichGO <- enrichGO(gene = topTADs_genes_scoreDT$entrezID,
                                                         OrgDb = org.Hs.eg.db,
                                                         # change default min gene set size (10)
                                                         minGSSize = 1,
                                                         keytype = "ENTREZID",
                                                         pvalueCutoff = 1,
                                                         pAdjustMethod = "BH",
                                                         ont="ALL")

                      outFile <- paste0(curr_outFold, "/", "topTADs_genes_enrichGO_dotplot.", plotType)
                      # do.call(plotType, list(file=outFile, height=myHeight, width=myWidth*2))
                      # dotplot(topTADs_genes_enrichGO, showCategory=nShowPlot)
                      # foo <- dev.off()
                      p <- dotplot(topTADs_genes_enrichGO, showCategory=nShowPlot)
                      ggsave(p, filename = outFile, height=myHeight, width=myWidth*2)
                      cat(paste0("... written: ", outFile, "\n"))

                      outFile <- paste0(curr_outFold, "/", "topTADs_genes_enrichGO_barplot.", plotType)
                      # do.call(plotType, list(file=outFile, height=myHeight, width=myWidth*2))
                      # barplot(topTADs_genes_enrichGO, showCategory=nShowPlot, order=T)
                      # foo <- dev.off()
                      p <- barplot(topTADs_genes_enrichGO, showCategory=nShowPlot, order=T)
                      ggsave(p, filename = outFile, height=myHeight, width=myWidth*2)
                      cat(paste0("... written: ", outFile, "\n"))
                      
                      topTADs_genes_enrichGO_resultDT <- topTADs_genes_enrichGO@result
                      topTADs_genes_enrichGO_resultDT <- topTADs_genes_enrichGO_resultDT[order(topTADs_genes_enrichGO_resultDT$p.adjust, topTADs_genes_enrichGO_resultDT$pvalue),]

                      topTADs_signifGO <- topTADs_genes_enrichGO_resultDT$ID[topTADs_genes_enrichGO_resultDT$p.adjust <= p_adj_topGOthresh]
                      txt <- paste0(toupper(script_name), "> found signif enrichGO for topTADs genes: ", length(topTADs_signifGO), "\n")
                      printAndLog(txt, pipLogFile)

                      topTADs_nTopGO <- topTADs_genes_enrichGO_resultDT$ID[1:min(nTopGO, nrow(topTADs_genes_enrichGO_resultDT))]
                      txt <- paste0(toupper(script_name), "> found enough nTop enrichGO for topTADs genes: ", length(topTADs_nTopGO), "/", nTopGO, "\n")
                      printAndLog(txt, pipLogFile)

                      topTADs_genes_enrichGO_topResultDT <- topTADs_genes_enrichGO_resultDT[1:min(nTopGO, nrow(topTADs_genes_enrichGO_resultDT)),]

                      outFile <- paste0(curr_outFold, "/", "topTADs_genes_enrichGO_topResultDT.Rdata")
                      save(topTADs_genes_enrichGO_topResultDT, file = outFile)
                      cat(paste0("... written: ", outFile, "\n"))


                      ### FOR THE topTopTable_genes
                      topTopTable_genes_scoreDT <- curr_topTable_DT[curr_topTable_DT$genes %in% topTopTable_genes_symb, c("genes","adj.P.Val", "entrezID")]

                      topTopTable_genes_enrichGO <- enrichGO(gene = topTopTable_genes_scoreDT$entrezID,
                                                             OrgDb = org.Hs.eg.db,
                                                             # change default min gene set size (10)
                                                             minGSSize = 1,
                                                             keytype = "ENTREZID",
                                                             pvalueCutoff = 1,
                                                             pAdjustMethod = "BH",
                                                             ont="ALL")

                      outFile <- paste0(curr_outFold, "/", "topTopTable_genes_enrichGO_dotplot.", plotType)
                      # do.call(plotType, list(file=outFile, height=myHeight, width=myWidth*2))
                      # dotplot(topTopTable_genes_enrichGO, showCategory=nShowPlot)
                      # foo <- dev.off()
                      p <- dotplot(topTopTable_genes_enrichGO, showCategory=nShowPlot)
                      ggsave(p, filename = outFile, height=myHeight, width=myWidth*2)
                      cat(paste0("... written: ", outFile, "\n"))
                      
                      outFile <- paste0(curr_outFold, "/", "topTopTable_genes_enrichGO_barplot.", plotType)
                      # do.call(plotType, list(file=outFile, height=myHeight, width=myWidth*2))
                      # barplot(topTopTable_genes_enrichGO, showCategory=nShowPlot, order=T)
                      # foo <- dev.off()
                      p <- barplot(topTopTable_genes_enrichGO, showCategory=nShowPlot, order=T)
                      ggsave(p, filename = outFile, height=myHeight, width=myWidth*2)
                      cat(paste0("... written: ", outFile, "\n"))
                      
                      topTopTable_genes_enrichGO_resultDT <- topTopTable_genes_enrichGO@result
                      topTopTable_genes_enrichGO_resultDT <- topTopTable_genes_enrichGO_resultDT[order(topTopTable_genes_enrichGO_resultDT$p.adjust, topTopTable_genes_enrichGO_resultDT$pvalue),]


                      topTopTable_signifGO <- topTopTable_genes_enrichGO_resultDT$ID[topTopTable_genes_enrichGO_resultDT$p.adjust <= p_adj_topGOthresh]
                      txt <- paste0(toupper(script_name), "> found signif enrichGO for topTopTable genes: ", length(topTopTable_signifGO), "\n")
                      printAndLog(txt, pipLogFile)

                      topTopTable_nTopGO <- topTopTable_genes_enrichGO_resultDT$ID[1:min(nTopGO, nrow(topTopTable_genes_enrichGO_resultDT))]
                      txt <- paste0(toupper(script_name), "> found enough nTop enrichGO for topTADs genes: ", length(topTADs_nTopGO), "/", nTopGO, "\n")
                      printAndLog(txt, pipLogFile)


                      topTopTable_genes_enrichGO_topResultDT <- topTopTable_genes_enrichGO_resultDT[1:min(nTopGO, nrow(topTopTable_genes_enrichGO_resultDT)),]
                      outFile <- paste0(curr_outFold, "/", "topTopTable_genes_enrichGO_topResultDT.Rdata")
                      save(topTopTable_genes_enrichGO_topResultDT, file = outFile)
                      cat(paste0("... written: ", outFile, "\n"))


                      intersectSignifGO_enrichGO <- intersect(topTADs_signifGO, topTopTable_signifGO)
                      intersectTopGO_enrichGO <- intersect(topTADs_nTopGO, topTopTable_nTopGO)

                      outFile <- paste0(curr_outFold, "/", "intersectSignifGO_enrichGO.Rdata")
                      save(intersectSignifGO_enrichGO, file = outFile)
                      cat(paste0("... written: ", outFile, "\n"))

                      outFile <- paste0(curr_outFold, "/", "intersectTopGO_enrichGO.Rdata")
                      save(intersectTopGO_enrichGO, file = outFile)
                      cat(paste0("... written: ", outFile, "\n"))

                      txt <- paste0(toupper(script_name), "> found intersect signif GO topTADs/topTopTable: ", length(intersectSignifGO_enrichGO), " (adj.pval <= ", p_adj_topGOthresh, ")\n")
                      printAndLog(txt, pipLogFile)

                      txt <- paste0(toupper(script_name), "> found intersect top GO topTADs/topTopTable: ", length(intersectTopGO_enrichGO), " (nTop = ", nTopGO, ")\n")
                      printAndLog(txt, pipLogFile)


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

                      # HERE USE OR NOT MIN SIZE GENE SET ???
                      all_genes_gseGO <- gseGO(geneList = all_genes_scores,
                                               OrgDb = org.Hs.eg.db,
                                               nPerm = nPerm_gseGO,
                                               pvalueCutoff = 1,
                                               verbose=F,
                                               ont="ALL")

                      outFile <- paste0(curr_outFold, "/", "all_genes_gseGO_dotplot.", plotType)
                      # do.call(plotType, list(file=outFile, height=myHeight, width=myWidth*2))
                      # dotplot(all_genes_gseGO, showCategory=nShowPlot)
                      # foo <- dev.off()
                      p <- dotplot(all_genes_gseGO, showCategory=nShowPlot)
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

                      all_genes_gseGO_topResultDT <- all_genes_gseGO_resultDT[1:min(nTopGO, nrow(all_genes_gseGO_resultDT)),]

                      outFile <- paste0(curr_outFold, "/", "all_genes_gseGO_topResultDT.Rdata")
                      save(all_genes_gseGO_topResultDT, file = outFile)
                      cat(paste0("... written: ", outFile, "\n"))

                      intersectSignifGO_gseGO <- intersect(topTADs_signifGO, all_genes_signifGO)
                      intersectTopGO_gseGO <- intersect(topTADs_nTopGO, all_genes_nTopGO)

                      outFile <- paste0(curr_outFold, "/", "intersectSignifGO_gseGO.Rdata")
                      save(intersectSignifGO_gseGO, file = outFile)
                      cat(paste0("... written: ", outFile, "\n"))

                      outFile <- paste0(curr_outFold, "/", "intersectTopGO_gseGO.Rdata")
                      save(intersectTopGO_gseGO, file = outFile)
                      cat(paste0("... written: ", outFile, "\n"))

                      txt <- paste0(toupper(script_name), "> found intersect signif GO topTADs/all genes gseGO: ", length(intersectSignifGO_gseGO), " (adj.pval <= ", p_adj_topGOthresh, ")\n")
                      printAndLog(txt, pipLogFile)

                      txt <- paste0(toupper(script_name), "> found intersect top GO topTADs/all genes gseGO: ", length(intersectTopGO_gseGO), " (nTop = ", nTopGO, ")\n")
                      printAndLog(txt, pipLogFile)



                      ################################################################################################################
                      ######################################################## KEGG Enrichment Analysis of a gene set (enrichKEGG)
                      ################################################################################################################

                      # use an over-representation test rather than an enrichment, because here I have already kind of thresholded list
                      # of genes as I take only genes from the topTADs (I don't care about their scores)


                      ### FOR THE topTADs_genes
                      topTADs_genes_scoreDT <- curr_topTable_DT[curr_topTable_DT$genes %in% topTADs_genes_symb, c("genes","adj.P.Val", "entrezID")]

                      # for the moment, background genes set to universe
                      topTADs_genes_enrichKEGG <- enrichKEGG(gene = topTADs_genes_scoreDT$entrezID,
                                                             organism = "hsa",
                                                             keyType = "ncbi-geneid",
                                                             pvalueCutoff = 1,
                                                             minGSSize = 1,
                                                             pAdjustMethod = "BH"
                                                            )

                      outFile <- paste0(curr_outFold, "/", "topTADs_genes_enrichKEGG_dotplot.", plotType)
                      # do.call(plotType, list(file=outFile, height=myHeight, width=myWidth*2))
                      # dotplot(topTADs_genes_enrichKEGG, showCategory=nShowPlot)
                      # foo <- dev.off()
                      p <- dotplot(topTADs_genes_enrichKEGG, showCategory=nShowPlot)
                      ggsave(p, filename = outFile, height=myHeight, width=myWidth*2)
                      cat(paste0("... written: ", outFile, "\n"))
                      
                      outFile <- paste0(curr_outFold, "/", "topTADs_genes_enrichKEGG_barplot.", plotType)
                      # do.call(plotType, list(file=outFile, height=myHeight, width=myWidth*2))
                      # barplot(topTADs_genes_enrichKEGG, showCategory=nShowPlot, order=T)
                      # foo <- dev.off()
                      p <- barplot(topTADs_genes_enrichKEGG, showCategory=nShowPlot, order=T)
                      ggsave(p, filename = outFile, height=myHeight, width=myWidth*2)
                      cat(paste0("... written: ", outFile, "\n"))
                      
                      topTADs_genes_enrichKEGG_resultDT <- topTADs_genes_enrichKEGG@result
                      topTADs_genes_enrichKEGG_resultDT <- topTADs_genes_enrichKEGG_resultDT[order(topTADs_genes_enrichKEGG_resultDT$p.adjust, topTADs_genes_enrichKEGG_resultDT$pvalue),]

                      topTADs_signifKEGG <- topTADs_genes_enrichKEGG_resultDT$ID[topTADs_genes_enrichKEGG_resultDT$p.adjust <= p_adj_topKEGGthresh]
                      txt <- paste0(toupper(script_name), "> found signif enrichKEGG for topTADs genes: ", length(topTADs_signifKEGG), "\n")
                      printAndLog(txt, pipLogFile)

                      topTADs_nTopKEGG <- topTADs_genes_enrichKEGG_resultDT$ID[1:min(nTopKEGG, nrow(topTADs_genes_enrichKEGG_resultDT))]
                      txt <- paste0(toupper(script_name), "> found enough nTop enrichKEGG for topTADs genes: ", length(topTADs_nTopKEGG), "/", nTopKEGG, "\n")
                      printAndLog(txt, pipLogFile)

                      topTADs_genes_enrichKEGG_topResultDT <- topTADs_genes_enrichKEGG_resultDT[1:min(nTopKEGG, nrow(topTADs_genes_enrichKEGG_resultDT)),]

                      outFile <- paste0(curr_outFold, "/", "topTADs_genes_enrichKEGG_topResultDT.Rdata")
                      save(topTADs_genes_enrichKEGG_topResultDT, file = outFile)
                      cat(paste0("... written: ", outFile, "\n"))


                      ### FOR THE topTopTable_genes
                      topTopTable_genes_scoreDT <- curr_topTable_DT[curr_topTable_DT$genes %in% topTopTable_genes_symb, c("genes","adj.P.Val", "entrezID")]

                      topTopTable_genes_enrichKEGG <- enrichKEGG(gene = topTopTable_genes_scoreDT$entrezID,
                                                                 organism = "hsa",
                                                                 keyType = "ncbi-geneid",
                                                                 pvalueCutoff = 1,
                                                                 minGSSize = 1,
                                                                 pAdjustMethod = "BH"
                                                                 )

                      outFile <- paste0(curr_outFold, "/", "topTopTable_genes_enrichKEGG_dotplot.", plotType)
                      # do.call(plotType, list(file=outFile, height=myHeight, width=myWidth*2))
                      # dotplot(topTopTable_genes_enrichKEGG, showCategory=nShowPlot)
                      # foo <- dev.off()
                      p <- dotplot(topTopTable_genes_enrichKEGG, showCategory=nShowPlot)
                      ggsave(p, filename = outFile, height=myHeight, width=myWidth*2)
                      cat(paste0("... written: ", outFile, "\n"))
                      
                      outFile <- paste0(curr_outFold, "/", "topTopTable_genes_enrichKEGG_barplot.", plotType)
                      # do.call(plotType, list(file=outFile, height=myHeight, width=myWidth*2))
                      # barplot(topTopTable_genes_enrichKEGG, showCategory=nShowPlot, order=T)
                      # foo <- dev.off()
                      p <- barplot(topTopTable_genes_enrichKEGG, showCategory=nShowPlot, order=T)
                      ggsave(p, filename = outFile, height=myHeight, width=myWidth*2)
                      cat(paste0("... written: ", outFile, "\n"))
                      
                      topTopTable_genes_enrichKEGG_resultDT <- topTopTable_genes_enrichKEGG@result
                      topTopTable_genes_enrichKEGG_resultDT <- topTopTable_genes_enrichKEGG_resultDT[order(topTopTable_genes_enrichKEGG_resultDT$p.adjust, topTopTable_genes_enrichKEGG_resultDT$pvalue),]

                      topTopTable_signifKEGG <- topTopTable_genes_enrichKEGG_resultDT$ID[topTopTable_genes_enrichKEGG_resultDT$p.adjust <= p_adj_topKEGGthresh]
                      txt <- paste0(toupper(script_name), "> found signif enrichKEGG for topTopTable genes: ", length(topTopTable_signifKEGG), "\n")
                      printAndLog(txt, pipLogFile)

                      topTopTable_nTopKEGG <- topTopTable_genes_enrichKEGG_resultDT$ID[1:min(nTopKEGG, nrow(topTopTable_genes_enrichKEGG_resultDT))]
                      txt <- paste0(toupper(script_name), "> found enough nTop enrichKEGG for topTADs genes: ", length(topTADs_nTopKEGG), "/", nTopKEGG, "\n")
                      printAndLog(txt, pipLogFile)

                      topTopTable_genes_enrichKEGG_topResultDT <- topTopTable_genes_enrichKEGG_resultDT[1:min(nTopKEGG, nrow(topTopTable_genes_enrichKEGG_resultDT)),]
                      outFile <- paste0(curr_outFold, "/", "topTopTable_genes_enrichKEGG_topResultDT.Rdata")
                      save(topTopTable_genes_enrichKEGG_topResultDT, file = outFile)
                      cat(paste0("... written: ", outFile, "\n"))


                      intersectSignifKEGG_enrichKEGG <- intersect(topTADs_signifKEGG, topTopTable_signifKEGG)
                      intersectTopKEGG_enrichKEGG <- intersect(topTADs_nTopKEGG, topTopTable_nTopKEGG)

                      outFile <- paste0(curr_outFold, "/", "intersectSignifKEGG_enrichKEGG.Rdata")
                      save(intersectSignifKEGG_enrichKEGG, file = outFile)
                      cat(paste0("... written: ", outFile, "\n"))

                      outFile <- paste0(curr_outFold, "/", "intersectTopKEGG_enrichKEGG.Rdata")
                      save(intersectTopKEGG_enrichKEGG, file = outFile)
                      cat(paste0("... written: ", outFile, "\n"))


                      txt <- paste0(toupper(script_name), "> found intersect signif KEGG topTADs/topTopTable: ", length(intersectSignifKEGG_enrichKEGG), " (adj.pval <= ", p_adj_topKEGGthresh, ")\n")
                      printAndLog(txt, pipLogFile)

                      txt <- paste0(toupper(script_name), "> found intersect top KEGG topTADs/topTopTable: ", length(intersectTopKEGG_enrichKEGG), " (nTop = ", nTopKEGG, ")\n")
                      printAndLog(txt, pipLogFile)


                      ################################################################################################################
                      ######################################################## Gene Set Enrichment Analysis of KEGG pathways (gseKEGG)
                      ################################################################################################################
                      # not sure but I think this is not appropriate
                      # this is useful when using the full set of genes, so that we don't need thresholding the genes of interest
                      # but in my case I already have a kind of threshold
                      # UPDATE: do gseKEGG for the full list of genes and then compare to the enrichKEGG for the reduced list of genes from top-ranking TADs/TopTable

                      ### PERFORM GENE SET ON ALL THE GENES WITH ADJ. PVAL AS SCORES
                      all_genes_scoreDT <- curr_topTable_DT[, c("genes","P.Value", "adj.P.Val", "entrezID")]
                      all_genes_scoreDT <- all_genes_scoreDT[order(all_genes_scoreDT$adj.P.Val, all_genes_scoreDT$P.Value),]
                      all_genes_scores <- setNames(-log10(all_genes_scoreDT$adj.P.Val), all_genes_scoreDT$entrezID)
                      stopifnot(!is.na(names(all_genes_scores)))

                      # HERE USE OR NOT MIN SIZE GENE SET ???
                      all_genes_gseKEGG <- gseKEGG(geneList = all_genes_scores,
                                                   organism = "hsa",
                                                   keyType = "ncbi-geneid",
                                                   nPerm = nPerm_gseKEGG,
                                                   pvalueCutoff = 1,
                                                   verbose=F
                                                   )

                      outFile <- paste0(curr_outFold, "/", "all_genes_gseKEGG_dotplot.", plotType)
                      # do.call(plotType, list(file=outFile, height=myHeight, width=myWidth*2))
                      # dotplot(all_genes_gseKEGG, showCategory=nShowPlot)
                      # foo <- dev.off()
                      p <- dotplot(all_genes_gseKEGG, showCategory=nShowPlot)
                      ggsave(p, filename = outFile, height=myHeight, width=myWidth*2)
                      cat(paste0("... written: ", outFile, "\n"))
                      
                      all_genes_gseKEGG_resultDT <- all_genes_gseKEGG@result
                      all_genes_gseKEGG_resultDT <- all_genes_gseKEGG_resultDT[order(all_genes_gseKEGG_resultDT$p.adjust, all_genes_gseKEGG_resultDT$pvalue),]

                      all_genes_signifKEGG <- all_genes_gseKEGG_resultDT$ID[all_genes_gseKEGG_resultDT$p.adjust <= p_adj_topKEGGthresh]
                      txt <- paste0(toupper(script_name), "> found signif KEGG for all genes (gseKEGG): ", length(all_genes_signifKEGG), "\n")
                      printAndLog(txt, pipLogFile)

                      all_genes_nTopKEGG <- all_genes_gseKEGG_resultDT$ID[1:min(nTopKEGG, nrow(all_genes_gseKEGG_resultDT))]
                      txt <- paste0(toupper(script_name), "> found enough nTop gseKEGG for all genes: ", length(all_genes_nTopKEGG), "/", nTopKEGG, "\n")
                      printAndLog(txt, pipLogFile)

                      all_genes_gseKEGG_topResultDT <- all_genes_gseKEGG_resultDT[1:min(nTopKEGG, nrow(all_genes_gseKEGG_resultDT)),]

                      outFile <- paste0(curr_outFold, "/", "all_genes_gseKEGG_topResultDT.Rdata")
                      save(all_genes_gseKEGG_topResultDT, file = outFile)
                      cat(paste0("... written: ", outFile, "\n"))


                      intersectSignifKEGG_gseKEGG <- intersect(topTADs_signifKEGG, all_genes_signifKEGG)
                      intersectTopKEGG_gseKEGG <- intersect(topTADs_nTopKEGG, all_genes_nTopKEGG)


                      outFile <- paste0(curr_outFold, "/", "intersectSignifKEGG_gseKEGG.Rdata")
                      save(intersectSignifKEGG_gseKEGG, file = outFile)
                      cat(paste0("... written: ", outFile, "\n"))

                      outFile <- paste0(curr_outFold, "/", "intersectTopKEGG_gseKEGG.Rdata")
                      save(intersectTopKEGG_gseKEGG, file = outFile)
                      cat(paste0("... written: ", outFile, "\n"))

                      txt <- paste0(toupper(script_name), "> found intersect signif KEGG topTADs/all genes gseKEGG: ", length(intersectSignifKEGG_gseKEGG), " (adj.pval <= ", p_adj_topKEGGthresh, ")\n")
                      printAndLog(txt, pipLogFile)

                      txt <- paste0(toupper(script_name), "> found intersect top KEGG topTADs/all genes gseKEGG: ", length(intersectTopKEGG_gseKEGG), " (nTop = ", nTopKEGG, ")\n")
                      printAndLog(txt, pipLogFile)


############################################################################################################################
############################################################################################################################
############################################################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))



