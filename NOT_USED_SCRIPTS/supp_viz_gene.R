startTime <- Sys.time()
script_name <- "supp_vize_gene"
script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script3_name <- "3_runMeanTADLogFC"
script11_name <- "11_runEmpPvalCombined"

# Rscript supp_viz_gene.R run_settings_GSE73765_noninf_salm.R FGF5
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(1)

#### FOR DEBUG  # !!!!!!!!!!!!!!!!!!!!!!!
# run_setting_file <- "run_settings_GSE73765_noninf_salm.R" 
# all_geneSymbol <- "FGF5"
# run_setting_file <-  "run_settings_GSE94736_old_young.R"
# all_geneSymbol <- "SLC16A14"
####################################

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 2)
  stop("ERROR: missing some command line arguments\n")

run_setting_file <- args[1]
all_geneSymbol <- args[2:length(args)]

# cat(paste0("setDir = ", setDir, "\n"))
source("main_settings.R") # setDir is the main_settings not in run_settings
source(run_setting_file)
source("TAD_DE_utils.R")

curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

plotType <- "svg"
plotWidth <- 10
plotHeight <- 7

all_emp_pval_combined <- eval(parse(text = load(paste0(pipOutFold, "/", script11_name, "/emp_pval_combined.Rdata"))))
all_emp_pval_combined <- rank(all_emp_pval_combined)

######################################################################
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

symbDT <- read.delim(symbolDT_file, header=T, stringsAsFactors = F, col.names=c("entrezID", "symbol"))
symbDT$entrezID <- as.character(symbDT$entrezID)

entrez2symbDT <- read.delim(entrezDT_file, header=T, stringsAsFactors=F)
entrez2symbDT <- entrez2symbDT[,c("entrezID", "symbol")]
colnames(entrez2symbDT) <- c("entrezID", "geneName")
entrez2symbDT$entrezID <- as.character(entrez2symbDT$entrezID)

######################################################################
samp1 <- eval(parse(text = load(paste0(setDir, "/", sample1_file))))
samp2 <- eval(parse(text = load(paste0(setDir, "/", sample2_file))))

# UPDATE SELECT THE GENES ACCORDING TO THE SETTINGS PREPARED IN 0_PREPGENEDATA
rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_rnaseqDT.Rdata"))))
geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))

rnaseqDT <- rnaseqDT[names(geneList),]    
stopifnot(all(rownames(rnaseqDT) == names(geneList)))
stopifnot(!any(duplicated(names(geneList))))

DE_topTable <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_topTable.Rdata"))))
stopifnot(!any(duplicated(names(geneList))))

gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% geneList,]

DE_topTable$genes <- unlist(sapply(DE_topTable$genes, function(x) geneList[x]))
rownames(DE_topTable) <- NULL
# ! duplicated row.names are not allowed !
dupEntrez <- unique(geneList[duplicated(geneList)])
geneList <- geneList[! geneList %in% dupEntrez]
stopifnot(!any(duplicated(geneList)))

DE_topTable <- DE_topTable[!DE_topTable$genes %in% dupEntrez,]
stopifnot(!any(duplicated(DE_topTable$genes)))

rnaseqDT <- rnaseqDT[which(rownames(rnaseqDT) %in% names(geneList)),]

stopifnot(all(rownames(rnaseqDT) == names(geneList)))
stopifnot(is.numeric(rnaseqDT[1,1]))
log2_rnaseqDT <- log2(rnaseqDT + 0.0001)
meanExpr <- rowMeans(log2_rnaseqDT, na.rm=T)
names(meanExpr) <- unlist(sapply(names(meanExpr), function(x) geneList[x]))

meanTADlogFC <- eval(parse(text = load(paste0(pipOutFold, "/", script3_name, "/", "all_meanLogFC_TAD.Rdata"))))

# retrieve which condition is cond1, i.e. the one that is more expressed when logFC is positive
geneHighestLogFC <- names(geneList[geneList == DE_topTable$genes[which.max(DE_topTable$logFC)] ])

samp1_vect <- rnaseqDT[geneHighestLogFC,samp1, drop=F]
stopifnot(dim(samp1_vect) == c(1,length(samp1)))
samp2_vect <- rnaseqDT[geneHighestLogFC,samp2, drop=F]
stopifnot(dim(samp2_vect) == c(1,length(samp2)))

plot_cond1 <- ifelse(as.numeric(rowMeans(samp1_vect, na.rm = T)) > as.numeric(rowMeans(samp2_vect, na.rm = T)), cond1, cond2)
plot_cond2 <- ifelse(as.numeric(rowMeans(samp1_vect, na.rm = T)) > as.numeric(rowMeans(samp2_vect, na.rm = T)), cond2, cond1)


foo <- foreach(geneSymb = all_geneSymbol) %dopar% {
  # retrieve the corresponding entrezID
  curr_entrez <- as.character(na.omit(symbDT$entrezID[symbDT$symbol == geneSymb]))
  if(length(curr_entrez) == 0) {
    cat(paste0("... > For gene:\t", geneSymb, "\t - found no corresponding entrezID, skip\n"))  
    return(NULL)
  }
  if(length(curr_entrez) > 1) {
    cat(paste0("... > For gene:\t", geneSymb, "\t - found multiple corresponding entrezID, skip\n"))  
    return(NULL)
  }
  # retrieve the corresponding TAD
  curr_tad <- gene2tadDT$region[gene2tadDT$entrezID == curr_entrez]
  if(length(curr_tad) == 0) {
    cat(paste0("... > For gene:\t", geneSymb, "\t - found no TAD, skip\n"))  
    return(NULL)
  }
  cat(paste0("... > For gene:\t", geneSymb, "\t - found:\t", curr_entrez, "\t", curr_tad, " (rank: ", all_emp_pval_combined[curr_tad], ")\n"))  
  
  p <- plot_lolliTAD(TAD_to_plot = curr_tad,
                g2t_table = gene2tadDT,
                id2name_table=entrez2symbDT, 
                DE_table = DE_topTable,
                meanExpr_vect = meanExpr, 
                geneList = geneList,
                textLeft =  meanTADlogFC[curr_tad] > 0,
                orderBy = "logFC", 
                cond1=plot_cond1, cond2=plot_cond2)
  
  outFile <- paste0(curr_outFold, "/plot_gene_TAD_", sub("-", "", geneSymb), "_", curr_entrez, "_", curr_tad, ".", plotType )
  ggsave(plot=p, filename =  outFile, height = plotHeight, width = plotWidth)
  cat(paste0("... written: ", outFile, "\n"))
}
