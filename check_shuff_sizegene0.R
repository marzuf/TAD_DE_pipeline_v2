require(foreach)
require(doMC)
require(GenomicRanges)

set.seed(20180102)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

registerDoMC(ifelse(SSHFS, 2, 30))

outFold <- "check_shuff_sizegene0"
system(paste0("mkdir -p ", outFold))

### generate randomTADs using the shuffle randomization
shuffle_sourceFile <- paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/ezh2_utils_fct.R")
shuffle_chromoPartition_v1 <- local({
  source(shuffle_sourceFile, local = TRUE)
  environment(shuffle_chromoPartition_v1) <- .GlobalEnv
  shuffle_chromoPartition_v1
})

assign_sourceFile <-  paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
assignGene2TADs <- local({
  source(assign_sourceFile, local = TRUE)
  environment(assignGene2TADs) <- .GlobalEnv
  assignGene2TADs
})

plot_multiDens <- function(size_list, plotTit="", legTxt=NULL, legPos="topright", my_ylab="density", my_xlab="") {
  dens <- lapply(size_list, function(x) density(na.omit(x)))
  names(dens) <- names(size_list)
  lengthDens <- unlist(lapply(size_list, function(x) length(na.omit(x))))
  plot(NA, xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")), 
       main=plotTit, xlab=my_xlab, ylab=my_ylab)
  foo <- mapply(lines, dens, col=1:length(dens))
  if(is.null(legTxt)){
    # legTxt <- names(dens)
    legTxt <- paste0(names(dens), " (n=", lengthDens, ")")
  }
  legend(legPos, legend=legTxt, fill=1:length(dens), bty='n')
}


minThresh <- 3
maxQuantile <- 0.99

### TAD data

TADpos_file <- paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/06_12_50kb_MAPQFILTER/consensus/gene2tad/KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1/TopDom/TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1_TADposDT.txt")
TADdt <- read.delim(TADpos_file, header=F, col.names=c("chromo", "region", "start", "end"), stringsAsFactors = F)
#chr1    chr1_TAD1       750001  1300000
#chr1    chr1_TAD2       2750001 3650000
#chr1    chr1_TAD3       3650001 4150000

### gene file
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/GENE_DATA_MZ/gencode.v19.gene.minimal.chr_notAmbiguous_removeDupMZ_changedHeader.txt")
entrezDT <- read.delim(entrezDT_file, header=T, stringsAsFactors=F)


### gene2tad file => for check
gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/06_12_50kb_MAPQFILTER/consensus/gene2tad/KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1/TopDom/TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1_gene2tadDT.txt")
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names=c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
#LINC00115       chr1    761586  762902  chr1_TAD1
#FAM41C  chr1    803451  812283  chr1_TAD1
#SAMD11  chr1    860260  879955  chr1_TAD1
#NOC2L   chr1    879584  894689  chr1_TAD1
stopifnot(!any(duplicated(gene2tadDT$entrezID)))
tmp_gene2tadDT <- gene2tadDT[,c("entrezID", "region")]
rownames(tmp_gene2tadDT) <- tmp_gene2tadDT$entrezID

### expression table

rnaseqDT_file <- paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/GENE_DATA_MZ/CL/cell.lines.byGene_noDup_MZ.RData")
exprDT <- eval(parse(text = load(rnaseqDT_file)))
expr_genes <- as.character(unlist(sapply(rownames(exprDT), function(x) unlist(strsplit(x,split="///")))))


### take only the genes for which I have expression
entrezDT <- entrezDT[entrezDT$symbol %in% expr_genes,]

######### assign the genes
g2t_dt <- assignGene2TADs(regionDT=TADdt, geneDT = entrezDT, assignMethod = "maxOverlap")


### check that with the true data, I should have the same gene2tad assignment
innerGenes <- intersect(g2t_dt$entrezID, gene2tadDT$entrezID)
length(innerGenes)

tmp_g2t_dt <- g2t_dt
stopifnot(!any(duplicated(tmp_g2t_dt$entrezID)))
rownames(tmp_g2t_dt) <- tmp_g2t_dt$entrezID

inter_true_g2t_dt <- tmp_gene2tadDT[innerGenes,]
inter_assigned_g2t_dt <- tmp_g2t_dt[innerGenes,]

stopifnot(all(inter_true_g2t_dt$region == inter_assigned_g2t_dt$region))

##################################################################################################################################################################
##################################################################################################################################################################
# FROM NOW, WORK ONLY WITH TADs
gene2tadDT <- gene2tadDT[grep("_TAD", gene2tadDT$region),]
obs_genes_TAD <- setNames(as.numeric(table(gene2tadDT$region)), names(table(gene2tadDT$region)))
head(obs_genes_TAD)

maxThresh <- as.numeric(quantile(obs_genes_TAD, probs=maxQuantile))

sizeFilter_obs <- obs_genes_TAD[obs_genes_TAD >= minThresh & obs_genes_TAD <= maxThresh]
length(sizeFilter_obs)

maxEndDT <- aggregate(end~chromo, data=TADdt, FUN=max, na.rm=TRUE)
TADdt <- TADdt[grep("_TAD", TADdt$region),]
domainDT <- TADdt[,c("chromo", "start", "end")]

all_chromo <- intersect(domainDT$chromo, gene2tadDT$chromo)

################################################################################# NOW SHUFFLE AND ASSIGN
nPerm=1000

outFile <- paste0(outFold, "/", "shuff_genesByTAD_withSize.Rdata")

# shuff_genesByTAD_withSize <- foreach(i_perm=1:nPerm) %dopar% {
#   allChr_randTADdt <- foreach(chromo = all_chromo, .combine="rbind") %do% {
#     cat(paste0("... perm ", i_perm, " - ", chromo, "\n"))
#     chrEnd <- maxEndDT$end[maxEndDT$chromo == chromo]
#     # select the initial seet of TADs for this chromosome
#     chromo_domainDT <- domainDT[domainDT$chromo == chromo,]
#     randTADdt <- shuffle_chromoPartition_v1(domainDT=chromo_domainDT, chrSize = chrEnd , preservePattern=FALSE)
#     randTADdt$region <- paste0(chromo, "_TAD", 1:nrow(randTADdt))
#     # the area of the chromo covered by TADs should be the same before/after shuffling !
#     stopifnot(abs(sum(chromo_domainDT$end-chromo_domainDT$start) - sum(randTADdt$end-randTADdt$start)) < 1e-10)
#     # ensure not overlapping
#     if(nrow(randTADdt) > 1) {
#       for(i in 2:nrow(randTADdt))
#         stopifnot(randTADdt$start[i] > randTADdt$end[i-1])
#     }
#     # ensure all starts smaller than ends
#     stopifnot(randTADdt$start < randTADdt$end)
#     randTADdt
#   } # end building new
# 
#   shuff_g2t_dt <- assignGene2TADs(regionDT=allChr_randTADdt, geneDT = entrezDT, assignMethod = "maxOverlap")
#   shuff_g2t_dt$region <- factor(as.character(shuff_g2t_dt$region), levels = as.character(allChr_randTADdt$region))
#   stopifnot(!any(is.na(shuff_g2t_dt$region)))
#   
#   tmp <- setNames(as.numeric(table(shuff_g2t_dt$region)), names(table(shuff_g2t_dt$region)))
#   allChr_randTADdt$size <- allChr_randTADdt$end-allChr_randTADdt$start+1
#   stopifnot(!any(duplicated(allChr_randTADdt$region)))
#   stopifnot(all(names(tmp) %in% allChr_randTADdt$region))
#   rownames(allChr_randTADdt) <- allChr_randTADdt$region
#   data.frame(region = names(tmp),
#              nbrGenes = tmp,
#              size = allChr_randTADdt[names(tmp),]$size,
#              stringsAsFactors = F)
#   
# }
# save(shuff_genesByTAD_withSize, file = outFile)
# cat(paste0("... written: ", outFile, "\n"))


load(paste0(outFold, "/", "shuff_genesByTAD_withSize.Rdata"))

n_TADs <- sapply(shuff_genesByTAD_withSize,  length)





################################################################################################################################################################## 
###################################################################################################### size of the regions without genes 
################################################################################################################################################################## 

TADpos_file <- paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/06_12_50kb_MAPQFILTER/consensus/gene2tad/KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1/TopDom/TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1_TADposDT.txt")
TADdt <- read.delim(TADpos_file, header=F, col.names=c("chromo", "region", "start", "end"), stringsAsFactors = F)
#chr1    chr1_TAD1       750001  1300000
#chr1    chr1_TAD2       2750001 3650000
#chr1    chr1_TAD3       3650001 4150000
TADdt <- TADdt[grepl("_TAD", TADdt$region),]
TADdt$size <- TADdt$end - TADdt$start + 1

### gene file
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/GENE_DATA_MZ/gencode.v19.gene.minimal.chr_notAmbiguous_removeDupMZ_changedHeader.txt")
entrezDT <- read.delim(entrezDT_file, header=T, stringsAsFactors=F)


### gene2tad file => for check
gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/06_12_50kb_MAPQFILTER/consensus/gene2tad/KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1/TopDom/TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1_gene2tadDT.txt")
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names=c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
#LINC00115       chr1    761586  762902  chr1_TAD1
gene2tadDT <- gene2tadDT[grep("_TAD", gene2tadDT$region),]
gene2tadDT$region <- factor(as.character(gene2tadDT$region), levels = as.character(TADdt$region))
obs_genes_TAD <- setNames(as.numeric(table(gene2tadDT$region)), names(table(gene2tadDT$region)))
head(obs_genes_TAD)

stopifnot(!any(duplicated(TADdt$region)))
stopifnot(all(names(obs_genes_TAD) %in% TADdt$region))
rownames(TADdt) <- TADdt$region

obs_regions_size <- data.frame(region = names(obs_genes_TAD),
           nbrGenes = obs_genes_TAD,
           size = TADdt[names(obs_genes_TAD),]$size,
           stringsAsFactors = F)

################################################################################################################################################################## 
obs_0_length_DT <- obs_regions_size[obs_regions_size$nbrGenes == 0,]
obs_0_nbr <- sum(obs_regions_size$nbrGenes == 0)

load(paste0(outFold, "/", "shuff_genesByTAD_withSize.Rdata"))


shuff_0_length_list <- lapply(shuff_genesByTAD_withSize, function(x) {
  x$size[x$nbrGenes == 0]
})

shuff_0_nbr_list <- lapply(shuff_genesByTAD_withSize, function(x) {
  sum(x$nbrGenes == 0)
})


################################################################################################################################################################## 
###################################################################################################### number of the regions without genes 
################################################################################################################################################################## 


outFile <- paste0(outFold, "/", "nbrTADs_0genes.png")
png(outFile, width=600, height=600)
plot_multiDens(list(shuff = unlist(shuff_0_nbr_list)), 
               my_xlab = "nbr TADs",
               plotTit = "Nbr of regions with 0 genes")
legend("topleft", legend = paste0("nbr obs.: ", obs_0_nbr), bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


################################################################################################################################################################## 
###################################################################################################### size of the regions without genes 



outFile <- paste0(outFold, "/", "sizeBpTADs_0genes.png")
png(outFile, width=600, height=600)
plot_multiDens(list(obs = obs_0_length_DT$size, shuff = unlist(shuff_0_length_list)), 
               my_xlab = "size of the TAD (bp)",
               plotTit = "size of the regions with 0 genes")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



outFile <- paste0(outFold, "/", "sizeBpTADs_0genes_log10.png")
png(outFile, width=600, height=600)
plot_multiDens(list(obs = log10(obs_0_length_DT$size), shuff = log10(unlist(shuff_0_length_list))), 
               my_xlab = "size of the TADs (bp) [log10]",
               plotTit = "size of the regions with 0 genes")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


