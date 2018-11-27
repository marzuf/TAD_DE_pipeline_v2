require(dplyr)
require(foreach)

setDir <- ""
# setDir <- "/media/electron"


load(paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/gene_expression/RNAseq/tads.RData"))
tads$genes <- as.character(tads$genes)
tads_genes <- unlist(sapply(tads$genes, function(x) strsplit(x=x, split=",")))

load(paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/gene_expression/RNAseq/ly19.notInf.RData"))
ly19.notInf$gene_short_name <- as.character(ly19.notInf$gene_short_name)


fpkmFile <- paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/gene_expression/RNAseq/LYWT_LY19Y641F_CuffDiff_FPKMs_fc.txt")
stopifnot(file.exists(fpkmFile))
fpkmDT <- read.delim(fpkmFile, header=T, stringsAsFactors = FALSE)
fpkmDT <- fpkmDT[,c("gene_short_name", "chr", "start", "end","log2_fc_y_vs_x")]
tmpDT <- inner_join(fpkmDT, ly19.notInf, by=c("gene_short_name"))
stopifnot(all(tmpDT$log2_fc_y_vs_x == tmpDT$log2_fc))


genesFC_not_in_TADs <- ly19.notInf$gene_short_name[!ly19.notInf$gene_short_name %in% tads_genes]

# check the observed fold change
# split the genes from the TADs
tads_genes_DT <- foreach(i = 1:nrow(tads), .combine='rbind') %dopar% {
  data.frame(gene_short_name = unlist(strsplit(x=tads$genes[i], split = ",")),
             tad = as.character(tads$region[i]),
             stringsAsFactors = F)
}

exprTAD_dt <- inner_join(ly19.notInf[,c("gene_short_name", "log2_fc")], tads_genes_DT, by="gene_short_name")

fc_DT <- aggregate(log2_fc ~ tad, data=exprTAD_dt, function(x) sum(x<0)/length(x))
obs_ratio <- setNames(fc_DT$log2_fc, fc_DT$tad)
obs_ratio <- obs_ratio[order(names(obs_ratio))]

ly19_tad_genes <- ly19.notInf$gene_short_name[ly19.notInf$gene_short_name %in% tads_genes_DT$gene_short_name]

nGenes_TAD <- setNames(as.numeric(table(exprTAD_dt$tad)), names(as.numeric(table(exprTAD_dt$tad))))
stopifnot(length(nGenes_TAD) == 1809)
stopifnot(length(nGenes_TAD[nGenes_TAD >=3 & nGenes_TAD <=9]) == 900)
# => OK, values retrieved from wave plots

##########################################################################
########################################################################## 0 data
##########################################################################

load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/OUTPUT_FOLDER_EZH2_MAPQ_v2_RNAseq/RNAseq_CL/5d_runPermutationsRandomTADsShuffle/all_gene_posDT_shuffle_1.Rdata"))
head(all_gene_posDT)
## => problem from the shuffling: only genes from TADs were taken !
all(all_gene_posDT$entrezID %in% tads_genes)
# TRUE
stopifnot(all(all_gene_posDT$entrezID %in% tads_genes))
# ... but should have taken only genes for which I have FC values
stopifnot(all(all_gene_posDT$entrezID %in% ly19.notInf$gene_short_name))

load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/OUTPUT_FOLDER_EZH2_MAPQ_v2_RNAseq/RNAseq_CL/5d_runPermutationsRandomTADsShuffle/permutationsList_shuffle.Rdata"))
foo <- lapply(permutationsList, function(x) {
  stopifnot(all(names(x) %in% tads_genes))
  # ... but should have taken only genes for which I have FC values
  stopifnot(all(names(x) %in% ly19.notInf$gene_short_name))
  })

## => I kept DE (FC) data only from the genes in TADs...
load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/OUTPUT_FOLDER_EZH2_MAPQ_v2_RNAseq/RNAseq_CL/1_runGeneDE/DE_rnaseqDT.Rdata"))
all(rownames(DE_rnaseqDT) %in% tads_genes)
# FALSE

# check I obtain the same ratio down as from starting from ly19dt
load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/OUTPUT_FOLDER_EZH2_MAPQ_v2_RNAseq/RNAseq_CL/8c_runAllDown/all_obs_ratioDown.Rdata"))
stopifnot(length(all_obs_ratio) == length(obs_ratio))
all_obs_ratio <- all_obs_ratio[order(names(all_obs_ratio))]
stopifnot(names(all_obs_ratio) == names(obs_ratio))
stopifnot(all_obs_ratio == obs_ratio)

# check I use all the genes for gene permutation
load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/OUTPUT_FOLDER_EZH2_MAPQ_v2_RNAseq/RNAseq_CL/5_runPermutationsMedian/permutationsDT.Rdata"))
permutationsDT[1:5,1:5]
stopifnot(nrow(permutationsDT) == length(ly19_tad_genes))
stopifnot(sort(rownames(permutationsDT)) == sort(ly19_tad_genes))

##########################################################################
########################################################################## vpip data
##########################################################################

load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/OUTPUT_FOLDER_EZH2_MAPQ_v2_RNAseq/RNAseq_CL_vpip/5d_runPermutationsRandomTADsShuffle/all_gene_posDT_shuffle_1.Rdata"))
head(all_gene_posDT)
## => problem from the shuffling: only genes from TADs were taken !
all(all_gene_posDT$entrezID %in% tads_genes)
# TRUE
stopifnot(all(all_gene_posDT$entrezID %in% tads_genes))
# ... but should have taken only genes for which I have FC values
stopifnot(all(all_gene_posDT$entrezID %in% ly19.notInf$gene_short_name))

load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/OUTPUT_FOLDER_EZH2_MAPQ_v2_RNAseq/RNAseq_CL_vpip/5d_runPermutationsRandomTADsShuffle/permutationsList_shuffle.Rdata"))
foo <- lapply(permutationsList, function(x) {
  stopifnot(all(names(x) %in% tads_genes))
  # ... but should have taken only genes for which I have FC values
  stopifnot(all(names(x) %in% ly19.notInf$gene_short_name))
})

## => I kept DE (FC) data only from the genes in TADs...
load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/OUTPUT_FOLDER_EZH2_MAPQ_v2_RNAseq/RNAseq_CL_vpip/1_runGeneDE/DE_rnaseqDT.Rdata"))
all(rownames(DE_rnaseqDT) %in% tads_genes)
# TRUE

# check I obtain the same ratio down as from starting from ly19dt
load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/OUTPUT_FOLDER_EZH2_MAPQ_v2_RNAseq/RNAseq_CL_vpip/8c_runAllDown/all_obs_ratioDown.Rdata"))
stopifnot(length(all_obs_ratio) == length(obs_ratio))
all_obs_ratio <- all_obs_ratio[order(names(all_obs_ratio))]
stopifnot(names(all_obs_ratio) == names(obs_ratio))
stopifnot(all_obs_ratio == obs_ratio)


# check I use all the genes for gene permutation
load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/OUTPUT_FOLDER_EZH2_MAPQ_v2_RNAseq/RNAseq_CL_vpip/5_runPermutationsMedian/permutationsDT.Rdata"))
permutationsDT[1:5,1:5]
stopifnot(nrow(permutationsDT) == length(ly19_tad_genes))
stopifnot(sort(rownames(permutationsDT)) == sort(ly19_tad_genes))

##########################################################################
########################################################################## v4 data
##########################################################################


# check I obtain the same ratio down as from starting from ly19dt
load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/OUTPUT_FOLDER_EZH2_MAPQ_v2_RNAseq/RNAseq_CL_v4/5d_runPermutationsRandomTADsShuffle/all_gene_posDT_shuffle_1.Rdata"))
head(all_gene_posDT)
all(all_gene_posDT$entrezID %in% tads_genes)
#  FALSE
stopifnot(!all(all_gene_posDT$entrezID %in% tads_genes))
# ... but should have taken only genes for which I have FC values
stopifnot(all(all_gene_posDT$entrezID %in% ly19.notInf$gene_short_name))

load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/OUTPUT_FOLDER_EZH2_MAPQ_v2_RNAseq/RNAseq_CL_v4/5d_runPermutationsRandomTADsShuffle/permutationsList_shuffle.Rdata"))
foo <- lapply(permutationsList, function(x) {
  stopifnot(!all(names(x) %in% tads_genes))
  # ... but should have taken only genes for which I have FC values
  stopifnot(all(names(x) %in% ly19.notInf$gene_short_name))
})

## => I kept DE (FC) data only from the genes in TADs...
load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/OUTPUT_FOLDER_EZH2_MAPQ_v2_RNAseq/RNAseq_CL_v4/1_runGeneDE/DE_rnaseqDT.Rdata"))
all(rownames(DE_rnaseqDT) %in% tads_genes)
# FALSE
stopifnot(!all(rownames(DE_rnaseqDT) %in% tads_genes))

load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/OUTPUT_FOLDER_EZH2_MAPQ_v2_RNAseq/RNAseq_CL_v4/5_runPermutationsMedian/permutationsDT.Rdata"))
permutationsDT[1:5,1:5]
stopifnot(nrow(permutationsDT) == length(ly19_tad_genes))
stopifnot(sort(rownames(permutationsDT)) == sort(ly19_tad_genes))

##########################################################################
########################################################################## v4 sizeLim
##########################################################################


# check I obtain the same ratio down as from starting from ly19dt
load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/OUTPUT_FOLDER_EZH2_MAPQ_v2_RNAseq/RNAseq_CL_v4_sizeLim/5d_runPermutationsRandomTADsShuffle/all_gene_posDT_shuffle_1.Rdata"))
head(all_gene_posDT)
all(all_gene_posDT$entrezID %in% tads_genes)
#  FALSE
stopifnot(!all(all_gene_posDT$entrezID %in% tads_genes))
# ... but should have taken only genes for which I have FC values
stopifnot(all(all_gene_posDT$entrezID %in% ly19.notInf$gene_short_name))

load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/OUTPUT_FOLDER_EZH2_MAPQ_v2_RNAseq/RNAseq_CL_v4_sizeLim/5d_runPermutationsRandomTADsShuffle/permutationsList_shuffle.Rdata"))
foo <- lapply(permutationsList, function(x) {
  stopifnot(!all(names(x) %in% tads_genes))
  # ... but should have taken only genes for which I have FC values
  stopifnot(all(names(x) %in% ly19.notInf$gene_short_name))
})

## => I kept DE (FC) data only from the genes in TADs...
load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/OUTPUT_FOLDER_EZH2_MAPQ_v2_RNAseq/RNAseq_CL_v4_sizeLim/1_runGeneDE/DE_rnaseqDT.Rdata"))
all(rownames(DE_rnaseqDT) %in% tads_genes)
# FALSE
stopifnot(!all(rownames(DE_rnaseqDT) %in% tads_genes))

load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/OUTPUT_FOLDER_EZH2_MAPQ_v2_RNAseq/RNAseq_CL_v4_sizeLim/5_runPermutationsMedian/permutationsDT.Rdata"))
permutationsDT[1:5,1:5]
stopifnot(nrow(permutationsDT) == length(ly19_tad_genes))
stopifnot(sort(rownames(permutationsDT)) == sort(ly19_tad_genes))