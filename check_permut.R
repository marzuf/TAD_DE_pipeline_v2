# load("OUTPUT_FOLDER/TCGAbrca_lum_bas/5c_runPermutationsRandomTADsGaussian/permutationsList_gaussian_1.Rdata")

### FIRST CHECK: THE FUNCTION THAT ASSIGNS GENES TO TADs WITH MAX OVERLAP WORKS CORRECTLY
# check: when using true set of TADs should get the same gene2tadDT

myfold <- paste0("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")
source(paste0(myfold, "/", "TAD_DE_utils.R"))
require(dplyr)
require(GenomicRanges)

setDir <- "/media/electron"

# consensus TADs
TADpos_file <- paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/06_12_50kb_MAPQFILTER/consensus/gene2tad/KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1/TopDom/TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1_TADposDT.txt")
tadDT <- read.delim(TADpos_file, col.names=c("chromo", "region", "start", "end"), header=F, stringsAsFactors = F)
head(tadDT)

# as assigned earlier
gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/06_12_50kb_MAPQFILTER/consensus/gene2tad/KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1/TopDom/TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1_gene2tadDT.txt")
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names=c("gene", "chromo", "start", "end", "region"), stringsAsFactors = F)
head(gene2tadDT)

# full list of genes
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/GENE_DATA_MZ/gencode.v19.gene.minimal.chr_notAmbiguous_removeDupMZ_changedHeader.txt")
entrezDT <- read.delim(entrezDT_file, header=T, stringsAsFactors = F)


### PREPARE THE GENES I WILL ASSIGN TO THE RANDOM TADs
# do not take the genes on chrY
gene2assignDT <- entrezDT[entrezDT$entrezID %in% gene2tadDT$gene,]
stopifnot(nrow(gene2tadDT) == nrow(gene2assignDT))
# gene2assignDT <- gene2assignDT[,c("entrezID", "chromo", "start", "end", "strand")]
stopifnot(c("chromo", "start", "end", "entrezID") %in% colnames(gene2assignDT))
if(gene2tadAssignMethod == "startPos") {
  stopifnot("strand" %in% colnames(gene2assignDT))
}

# take only those for which I have expression data !!!
# gene2assignDT <- gene2assignDT[as.character(gene2assignDT$entrezID) %in% DE_topTable$genes,]


# assign with the function
assignDT <- assignGene2TADs(regionDT = tadDT, geneDT = entrezDT, assignMethod = "maxOverlap")

missingAssignGenes <- assignDT$entrezID[!assignDT$entrezID %in% gene2tadDT$gene]
length(missingAssignGenes)

missingG2Tgenes <- gene2tadDT$gene[!gene2tadDT$gene %in% assignDT$entrezID]
length(missingG2Tgenes)

all_DT <- inner_join(gene2tadDT, assignDT, by=c("gene"= "entrezID"))
all(all_DT$region.x == all_DT$region.y)

# the assignment seems to be ok

dt1 <- eval(parse(text=load("/media/electron/mnt/ed4/marie/scripts/EZH2_final_MAPQ/GENE_DATA_MZ/CL/cell.lines.byGene_noDup_MZ.RData")))
dt2 <- eval(parse(text = load( "/media/electron/mnt/ed4/marie/scripts/EZH2_final_MAPQ/GENE_DATA_MZ/human/GSE23501_matrixByGene_noDup_MZ.RData")))

dt1_genes <-  unlist(sapply(rownames(dt1), function(x) unlist(strsplit(x,split="///"))))
dt2_genes <-  unlist(sapply(rownames(dt2), function(x) unlist(strsplit(x,split="///"))))

any(missingAssignGenes %in% dt1_genes)
any(missingAssignGenes %in% dt2_genes)

all(gene2tadDT$gene %in% dt1_genes)
all(gene2tadDT$gene %in% dt2_genes)


##########
# all TAD positions
t2t_shuff <- eval(parse(text = load(paste0(shuffDir, "/", "allChr_randTADdt_shuffle_", i_shuff, ".Rdata"))))
t2t_shuff_TAD <- t2t_shuff[grepl("_TAD", t2t_shuff$region),]
tadsize_shuff <- log10(t2t_shuff_TAD$end-t2t_shuff_TAD$start+1)

txt <- paste0("> Number of TADs shuff. data: ", nrow(t2t_shuff_TAD), "\n")
cat(txt)

assignDT2 <- assignGene2TADs(regionDT = t2t_shuff_TAD, geneDT = entrezDT, assignMethod = "maxOverlap")
# for(i in 2:nrow(t2t_shuff_TAD)){
#   stopifnot(t2t_shuff_TAD$start[i] == t2t_shuff_TAD$end[i-1]+1)
# }

# all gene2tad
g2t_shuff <- eval(parse(text = load(paste0(shuffDir, "/", "all_gene_posDT_shuffle_", i_shuff, ".Rdata"))))
txt <- paste0("> Number of genes shuff ", i_shuff, " data: ", nrow(g2t_shuff), "\n")
cat(txt)

all(g2t_shuff %in% c(dt1_genes, dt2_genes))

missingGenes2 <- assignDT2$entrezID[!assignDT2$entrezID %in% g2t_shuff]