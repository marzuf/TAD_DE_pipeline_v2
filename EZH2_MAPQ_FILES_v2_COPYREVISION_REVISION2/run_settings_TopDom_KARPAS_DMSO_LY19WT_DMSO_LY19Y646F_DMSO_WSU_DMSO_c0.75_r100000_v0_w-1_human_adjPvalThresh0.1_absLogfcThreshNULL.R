# > file written: Wed, 20 Jun 2018 17:50:40 +0200 
# *** EXPECTED FORMAT:
# full path (starting with /mnt/...)
# colnames = samplesID
# rownames = geneID

# CONSENSUS TADs HAVE BEEN CALLED WITH FOLLOWING PARAMETERS:
# -c 0.75 -r100000 -v0 -w-1

# path to output folder:
pipOutFold <- "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/OUTPUT_FOLDER_EZH2_MAPQ_v2_REVISION_adjPvalThresh0.1_absLogfcThreshNULL_human/TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1"

# EZH2 data are microarray data:
microarray <- TRUE    # UPDATE 20.06.2018 FOR ICGC DATA !!!, stored in variable
normRNAseq <- FALSE    # UPDATE 20.06.2018 FOR ICGC DATA !!!, added variable, ICGC already normalized (no need for FPKM)

# OVERWRITE THE DEFAULT SETTINGS FOR INPUT FILES # there are the same as for EZH2
TADpos_file <- paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/06_12_50kb_MAPQFILTER/consensus/gene2tad/KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1/TopDom/TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1_TADposDT.txt")
#chr1    chr1_TAD1       750001  1300000
#chr1    chr1_TAD2       2750001 3650000
#chr1    chr1_TAD3       3650001 4150000

gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/06_12_50kb_MAPQFILTER/consensus/gene2tad/KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1/TopDom/TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1_gene2tadDT.txt")
#LINC00115       chr1    761586  762902  chr1_TAD1
#FAM41C  chr1    803451  812283  chr1_TAD1
#SAMD11  chr1    860260  879955  chr1_TAD1
#NOC2L   chr1    879584  894689  chr1_TAD1

# entrezDT_file <-  paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/06_12_50kb_MAPQFILTER/consensus/gene2tad/KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1/TopDom/TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1_entrezDT.txt")
#entrezID        symbol
#LINC00115       LINC00115
#FAM41C  FAM41C
#SAMD11  SAMD11
#NOC2L   NOC2L
            # CHANGED V2 !!!!
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/GENE_DATA_MZ/gencode.v19.gene.minimal.chr_notAmbiguous_removeDupMZ_changedHeader.txt")


# overwrite main_settings.R: nCpu <- 25
nCpu <- 15

# *************************************************************************************************************************
# ************************************ SETTINGS FOR 0_prepGeneData
# *************************************************************************************************************************

rnaseqDT_file <- "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/GENE_DATA_MZ/human/GSE23501_matrixByGene_noDup_MZ.RData"
my_sep <- "\t"

# EZH2 INPUT FORMAT IS Rdata
inRdata <- TRUE

# EZH2 FORMAT: THIS IS A HACK, SET entrezID
geneID_format <- "entrezID"
stopifnot(geneID_format %in% c("ensemblID", "entrezID", "geneSymbol"))

# EZH2: GENE ID ARE IN ROW NAMES
geneID_loc <- "rn"
stopifnot(geneID_loc == "rn" | is.numeric(geneID_loc))

# EZH2: not used for EZH2 (use match between rownames and gene ID, but rownames cannot be used in EZH2 data)
removeDupGeneID <- TRUE

# FOR TopDom: do not set max quantile limit (1), for DI set to 0.95
maxQuantGeneTAD <- 0.99


# update 18.06.2018: ADDED threshold for adjPval and absLogfc
gene_adjPvalThresh <- 0.1
gene_absLogfcThresh <- NULL

# *************************************************************************************************************************
# ************************************ SETTINGS FOR 1_runGeneDE
# *************************************************************************************************************************

# for the paired scenario, need a design file (3 columns, with header name, 1st col=sample IDs; 2nd col=pair IDs; 3dcol=conditions)
designFile <- NULL

# labels for conditions
cond1 <- "wt"
cond2 <- "mut"

# path to sampleID for each condition - should be Rdata
sample1_file <- "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/gene_expression/human/wt_ID.Rdata"
sample2_file <- "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/gene_expression/human/mut_ID.Rdata"




# *************************************************************************************************************************
# ************************************ SETTINGS FOR PERMUTATIONS (5#_, 8c_)
# *************************************************************************************************************************

# number of permutations
nRandomPermut <- 10000

gene2tadAssignMethod <- "maxOverlap"

nRandomPermutShuffle <- 10000

step8_for_permutGenes <- TRUE
step8_for_randomTADsFix <- FALSE
step8_for_randomTADsGaussian <- FALSE
step8_for_randomTADsShuffle <- FALSE   # set to FALSE for REVISION 19.06.2018

step14_for_randomTADsShuffle <- FALSE  # set to FALSE for REVISION 19.06.2018
            

