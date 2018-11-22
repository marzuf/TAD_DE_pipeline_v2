########################################################
### set paths to main files needed accross the workflow + general settings
########################################################

historyDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/data/gene_history_reordered.txt")

# file with mapping from entrez to chromosomic positions
#entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2position/filter_entrez_map.Rdata") # previous also held symbols
# holds only netrezID and positions
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")

# file with symbol synonyms to entrezID
#synoDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/all_entrez_syno_1.Rdata")
# mapping entrezID and all possible symbols
symbolDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_SYMBOL/final_entrez2syno.txt")

# mapping entrezID and ensemblID
ensemblDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_ENSEMBL/final_entrez2ensembl.txt")

# file with coordinates of all regions
TADpos_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_DI_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")    

# file with assignment from entrez to all regions
gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_DI_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 

# match dataset with kind of process
process_file <- paste0(setDir, "/mnt/ed4/marie/other_datasets/dataset_process.csv")
process_file_with_samp <- paste0(setDir, "/mnt/ed4/marie/other_datasets/dataset_process_with_samp_annotated.csv")

nCpu <- 15

gene2tadAssignMethod <- "startPos"
gene2tadAssignMethod <- "maxOverlap"

########################################################
### which kind of data to prepare and use across pipeline
########################################################

# if TRUE, from scripts 2_ to further, use only the genes that pass the DE analysis threshold for all analyses
useFilterCountData <- TRUE

# use only TADs (and corresponding genes)>= minSize et < maxQuantile and genes belonging to those TADs
useFilterSizeData <- TRUE

# if TRUE, from scripts 2_ to further, only focus on  TAD regions
useTADonly <- TRUE

########################################################
### Gene filter min CPM - for filter data (above) and 1_runGeneDE
########################################################

minCpmRatio <- 20/888 


########################################################
### TAD filter number of genes - for filter data (above) and used in 2_runWilcoxonTAD etc.
########################################################

# keep TAD with >= minNbrGeneTAD
minNbrGeneTAD <- 3

# keep TAD with <= maxQuant Genes
maxQuantGeneTAD <- 0.95


########################################################
### permutation settings - 5_runPermutationsMedian.R        => STEP 5
########################################################

withExprClass <- TRUE

# number of permutations to run
nRandomPermut <- 100

# number of classes of expression used for permutate expression data
permutExprClass <- 5

########################################################
### permutation settings (TAD permutation with fix size)- 5b_runPermutationsRandomTADsFixSize.R        => STEP 5b
########################################################
nRandomPermutFixSize <- 1
#percentTADsizeFix <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4, 1.5)
percentTADsizeFix <- c(0.9,1,1.1)

########################################################
### permutation settings (TAD permutation with Gaussian dist. sizes)- 5c_runPermutationsRandomTADsGaussian.R        => STEP 5c
########################################################
nRandomPermutGaussian <- 10
#percentTADsizeGaussian <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4, 1.5)
percentTADsizeGaussian <- c(0.9,1,1.1)

########################################################
### permutation settings (TAD permutation with shuffle)- 5d_runPermutationsRandomTADsShuffle.R        => STEP 5d
########################################################
nRandomPermutShuffle <- 10

########################################################
### ratio calculation also for the randomization of TADs => STEP 8
########################################################

step8_for_permutGenes <- FALSE
step8_for_randomTADsFix <- TRUE
step8_for_randomTADsGaussian <- TRUE
step8_for_randomTADsShuffle <- TRUE

########################################################
### number of TADs to select for plotting                   => STEP 13
########################################################

nTopLolliPlot <- 20
nTopVennPlot <- 50


########################################################
### for the plot with the area                              => STEP 14h
########################################################
step14_for_randomTADsFix <- TRUE
step14_for_randomTADsGaussian <- TRUE

########################################################
### for the plot with the area                              => STEP 14i
########################################################
step14_for_randomTADsShuffle <- TRUE

########################################################
### which ratios to compute                                 => STEPS 8c, 12c, 14c, 17c
########################################################


# allDown <- c("ratioDown", "FCdown", "meanRatioFC" , "prodConcord", "prodMeanConcord", "prodLogRatioNbr", "prodRatioSum", "prodSignedRatio")
# prodLogRatioNbr is not bounded => not used it anymore
# prodRatioSum is similar to prodSignedRatio
#allDown <- c("ratioDown", "FCdown", "meanConcordRatioFC" , "prodConcord", "prodMeanConcord", "prodSignedRatio")
#allDown <- c("ratioDown", "rescWeighted", "rescWeightedQQ", "prodSignedRatio")
#allDown <- c("ratioDown", "FCdown", "prodSignedRatio")
allDown <- c("ratioDown", "prodSignedRatio")
#allDown <- c("prodSignedRatio")
allDown_limited <- allDown
#allDown <- c("FCdown")
#allDown <- "rescWeightedQQ"



########################################################
### which quantile of the permutations to consider          => STEPS 8c, 12c, 14c, 17c
########################################################


permThresh <- 0.95



########################################################
### which scores to plant ranked                             => STEP 18
########################################################

# toPlotRanked <- c("ratioDown", "FCdown") in run settings NOW



########################################################
### to run without voom                             => STEP 1
########################################################

runWithoutVoom <- c("GSE102073", "GSE71119")

