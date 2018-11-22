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

library(foreach)
library(doMC)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

registerDoMC(ifelse(SSHFS, 2, 20))


all_datasets <- c("CL", "human", "GSK", "EPZ")
all_shuff <- c(1:10)

plotType <- "png"
myHeight <- 600
myWidth <- 600

outFold <- paste0("CMP_trueTADs_shuffTADs")
system(paste0("mkdir -p ", outFold))

# dataset <- "CL"

activeLogfile <- paste0(outFold, "/", "activity_logFile.txt")
system(paste0("rm -f ", activeLogfile))

for(dataset in all_datasets){
  
  prepDir <- paste0("OUTPUT_FOLDER_EZH2_MAPQ_v2_", dataset, "/TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1/0_prepGeneData")
  shuffDir <- paste0("OUTPUT_FOLDER_EZH2_MAPQ_v2_", dataset, "/TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1/5d_runPermutationsRandomTADsShuffle")
  deDir <-  paste0("OUTPUT_FOLDER_EZH2_MAPQ_v2_", dataset, "/TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1//1_runGeneDE")
  
  rangeSize <- eval(parse(text = load(paste0(prepDir, "/", "rangeTADgenes.Rdata"))))
  
  settingF <- paste0("EZH2_MAPQ_FILES_v2/run_settings_TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1_", dataset, ".R")
  source(settingF)
  
  
  ratioShuffFile <-  paste0("OUTPUT_FOLDER_EZH2_MAPQ_v2_", dataset, "/TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1/8c_runAllDown/ratioDown_randomShuffleList.Rdata")
  ratioObsFile <-  paste0("OUTPUT_FOLDER_EZH2_MAPQ_v2_", dataset, "/TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1/8c_runAllDown/all_obs_ratioDown.Rdata")
  observedRatio <-  eval(parse(text =  load(ratioObsFile)))
  ratioList <- eval(parse(text =  load(ratioShuffFile)))
  
  ratioPermutFile <-  paste0("OUTPUT_FOLDER_EZH2_MAPQ_v2_", dataset, "/TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1/8c_runAllDown/ratioDown_permDT.Rdata")
  ratioPermut <- eval(parse(text =  load(ratioPermutFile)))
  
  all_ratios <- c(0.1,0.2)
  
  ### FOR THE PERMUTATION = TAD SHUFFLING
  
  for(ratioActive in all_ratios) {
    # inactive = ratioDown >= 0.9 
    inactive_obs <- observedRatio[observedRatio >= (1-ratioActive)]
    nInactive_obs <- length(inactive_obs)
    nInactive_shuff <- sapply(ratioList, function(x) length(x[x >= (1-ratioActive)]))
    
    outFile <- paste0(outFold, "/", "inactive_", gsub("\\.", "", ratioActive), "_", dataset, "_TADshuff.", plotType)
    do.call(plotType, list(outFile, width=myWidth, height=myHeight))
    plot_multiDens(list(nbrInactive_shuff=nInactive_shuff), plotTit = paste0(dataset, " - nbr of TADs with ratioDown >= ", (1-ratioActive)))
    abline(v=length(nInactive_obs))
    
    legTxt <- paste0("obs. # inactive TADs= ", nInactive_obs,  
                     "\n(range shuffl.: ", paste0(range(nInactive_shuff), collapse="-"), ")",
                     "\n(mean shuffl.: ", round(mean(nInactive_shuff), 2), ")",
                     "\n(sd shuffl.: ", round(sd(nInactive_shuff), 2), 
                     ")")
    
    cat(paste0("> ", dataset, " - ", ratioActive, " - permutations: TAD shuffling\n"), file = activeLogfile, append=T)
    cat(paste0(legTxt, "\n"), file = activeLogfile, append=T)
    
    legend("topleft", bty="n", legend = legTxt)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
     
   # active = ratioDown <= 0.1
   active_obs <- observedRatio[observedRatio <= ratioActive]
   nActive_obs <- length(active_obs)
   nActive_shuff <- sapply(ratioList, function(x) length(x[x <= ratioActive]))
   
   outFile <- paste0(outFold, "/", "active_", gsub("\\.", "", ratioActive), "_", dataset, "_TADshuff.", plotType)
   do.call(plotType, list(outFile, width=myWidth, height=myHeight))
   plot_multiDens(list(nbrActive_shuff=nActive_shuff), plotTit = paste0(dataset, "- nbr of TADs with ratioDown <= ", ratioActive))
   abline(v=length(nActive_obs))
   
   legTxt <- paste0("obs. # active TADs= ", nActive_obs, 
                    "\n(range shuffl.: ", paste0(range(nActive_shuff), collapse="-"), ")",
                    "\n(mean shuffl.: ", round(mean(nActive_shuff), 2), ")",
                    "\n(sd shuffl.: ", round(sd(nActive_shuff), 2),
                    ")")
   cat(paste0(legTxt, "\n"), file = activeLogfile, append=T)
   
   legend("topleft", bty="n", legend = legTxt)
   
   foo <- dev.off()
   cat(paste0("... written: ", outFile, "\n"))
  }
 
  
  
  ### FOR THE PERMUTATION = TAD SHUFFLING
  
  for(ratioActive in all_ratios) {
    # inactive = ratioDown >= 0.9 
    inactive_obs <- observedRatio[observedRatio >= (1-ratioActive)]
    nInactive_obs <- length(inactive_obs)
    nInactive_permut <- apply(ratioPermut, 2, function(x) length(x[x >= (1-ratioActive)]))
    
    outFile <- paste0(outFold, "/", "inactive_", gsub("\\.", "", ratioActive), "_", dataset, "_genePermut.", plotType)
    do.call(plotType, list(outFile, width=myWidth, height=myHeight))
    plot_multiDens(list(nbrInactive_permut=nInactive_permut), plotTit = paste0(dataset, " - nbr of TADs with ratioDown >= ", (1-ratioActive)))
    abline(v=length(nInactive_obs))
    
    legTxt <- paste0("obs. # inactive TADs= ", nInactive_obs,  
                     "\n(range permut.: ", paste0(range(nInactive_permut), collapse="-"), ")",
                     "\n(mean permut.: ", round(mean(nInactive_permut), 2), ")",
                     "\n(sd permut.: ", round(sd(nInactive_permut), 2), 
                     ")")
    
    cat(paste0("> ", dataset, " - ", ratioActive, " - permutations: gene permut.\n"), file = activeLogfile, append=T)
    cat(paste0(legTxt, "\n"), file = activeLogfile, append=T)
    
    legend("topleft", bty="n", legend = legTxt)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    # active = ratioDown <= 0.1
    active_obs <- observedRatio[observedRatio <= ratioActive]
    nActive_obs <- length(active_obs)
    nActive_permut <- apply(ratioPermut, 2, function(x) length(x[x <= ratioActive]))
    
    outFile <- paste0(outFold, "/", "active_", gsub("\\.", "", ratioActive), "_", dataset, "_genePermut.", plotType)
    do.call(plotType, list(outFile, width=myWidth, height=myHeight))
    plot_multiDens(list(nbrActive_permut=nActive_permut), plotTit = paste0(dataset, "- nbr of TADs with ratioDown <= ", ratioActive))
    abline(v=length(nActive_obs))
    
    legTxt <- paste0("obs. # active TADs= ", nActive_obs, 
                     "\n(range permut.: ", paste0(range(nActive_permut), collapse="-"), ")",
                     "\n(mean permut.: ", round(mean(nActive_permut), 2),")",
                     "\n(sd permut.: ", round(sd(nActive_permut), 2),
                     ")")
    cat(paste0(legTxt, "\n"), file = activeLogfile, append=T)
    
    legend("topleft", bty="n", legend = legTxt)
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
  
   
}

stop("---ok\n")

geneList <- eval(parse(text = load(paste0(prepDir, "/", "rna_geneList.Rdata"))))
pipeline_regionList <- eval(parse(text = load(paste0(prepDir, "/", "pipeline_regionList.Rdata"))))
DE_dt <- eval(parse(text = load(paste0(deDir, "/", "DE_topTable.Rdata"))))
stopifnot(all(DE_dt$genes %in% names(geneList)))
DE_dt$genes <- unlist(sapply(DE_dt$genes, function(x) geneList[x]))
stopifnot(!is.na(DE_dt$genes))

for(dataset in all_datasets){
  
  prepDir <- paste0("OUTPUT_FOLDER_EZH2_MAPQ_v2_", dataset, "/TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1/0_prepGeneData")
  shuffDir <- paste0("OUTPUT_FOLDER_EZH2_MAPQ_v2_", dataset, "/TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1/5d_runPermutationsRandomTADsShuffle")
  deDir <-  paste0("OUTPUT_FOLDER_EZH2_MAPQ_v2_", dataset, "/TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1//1_runGeneDE")
  
  rangeSize <- eval(parse(text = load(paste0(prepDir, "/", "rangeTADgenes.Rdata"))))
  
  settingF <- paste0("EZH2_MAPQ_FILES_v2/run_settings_TopDom_KARPAS_DMSO_LY19WT_DMSO_LY19Y646F_DMSO_WSU_DMSO_c0.75_r100000_v0_w-1_", dataset, ".R")
  source(settingF)
  
  geneList <- eval(parse(text = load(paste0(prepDir, "/", "rna_geneList.Rdata"))))
  pipeline_regionList <- eval(parse(text = load(paste0(prepDir, "/", "pipeline_regionList.Rdata"))))
  DE_dt <- eval(parse(text = load(paste0(deDir, "/", "DE_topTable.Rdata"))))
  stopifnot(all(DE_dt$genes %in% names(geneList)))
  DE_dt$genes <- unlist(sapply(DE_dt$genes, function(x) geneList[x]))
  stopifnot(!is.na(DE_dt$genes))
  
  
  # check that there are the same region that I have in the full permutation data
  cat("... load permutation data")
  all_list <- eval(parse(text = load(paste0(shuffDir, "/", "permutationsList_shuffle.Rdata"))))
  
  # number of TADs retained
  all_nbrTADs <- unlist(sapply(all_list, function(x) length(unique(x))))
  
  outFile <- paste0(outFold, "/", "nbr_retained_TADs_", dataset, "_shuffdata.", plotType)
  do.call(plotType, list(outFile, width=myWidth, height=myHeight))
  plot_multiDens(list(nbrTADs_shuff=all_nbrTADs), plotTit ="Nbr retained TADs all shuff")
  mtext("(nbr retained TADs during shufflings)")
  abline(v=length(pipeline_regionList))
  legend("topleft", bty="n", legend = paste0("obs. # = ", length(pipeline_regionList)), lty=1)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  foo <- foreach(i_shuff = all_shuff) %dopar% {
    
    logFile <- paste0(outFold, "/", dataset, "_shuff", i_shuff, "_logFile.txt")
    system(paste0("rm -f ", logFile))
    
    #***************************************************************************
    ######### TRUE DATA
    #***************************************************************************
    
    # all TAD positions
    # t2t_trueF <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_DI_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")
    t2t_trueF <- TADpos_file
    t2t_true <- read.delim(t2t_trueF, header=F, col.names=c("chromo", "region", "start","end"), stringsAsFactors = FALSE)
    head(t2t_true)
    t2t_true_TAD <- t2t_true[grepl("_TAD", t2t_true$region),]
    tadsize_true <- log10(t2t_true_TAD$end-t2t_true_TAD$start+1)
    txt <- paste0("> Number of TADs true data: ", nrow(t2t_true_TAD), "\n")
    cat(txt)
    cat(txt, file=logFile, append=TRUE)
    
    # all gene2tad positions !! for EZH2_ contains only the genes for which I have expression values !!!
    # g2t_trueF <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_DI_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt")
    g2t_trueF <- gene2tadDT_file
    g2t_true <- read.delim(g2t_trueF, header=F, col.names=c("gene", "chromo", "start","end", "region"), stringsAsFactors = FALSE)
    head(g2t_true)
    txt <- paste0("> Number of genes true data: ", nrow(g2t_true), "\n")
    cat(txt)
    cat(txt, file=logFile, append=TRUE)
    
    # take only the TADs
    g2t_true_TAD <- g2t_true[grepl("_TAD", g2t_true$region),]
    txt <- paste0("> Number of TADs with genes assigned true data: ", length(unique(g2t_true_TAD$region)), "\n")
    cat(txt)
    cat(txt, file=logFile, append=TRUE)
    txt <- paste0("> Number of TAD genes true data: ", nrow(g2t_true_TAD), "\n")
    cat(txt)
    cat(txt, file=logFile, append=TRUE)
    # ... how many genes by region:
    nGenesByTAD_true_TAD <- setNames(as.numeric(table(g2t_true_TAD$region)),names(table(g2t_true_TAD$region)) )
    
    # take only the genes with expression data
    g2t_true_TAD_DE <- g2t_true_TAD[g2t_true_TAD$gene %in% DE_dt$genes,]
    txt <- paste0("> Number of TAD genes with DE true data: ", nrow(g2t_true_TAD_DE), "\n")
    cat(txt)
    cat(txt, file=logFile, append=TRUE)
    # ... how many genes by region
    nGenesByTAD_true_TAD_DE <- setNames(as.numeric(table(g2t_true_TAD_DE$region)),names(table(g2t_true_TAD_DE$region)) )
    
    # take only regions that pass size thresh
    filterSizeRegions_true <- names(nGenesByTAD_true_TAD_DE)[nGenesByTAD_true_TAD_DE >= rangeTADgenes[1] & nGenesByTAD_true_TAD_DE <= rangeTADgenes[2]]
    txt <- paste0("> Number of TAD genes with enough and not too many genes true data: ", length(filterSizeRegions_true), "\n")
    cat(txt)
    cat(txt, file=logFile, append=TRUE)
    # ... find out their size
    t2t_true_TAD_sizefilter <- t2t_true_TAD[t2t_true_TAD$region %in% filterSizeRegions_true,]
    tadsize_true_sizeFilter <- log10(t2t_true_TAD_sizefilter$end-t2t_true_TAD_sizefilter$start+1)
    
    g2t_true_TAD_DE_sizeFilter <- g2t_true_TAD_DE[g2t_true_TAD_DE$region %in% filterSizeRegions_true,]
    # ...  how many genes by region
    nGenesByTAD_true_TAD_DE_sizeFilter <- setNames(as.numeric(table(g2t_true_TAD_DE_sizeFilter$region)),names(table(g2t_true_TAD_DE_sizeFilter$region)) )
    
    #***************************************************************************
    ######### SHUFFLED DATA
    #***************************************************************************
    
    # all TAD positions
    t2t_shuff <- eval(parse(text = load(paste0(shuffDir, "/", "allChr_randTADdt_shuffle_", i_shuff, ".Rdata"))))
    t2t_shuff_TAD <- t2t_shuff[grepl("_TAD", t2t_shuff$region),]
    tadsize_shuff <- log10(t2t_shuff_TAD$end-t2t_shuff_TAD$start+1)
    
    txt <- paste0("> Number of TADs shuff. data: ", nrow(t2t_shuff_TAD), "\n")
    cat(txt)
    cat(txt, file=logFile, append=TRUE)
    # for(i in 2:nrow(t2t_shuff_TAD)){
    #   stopifnot(t2t_shuff_TAD$start[i] == t2t_shuff_TAD$end[i-1]+1)
    # }
    
    # all gene2tad
    g2t_shuff <- eval(parse(text = load(paste0(shuffDir, "/", "all_gene_posDT_shuffle_", i_shuff, ".Rdata"))))
    txt <- paste0("> Number of genes shuff ", i_shuff, " data: ", nrow(g2t_shuff), "\n")
    cat(txt)
    cat(txt, file=logFile, append=TRUE)
    
    # take only the TADs
    g2t_shuff_TAD <- g2t_shuff[grepl("_TAD", g2t_shuff$region),]
    txt <- paste0("> Number of TADs with genes assigned shuff", i_shuff, " data: ", length(unique(g2t_shuff_TAD$region)), "\n")
    cat(txt)
    cat(txt, file=logFile, append=TRUE)
    
    txt <- paste0("> Number of TAD genes shuff ", i_shuff, " data: ", nrow(g2t_shuff_TAD), "\n")
    cat(txt)
    cat(txt, file=logFile, append=TRUE)
    # ... how many genes by region
    nGenesByTAD_shuff_TAD <- setNames(as.numeric(table(g2t_shuff_TAD$region)),names(table(g2t_shuff_TAD$region)) )
    
    # take only the genes with expression data
    g2t_shuff_TAD_DE <- g2t_shuff_TAD[g2t_shuff_TAD$entrezID %in% DE_dt$genes,]
    txt <- paste0("> Number of TAD genes with DE shuff ", i_shuff, " data: ", nrow(g2t_shuff_TAD_DE), "\n")
    cat(txt)
    cat(txt, file=logFile, append=TRUE)
    
    # how many genes by region
    nGenesByTAD_shuff_TAD_DE <- setNames(as.numeric(table(g2t_shuff_TAD_DE$region)),names(table(g2t_shuff_TAD_DE$region)) )
    
    # take only regions that pass size thresh
    filterSizeRegions_shuff <- names(nGenesByTAD_shuff_TAD_DE)[nGenesByTAD_shuff_TAD_DE >= rangeTADgenes[1] & nGenesByTAD_shuff_TAD_DE <= rangeTADgenes[2]]
    txt <- paste0("> Number of TAD genes with enough and not too many genes shuff", i_shuff, " data: ", length(filterSizeRegions_shuff), "\n")
    cat(txt)
    cat(txt, file=logFile, append=TRUE)
    
    # ... find out their size
    t2t_shuff_TAD_sizefilter <- t2t_shuff_TAD[t2t_shuff_TAD$region %in% filterSizeRegions_shuff,]
    tadsize_shuff_sizeFilter <- log10(t2t_shuff_TAD_sizefilter$end-t2t_shuff_TAD_sizefilter$start+1)
    
    g2t_shuff_TAD_DE_sizeFilter <- g2t_shuff_TAD_DE[g2t_shuff_TAD_DE$region %in% filterSizeRegions_shuff,]
    # ...  how many genes by region
    nGenesByTAD_shuff_TAD_DE_sizeFilter <- setNames(as.numeric(table(g2t_shuff_TAD_DE_sizeFilter$region)),names(table(g2t_shuff_TAD_DE_sizeFilter$region)) )
    
    stopifnot(length(unique(all_list[[i_shuff]])) == length(filterSizeRegions_shuff))
    
    ##################### # OF GENES BY TAD
    
    outFile <- paste0(outFold, "/", "all_genes_within_TADs_", dataset, "_shuff", i_shuff, ".", plotType)
    do.call(plotType, list(outFile, width=myWidth, height=myHeight))
    plot_multiDens(list(true=nGenesByTAD_true_TAD, shuff=nGenesByTAD_shuff_TAD), plotTit ="Number of genes/TAD", my_xlab = "Number of genes")
    mtext("(all genes within TADs)")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- paste0(outFold, "/", "withDE_genes_within_TADs_", dataset, "_shuff", i_shuff, ".", plotType)
    do.call(plotType, list(outFile, width=myWidth, height=myHeight))
    plot_multiDens(list(true=nGenesByTAD_true_TAD_DE, shuff=nGenesByTAD_shuff_TAD_DE), plotTit ="Number of genes/TAD", my_xlab = "Number of genes")
    mtext("(only genes with DE data)")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- paste0(outFold, "/", "sizeFilterDE_genes_within_TADs_", dataset, "_shuff", i_shuff, ".", plotType)
    do.call(plotType, list(outFile, width=myWidth, height=myHeight))
    plot_multiDens(list(true=nGenesByTAD_true_TAD_DE_sizeFilter, shuff=nGenesByTAD_shuff_TAD_DE_sizeFilter), plotTit ="Number of genes/TAD", my_xlab = "Number of genes")
    mtext("(only genes with DE data+sizeFilter)")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    ##################### TAD SIZE
    outFile <- paste0(outFold, "/", "all_TADs_tadsize_", dataset, "_shuff", i_shuff, ".", plotType)
    do.call(plotType, list(outFile, width=myWidth, height=myHeight))
    plot_multiDens(list(true=tadsize_true, shuff=tadsize_shuff), plotTit ="TAD size (bp) [log10]")
    mtext("(all regions that are TADs)")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- paste0(outFold, "/", "sizeFilter_TADs_tadsize_", dataset, "_shuff", i_shuff, ".", plotType)
    do.call(plotType, list(outFile, width=myWidth, height=myHeight))
    plot_multiDens(list(true=tadsize_true_sizeFilter, shuff=tadsize_shuff_sizeFilter), plotTit ="TAD size (bp) [log10]")
    mtext("(all regions that are TADs, after sizeFilter)")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    # number of TADs
    
    
    #***************************************************************************
    ######### look at some missing TADs
    #***************************************************************************
    
    # missingTADs <- t2t_shuff_TAD$region[!t2t_shuff_TAD$region %in% g2t_shuff_TAD$region]
    # exMissing <- missingTADs[1]
    # exMissingDT <- t2t_shuff_TAD[t2t_shuff_TAD$region == exMissing,]
    # 
    # g2t_chr <- gene2tadDT[gene2tadDT$chromo == as.character(exMissingDT$chromo[1]),]
    # g2t_chr [ (g2t_chr$start >= exMissingDT$start[1] & g2t_chr$start <= exMissingDT$end[1]  ) |
    #   (g2t_chr$end >= exMissingDT$start[1] & g2t_chr$end <= exMissingDT$end[1]  ) ,]
    
    
  } # end iterating over permutations
} #end iterating over datasets
zz



