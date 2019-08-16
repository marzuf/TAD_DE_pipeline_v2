#!/usr/bin/Rscript


cat(paste0("> START ", "1cleanInput",  "\n"))


startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: rna_rnaseqDT.Rdata
# - script0: rna_geneList.Rdata
# - script0: pipeline_geneList.Rdata
################################################################################

################  OUTPUT
# - DE_madnorm_rnaseqDT.Rdata or DE_qqnorm_rnaseqDT.Rdata
# - DE_topTable.Rdata
# - DE_geneList.Rdata
################################################################################

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

script0_name <- "0_prepGeneData"
script_name <- "1_runGeneDE"
stopifnot(file.exists(paste0(pipScriptDir, "/",script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(edgeR, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  
suppressPackageStartupMessages(library(limma, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  


# UPDATE 21.02.2018
# add specification for input data type
stopifnot(exists("inputDataType"))
stopifnot(inputDataType %in% c("raw", "RSEM", "FPKM", "DESeq2", "microarray"))

cat(paste0("> input data type: ", as.character(inputDataType), "\n"))

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

stopifnot(file.exists(paste0(pipOutFold, "/", script0_name, "/rna_rnaseqDT.Rdata")))
rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_rnaseqDT.Rdata"))))
if(ncol(rnaseqDT) >= 5 & nrow(rnaseqDT) >= 5)
    rnaseqDT[1:5,1:5]
initRowNbr <- nrow(rnaseqDT)
stopifnot(is.numeric(rnaseqDT[1,1]))

# ADDED 16.11.2018 to check using other files
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)

# TAKE ONLY THE GENES FOR WHICH I HAVE POSITIONS
stopifnot(file.exists(paste0(pipOutFold, "/", script0_name,  "/", "rna_geneList.Rdata")))
# init_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name,  "/", "rna_geneList.Rdata"))))
rna_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name,  "/", "rna_geneList.Rdata"))))
# => UPDATE: TAKE ONLY THE GENE LIST PREPARED IN 0_prepGeneData ACCORDING TO CURRENT SETTINGS
# UPDATE: compute RNA DE for all the genes, not the filtered ones !
# (in the previous version: DE analysis was done only for the filtered genes; i.e. e.g. those belonging to TADs)
# stopifnot(file.exists(paste0(pipOutFold, "/", script0_name,  "/", "pipeline_geneList.Rdata")))
# rna_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name,  "/", "pipeline_geneList.Rdata"))))
# txt <- paste0(toupper(script_name), "> Start with # genes: ", length(rna_geneList), "/", length(init_geneList), "\n")
# printAndLog(txt, pipLogFile)

txt <- paste0(toupper(script_name), "> Start with # genes: ", length(rna_geneList), "\n")
printAndLog(txt, pipLogFile)

rnaseqDT <- rnaseqDT[which(rownames(rnaseqDT) %in% names(rna_geneList)),]
# alread filtered in 0_geneData
txt <- paste0(toupper(script_name), "> Number of genes available (with position information) for DE analysis: ", nrow(rnaseqDT), "/", initRowNbr, "\n")
printAndLog(txt, pipLogFile)

stopifnot(length(rna_geneList) == nrow(rnaseqDT))

# RUN THE DE ANALYSIS
samp1 <- eval(parse(text=load(paste0(setDir, "/", sample1_file))))
samp2 <- eval(parse(text=load(paste0(setDir, "/", sample2_file))))
# ensure the samples are present in the column names
stopifnot(all(samp1 %in% colnames(rnaseqDT)))
stopifnot(all(samp2 %in% colnames(rnaseqDT)))

rnaseqDT <- rnaseqDT[,c(samp1, samp2)]

totSamples <- length(samp1) + length(samp2)

# FILTER THE EXPRESSION DATA TO MIN CPM FILTER  ################################################################ CPM FOR MICROARRAY OR NOT ????????????????????????????
# apply CPM only on raw counts or RSEM
                    #if(inputDataType == "raw" | inputDataType == "RSEM") {
                    #  cpm_exprDT <- cpm(rnaseqDT)
                    #  txt <- paste0(toupper(script_name), "> NA in cpm_exprDT: ", 
                    #                sum(is.na(cpm_exprDT)), "/", dim(cpm_exprDT)[1]*dim(cpm_exprDT)[2], " (",
                    #                round((sum(is.na(cpm_exprDT))/(dim(cpm_exprDT)[1]*dim(cpm_exprDT)[2]) * 100),2), "%)\n")
                    #  printAndLog(txt, pipLogFile)
                    #  rowsToKeep <- rowSums(cpm_exprDT, na.rm=T) >= (minCpmRatio * ncol(rnaseqDT))
                    #  
                    #  txt <- paste0(toupper(script_name), "> CPM filter, genes to retain: ", sum(rowsToKeep), "/", nrow(rnaseqDT), "\n")
                    #  printAndLog(txt, pipLogFile)
                    #} else if(inputDataType == "FPKM") {
                    #  txt <- paste0(toupper(script_name), "> !!! FPKM filter not applied !!!", "\n")
                    #  printAndLog(txt, pipLogFile)
                    #  cpm_exprDT <- rnaseqDT
                    #  txt <- paste0(toupper(script_name), "> NA in cpm_exprDT: ", 
                    #                sum(is.na(cpm_exprDT)), "/", dim(cpm_exprDT)[1]*dim(cpm_exprDT)[2], " (",
                    #                round((sum(is.na(cpm_exprDT))/(dim(cpm_exprDT)[1]*dim(cpm_exprDT)[2]) * 100),2), "%)\n")
                    #  printAndLog(txt, pipLogFile)
                    #  rowsToKeep <- rowSums(cpm_exprDT, na.rm=T) >= (minCpmRatio * ncol(rnaseqDT))
                    #} else if (inputDataType == "microarray" | inputDataType == "DESeq2") {
                    #  txt <- paste0(toupper(script_name), "> !!! CPM filter not applied !!!", "\n")
                    #  printAndLog(txt, pipLogFile)
                    #  rowsToKeep <- rep(TRUE, nrow(rnaseqDT))
                    #} else {
                    #  stop("ERROR\n")
                    #}
# UPDATE 15.08.2019
  countFilter_rnaseqDT <- rnaseqDT
  if(inputDataType == "raw" | inputDataType == "RSEM") {

    # CHANGE: at least min_counts reads in at least min_sampleRatio of samples
    # for each gene, how many samples have enough read (number)
    genes_nSamples_atLeastMinReads <- apply(countFilter_rnaseqDT, 1, function(x) sum(x >= min_counts)) 
    stopifnot(names(genes_nSamples_atLeastMinReads) == rownames(countFilter_rnaseqDT))
    # convert # of samples to ratio of samples
    genes_ratioSamples_atLeastMinReads <- genes_nSamples_atLeastMinReads/totSamples
    stopifnot(genes_ratioSamples_atLeastMinReads >= 0)
    stopifnot(genes_ratioSamples_atLeastMinReads <= 1)
    stopifnot(length(genes_ratioSamples_atLeastMinReads) == nrow(countFilter_rnaseqDT))
    # which genes have a ratio of samples with at least min counts above the min. ratio of samples
    rowsToKeep <- genes_ratioSamples_atLeastMinReads >= min_sampleRatio
    stopifnot(length(rowsToKeep) == nrow(countFilter_rnaseqDT))


    stopifnot(names(rowsToKeep) == rownames(countFilter_rnaseqDT))
    keptRatio <- sum(rowsToKeep)/length(rowsToKeep)
    txt <- paste0(toupper(script_name), "> useFilterCountData is TRUE -> minCount-filtered geneList; to keep: ", sum(rowsToKeep), "/", length(rowsToKeep), "(", round(keptRatio*100,2), "%)\n")
    printAndLog(txt, pipLogFile)
} else {
stop("error")
}



# go further only with the ones that passed the filter
exprDT <- rnaseqDT[rowsToKeep, ]
geneList <- rownames(rnaseqDT)[rowsToKeep]
stopifnot(all(rownames(exprDT) == geneList))
stopifnot(nrow(exprDT) == sum(rowsToKeep))

cat("... prepare design (model matrix)\n")

# for the design:
my_group <- unlist(sapply(colnames(exprDT), function(x) {
  ifelse(x %in% samp1, cond1, ifelse(x %in% samp2, cond2, NA))
}))
stopifnot(!any(is.na(my_group)))
stopifnot(length(my_group) == ncol(exprDT))

# design matrix, cond1 as reference (for voom())
my_group_design <- factor(my_group, levels = c(cond1, cond2))
my_design <- model.matrix( ~ my_group_design)

outFile <- paste0(curr_outFold, "/", "boxplot_MDS_replicates.png")

if(inputDataType == "raw" | inputDataType == "RSEM"){
  cpm_expr_tmp <- cpm(exprDT, log=T)
  cpm_expr_tmp[is.na(cpm_expr_tmp)] <- 0
}else{
  cpm_expr_tmp <- exprDT
}

labcol <- unlist(sapply(colnames(cpm_expr_tmp), function(x) ifelse(x %in% samp1, "blue", "red") ))

# cat("... plot MDS\n") # do not do this, too slow when doing DE for all the genes !!! (and not the TAD-filtered genes)
# png(paste0(outFile), width=1000)
# par(mfrow=c(1,2), oma=par()$oma + c(5,0,0,0))
# boxplot(cpm_expr_tmp,las=2)
# plotMDS(cpm_expr_tmp, col=labcol)
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))

# update 12.10 => cannot create DGEList object for microarray EZH2 data (log-transformed; contain negative values)
# if the input data are already log-normalized or microarray data:
# -> cannot create DGEList object (needs non-negative values)
# -> do not apply voom, that is used before linear modeling to 
#    "Transform count data to log2-counts per million (logCPM), estimate the mean-variance relationship and use this to compute appropriate observation-level weights"

cat(paste0("... start DE analysis", "\t", Sys.time(), "\n"))
# printAndLog(txt, pipLogFile)

if(inputDataType == "raw" | inputDataType == "RSEM") {
  ## NORMAL RNASEQ DE ANALYSIS WITH RAW COUNTS
  # create the DGE object with the alredy filtered data
                      seqDataTmp <- DGEList(exprDT, group=my_group, genes = rownames(exprDT))
                      #calcNormFactors:
                      #it is usual to apply scale normalization to RNA-seq read counts, e.g. TMM normalization
                      seqDataTmp <- try(calcNormFactors(seqDataTmp))
                      if(class(seqDataTmp) == "try-error") {
                        seqData <- DGEList(exprDT, group=my_group, genes = rownames(exprDT))
                        seqData <- calcNormFactors(seqData, method = "none")
                        txt <- paste0(toupper(script_name), "> could not compute calcNormFactors with default method, used method = \"none\" \n")
                        printAndLog(txt, pipLogFile)
                        save(exprDT, file="exprDT.Rdata") #example: GSE52166_prePf_postPf
                      } else{
                        seqData <- DGEList(exprDT, group=my_group, genes = rownames(exprDT))
                        seqData <- calcNormFactors(seqData)
                      }
  outFile <- paste0(curr_outFold, "/", "mean_variance_voom.png")
  png(paste0(outFile))
  voomData <- voom(seqData, design = my_design, plot=T)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  vfitData <- lmFit(voomData, voomData$design)
  efitData <- eBayes(vfitData)
} else if(inputDataType == "microarray") {
  ## RNASEQ ANALYSIS FOR MICROARRAY
  ## https://support.bioconductor.org/p/67590/
  #fit <- lmFit(y, design)
  #fit <- eBayes(fit, trend=TRUE)
  ## gives (for microarray data) essentially the same effect as:
  #v <- vooma(y, design)
  #fit <- lmFit(v, design)
  #fit <- eBayes(fit)
  ##We generally recommend the former pipeline over the second for microarray data.
  vfitData <- lmFit(exprDT, my_design)
  efitData <- eBayes(vfitData, trend=TRUE)
} else if(inputDataType == "FPKM"){
  ## RNASEQ DE ANALYSIS FOR FPKM DATA
  #https://stat.ethz.ch/pipermail/bioconductor/2013-November/056309.html
  #If FPKM is really all you have, then convert the values to a log2 scale 
  #and do an ordinary limma analysis as you would for microarray data, using 
  #eBayes() with trend=TRUE. 
  exprDT <- exprDT + 0.0001
  exprDT <- log2(exprDT)
  vfitData <- lmFit(exprDT, my_design) 
  efitData <- eBayes(vfitData, trend = TRUE)
} else if(inputDataType == "DESeq2"){
  ## RNASEQ FOR ALREADY NORMALIZED/LOG-TRANSFORMED COUNTS
  vfitData <- lmFit(exprDT, my_design) 
  efitData <- eBayes(vfitData)
} else {
  stop("error\n")
}

cat(paste0("... end DE analysis", "\t", Sys.time(), "\n"))
# printAndLog(txt, pipLogFile)

if(inputDataType == "DESeq2" | inputDataType == "microarray" | inputDataType == "FPKM"){
  DE_topTable <- topTable(efitData, coef=ncol(my_design), number=Inf, sort.by="p")
  DE_topTable$genes <- rownames(DE_topTable)
}else if(inputDataType == "raw" | inputDataType == "RSEM") {
  DE_topTable <- topTable(efitData, coef=ncol(voomData$design), number=Inf, sort.by="p") 
} else{
  stop("error\n")
}
stopifnot(all(DE_topTable$genes %in% rownames(exprDT)))
stopifnot(all(DE_topTable$genes %in% names(rna_geneList)))

png(paste0(curr_outFold, "/", "MA_plot.png"), width=500)
plotMA(efitData)
foo <- dev.off()

png(paste0(curr_outFold, "/", "volcano_plot.png"), width=1000)
plot(DE_topTable$logFC, -log10(DE_topTable$adj.P.Val), pch=16, xlab="logFC", ylab="-log10(adj. p-val)")
plot(DE_topTable$logFC, -log10(DE_topTable$adj.P.Val), pch=16, xlab="logFC", ylab="-log10(adj. p-val)")
text(x = DE_topTable$logFC[1:10], y= -log10(DE_topTable$adj.P.Val[1:10]), 
     col="red", pch=16, labels = DE_topTable$genes[1:10])  
foo <- dev.off()

# take one gene to retrieve which condition was used as reference (in principle it is condition1/condition2 
# where condition1 is the 2nd in alphabetical order, e.g. lum/basal, mss/msi where negative
# logFC indicates more expressed in condition2, basal, msi, etc.)

x <- DE_topTable$genes[1]
exprCond1 <- mean(as.numeric(exprDT[x, samp1]))
exprCond2 <- mean(as.numeric(exprDT[x, samp2]))
txt <- paste0(toupper(script_name), "> Gene ", x, " logFC = ",round(DE_topTable$logFC[1],2), "; mean_", cond1, " = ", round(exprCond1,2), "; mean_", cond2, " = ", round(exprCond2,2), "\n")
printAndLog(txt, pipLogFile)

if(DE_topTable$logFC[1] > 0) {
  if(exprCond1 > exprCond2)
    txt <- paste0("> logFC > 0 when ", cond1, " > " , cond2, " => direction is ", cond1, "/", cond2, "\n")
  if(exprCond2 > exprCond1)
    txt <- paste0("> logFC > 0 when ", cond2, " > " , cond1, " => direction is ", cond2, "/", cond1, "\n")
} else if(DE_topTable$logFC[1] < 0) {
  if(exprCond1 > exprCond2)
    txt <- paste0("> logFC < 0 when ", cond1, " > " , cond2, " => direction is ", cond2, "/", cond1, "\n")
  if(exprCond2 > exprCond1)
    txt <- paste0("> logFC > 0 when ", cond2, " > " , cond1, " => direction is ", cond1, "/", cond2, "\n")
} 
printAndLog(txt, pipLogFile)
cat("... end DE\n... prepare output data\n")
#### PREPARE THE QQNORM DATA FOR THE GENES I USED FOR DE ANALYSIS
if(inputDataType == "microarray") {
  cat("... madnorm the data for other analyses \n")
  madnorm_exprDT <- madNorm(exprDT)
  stopifnot(all(dim(madnorm_exprDT) == dim(exprDT)))
  rownames(madnorm_exprDT) <- rownames(exprDT)
  colnames(madnorm_exprDT) <- colnames(exprDT)
} else {
  cat("... qqnorm the data for other analyses \n")
  qqnorm_exprDT <- t(apply(exprDT, 1, quantNorm))
  stopifnot(all(dim(qqnorm_exprDT) == dim(exprDT)))
  rownames(qqnorm_exprDT) <- rownames(exprDT)
  colnames(qqnorm_exprDT) <- colnames(exprDT)
}
#### PREPARE DATA TO WRITE IN FILES
if(inputDataType == "raw" | inputDataType == "RSEM") stopifnot(all(dim(voomData$E) == dim(exprDT)))
stopifnot(all(rownames(exprDT) == geneList))

DE_rnaseqDT <- exprDT
if(inputDataType == "microarray") {
  DE_madnorm_rnaseqDT <- madnorm_exprDT
} else {
  DE_qqnorm_rnaseqDT <- qqnorm_exprDT 
}
DE_geneList <- rna_geneList[rownames(exprDT)]
stopifnot(length(DE_geneList) > 0)

##### WRITE DATA IN FILE
cat("... write data in files\n")
# the expression data used for the DE analysis (i.e. CAGE seq after filtering minimum cpm count)
save(DE_rnaseqDT, file = paste0(curr_outFold, "/", "DE_rnaseqDT.Rdata"))
cat(paste0("... written: ", paste0(curr_outFold, "/", "DE_rnaseqDT.Rdata"), "\n"))
# the same but qqnorm
if(inputDataType == "microarray") {
  save(DE_madnorm_rnaseqDT, file = paste0(curr_outFold, "/", "DE_madnorm_rnaseqDT.Rdata"))
  cat(paste0("... written: ", paste0(curr_outFold, "/", "DE_madnorm_rnaseqDT.Rdata"), "\n"))
} else{
  save(DE_qqnorm_rnaseqDT, file = paste0(curr_outFold, "/", "DE_qqnorm_rnaseqDT.Rdata"))
  cat(paste0("... written: ", paste0(curr_outFold, "/", "DE_qqnorm_rnaseqDT.Rdata"), "\n")) 
}
# the DE topTable
save(DE_topTable, file = paste0(curr_outFold, "/", "DE_topTable.Rdata"))
cat(paste0("... written: ", paste0(curr_outFold, "/", "DE_topTable.Rdata"), "\n"))
# gene list
save(DE_geneList, file = paste0(curr_outFold, "/", "DE_geneList.Rdata"))
cat(paste0("... written: ", paste0(curr_outFold, "/", "DE_geneList.Rdata"), "\n"))


##################################################################################### ADDED 20.06.2018

# STOPIFNOT IN SCRIPT 8: all(pipeline_geneList %in% DE_topTable$genes) 

check_pipeline_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name,  "/", "pipeline_geneList.Rdata"))))
check_DE_topTable <- eval(parse(text = load(paste0(curr_outFold, "/", "DE_topTable.Rdata"))))

#stopifnot(check_pipeline_geneList %in% check_DE_topTable$genes)
#update correct 17.08.2018
stopifnot(names(check_pipeline_geneList) %in% check_DE_topTable$genes)

#####################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))
