#!/usr/bin/Rscript


cat(paste0("> START ", "saveVoomOnly",  "\n"))


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

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_saveVoom_logFile.txt")
system(paste0("rm -f ", pipLogFile))

stopifnot(file.exists(paste0(pipOutFold, "/", script0_name, "/rna_rnaseqDT.Rdata")))
rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_rnaseqDT.Rdata"))))
if(ncol(rnaseqDT) >= 5 & nrow(rnaseqDT) >= 5)
    rnaseqDT[1:5,1:5]
initRowNbr <- nrow(rnaseqDT)
stopifnot(is.numeric(rnaseqDT[1,1]))

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

# FILTER THE EXPRESSION DATA TO MIN CPM FILTER  ################################################################ CPM FOR MICROARRAY OR NOT ????????????????????????????
# apply CPM only on raw counts or RSEM
if(inputDataType == "raw" | inputDataType == "RSEM") {
  cpm_exprDT <- cpm(rnaseqDT)
  txt <- paste0(toupper(script_name), "> NA in cpm_exprDT: ", 
                sum(is.na(cpm_exprDT)), "/", dim(cpm_exprDT)[1]*dim(cpm_exprDT)[2], " (",
                round((sum(is.na(cpm_exprDT))/(dim(cpm_exprDT)[1]*dim(cpm_exprDT)[2]) * 100),2), "%)\n")
  printAndLog(txt, pipLogFile)
  rowsToKeep <- rowSums(cpm_exprDT, na.rm=T) >= (minCpmRatio * ncol(rnaseqDT))
  
  txt <- paste0(toupper(script_name), "> CPM filter, genes to retain: ", sum(rowsToKeep), "/", nrow(rnaseqDT), "\n")
  printAndLog(txt, pipLogFile)
} else if(inputDataType == "FPKM") {
  txt <- paste0(toupper(script_name), "> !!! FPKM filter not applied !!!", "\n")
  printAndLog(txt, pipLogFile)
  cpm_exprDT <- rnaseqDT
  txt <- paste0(toupper(script_name), "> NA in cpm_exprDT: ", 
                sum(is.na(cpm_exprDT)), "/", dim(cpm_exprDT)[1]*dim(cpm_exprDT)[2], " (",
                round((sum(is.na(cpm_exprDT))/(dim(cpm_exprDT)[1]*dim(cpm_exprDT)[2]) * 100),2), "%)\n")
  printAndLog(txt, pipLogFile)
  rowsToKeep <- rowSums(cpm_exprDT, na.rm=T) >= (minCpmRatio * ncol(rnaseqDT))
} else if (inputDataType == "microarray" | inputDataType == "DESeq2") {
  txt <- paste0(toupper(script_name), "> !!! CPM filter not applied !!!", "\n")
  printAndLog(txt, pipLogFile)
  rowsToKeep <- rep(TRUE, nrow(rnaseqDT))
} else {
  stop("ERROR\n")
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

# outFile <- paste0(curr_outFold, "/", "boxplot_MDS_replicates.png")

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
                        # save(exprDT, file="exprDT.Rdata") #example: GSE52166_prePf_postPf
                      } else{
                        seqData <- DGEList(exprDT, group=my_group, genes = rownames(exprDT))
                        seqData <- calcNormFactors(seqData)
                      }
  # outFile <- paste0(curr_outFold, "/", "mean_variance_voom.png")
  # png(paste0(outFile))
  # voomData <- voom(seqData, design = my_design, plot=T)
  # foo <- dev.off()
  # cat(paste0("... written: ", outFile, "\n"))
  voomData <- voom(seqData, design = my_design, plot=F)
  # vfitData <- lmFit(voomData, voomData$design)
  # efitData <- eBayes(vfitData)
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
  exprDT <- exprDT
  # vfitData <- lmFit(exprDT, my_design)
  # efitData <- eBayes(vfitData, trend=TRUE)
} else if(inputDataType == "FPKM"){
  ## RNASEQ DE ANALYSIS FOR FPKM DATA
  #https://stat.ethz.ch/pipermail/bioconductor/2013-November/056309.html
  #If FPKM is really all you have, then convert the values to a log2 scale 
  #and do an ordinary limma analysis as you would for microarray data, using 
  #eBayes() with trend=TRUE. 
  exprDT <- exprDT + 0.0001
  exprDT <- log2(exprDT)
  # vfitData <- lmFit(exprDT, my_design) 
  # efitData <- eBayes(vfitData, trend = TRUE)
} else if(inputDataType == "DESeq2"){
  ## RNASEQ FOR ALREADY NORMALIZED/LOG-TRANSFORMED COUNTS
  exprDT <- exprDT
  # vfitData <- lmFit(exprDT, my_design) 
  # efitData <- eBayes(vfitData)
} else {
  stop("error\n")
}

voom_lmFitInputDT <- exprDT
save(voom_lmFitInputDT, file = paste0(curr_outFold, "/", "voom_lmFitInputDT.Rdata"))
cat(paste0("... written: ", paste0(curr_outFold, "/", "voom_lmFitInputDT.Rdata"), "\n"))

#####################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))
