startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: rna_rnaseqDT.Rdata
# - script0: rna_geneList.Rdata
# - script0: pipeline_geneList.Rdata
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
cat("> ! WARNING: for paired samples\n")

source("main_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(edgeR, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  
suppressPackageStartupMessages(library(limma, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  

if(!exists("designFile") | is.null(designFile)) {
  stop("! designFile not provided, should maybe run 1_runGeneDE.R (not paired)\n")
}
stopifnot(file.exists(designFile))
designDT <- read.delim(designFile, header=TRUE, stringsAsFactors = FALSE, col.names=c("samples", "pair", "condition"))

# if microarray was not set in the settings file -> by default set to  FALSE
if(!exists("microarray")) microarray <- FALSE

# if FPKM was not set in the settings file -> by default set to  FALSE
if(!exists("FPKM")) FPKM <- FALSE


# for some datasets, the count matrix contains negative values => run without voom and no cpm
if(any(unlist(sapply(runWithoutVoom, function(x) regexpr(x, settingF) > 0)))) {
  applyVoomAndCPM <- FALSE
} else{
  applyVoomAndCPM <- TRUE
}

if(microarray) applyVoomAndCPM <- FALSE

cat(paste0("> applyVoomAndCPM: ", as.character(applyVoomAndCPM), "\n"))
cat(paste0("> microarray data: ", as.character(microarray), "\n"))
cat(paste0("> FPKM data: ", as.character(FPKM), "\n"))

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

# ADDED 27.11.2018 to check using other files
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)


stopifnot(file.exists(paste0(pipOutFold, "/", script0_name, "/rna_rnaseqDT.Rdata")))
rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_rnaseqDT.Rdata"))))
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
txt <- paste0(toupper(script_name), "> Number of genes available for DE analysis: ", nrow(rnaseqDT), "/", initRowNbr, "\n")
printAndLog(txt, pipLogFile)

# RUN THE DE ANALYSIS
samp1 <- eval(parse(text=load(paste0(setDir, "/", sample1_file))))
samp2 <- eval(parse(text=load(paste0(setDir, "/", sample2_file))))
# ensure the samples are present in the column names
stopifnot(all(samp1 %in% colnames(rnaseqDT)))
stopifnot(all(samp2 %in% colnames(rnaseqDT)))

## =>>>>> PAIRED GENE DE ANALYSIS added here check paired info is available for the samples !!!
stopifnot(all(samp1 %in% designDT$samples))
stopifnot(all(samp2 %in% designDT$samples))

if(! all(designDT$samples %in% c(samp1,samp2))) {
  cat("!!! WARNING: some samples from designFile are not available in expression table !\n")
}

designDT <- designDT[designDT$samples %in% c(samp1,samp2),]
stopifnot(nrow(designDT) > 0)

stopifnot(designDT$condition %in% c(cond1,cond2))
stopifnot(c(cond1, cond2) %in% designDT$condition)

designDT$condition <- factor(as.character(designDT$condition), levels=c(cond1, cond2))
designDT$pair <- factor(as.character(designDT$pair), levels = as.character(unique(designDT$pair)))
designDT <- designDT[order(as.numeric(designDT$condition), designDT$pair),]

# reorder samp1 and samp2 so that they match the design (not sure it is necessary, but seems cleaner)
samp1b <- samp1[match(designDT$samples[designDT$samples %in% samp1], samp1)]
stopifnot(setequal(samp1b, samp1))
samp2b <- samp2[match(designDT$samples[designDT$samples %in% samp2], samp2)]
stopifnot(setequal(samp2b, samp2))
samp1 <- samp1b
samp2 <- samp2b
## <<<<<=
rnaseqDT <- rnaseqDT[,c(samp1, samp2)]

# FILTER THE EXPRESSION DATA TO MIN CPM FILTER  ################################################################ CPM FOR MICROARRAY OR NOT ????????????????????????????
if(applyVoomAndCPM & !microarray) {
  cpm_exprDT <- cpm(rnaseqDT)
  txt <- paste0(toupper(script_name), "> NA in cpm_exprDT: ", 
                sum(is.na(cpm_exprDT)), "/", dim(cpm_exprDT)[1]*dim(cpm_exprDT)[2], " (",
                round((sum(is.na(cpm_exprDT))/(dim(cpm_exprDT)[1]*dim(cpm_exprDT)[2]) * 100),2), "%)\n")
  printAndLog(txt, pipLogFile)
  rowsToKeep <- rowSums(cpm_exprDT, na.rm=T) >= (minCpmRatio * ncol(rnaseqDT))
  
  txt <- paste0(toupper(script_name), "> CPM filter, genes to retain: ", sum(rowsToKeep), "/", nrow(rnaseqDT), "\n")
  printAndLog(txt, pipLogFile)
} else{
  txt <- paste0(toupper(script_name), "> !!! CPM filter not applied !!!", "\n")
  printAndLog(txt, pipLogFile)
  rowsToKeep <- rep(TRUE, nrow(rnaseqDT))
}


# go further only with the ones that passed the filter
exprDT <- rnaseqDT[rowsToKeep, ]
geneList <- rownames(rnaseqDT)[rowsToKeep]
stopifnot(all(rownames(exprDT) == geneList))

## =>>>>> PAIRED GENE DE ANALYSIS 
########################### !!!! HERE DESIGN CHANGED FOR PAIRED SAMPLES
# Limma Paired
# > targets[1:6,]
#          cell pair treatment
# 1  GSM1196699    1      DMSO
# 2  GSM1196700    2      DMSO
# 3  GSM1196701    3      DMSO
# 4  GSM1196705    1       EPZ
# 5  GSM1196706    2       EPZ
# 6  GSM1196707    3       EPZ
# 
# Pair<-factor(targets$pair)
# Treat<-factor(targets$treatment,levels=c("DMSO","EPZ"))
# design<-model.matrix(~Pair+Treat)
# 
# 
# fit_pair<-lmFit(m,design)
# fit_pair <- eBayes(fit_pair)
# results = topTable(fit_pair, coef="TreatEPZ", number = nrow(m))

Pair <- designDT$pair
Treat <- designDT$condition
my_design <- model.matrix(~Pair + Treat)

                # # for the design:
                # my_group <- unlist(sapply(colnames(exprDT), function(x) {
                #   ifelse(x %in% samp1, cond1, ifelse(x %in% samp2, cond2, NA))
                # }))
                # stopifnot(!any(is.na(my_group)))
                # stopifnot(length(my_group) == ncol(exprDT))

                # # design matrix, cond1 as reference (for voom())
                # my_group_design <- factor(my_group, levels = c(cond1, cond2))
                # my_design <- model.matrix( ~ my_group_design)

outFile <- paste0(curr_outFold, "/", "boxplot_MDS_replicates.png")

if(applyVoomAndCPM & !microarray){
  cpm_expr_tmp <- cpm(exprDT, log=T)
  cpm_expr_tmp[is.na(cpm_expr_tmp)] <- 0
}else{
  cpm_expr_tmp <- exprDT
}

labcol <- unlist(sapply(colnames(cpm_expr_tmp), function(x) ifelse(x %in% samp1, "blue", "red") ))

png(paste0(outFile), width=1000)
par(mfrow=c(1,2), oma=par()$oma + c(5,0,0,0))
boxplot(cpm_expr_tmp,las=2)
plotMDS(cpm_expr_tmp, col=labcol)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

# update 12.10 => cannot create DGEList object for microarray EZH2 data (log-transformed; contain negative values)
# if the input data are already log-normalized or microarray data:
# -> cannot create DGEList object (needs non-negative values)
# -> do not apply voom, that is used before linear modeling to 
#    "Transform count data to log2-counts per million (logCPM), estimate the mean-variance relationship and use this to compute appropriate observation-level weights"

if(!microarray & applyVoomAndCPM) {
  ## NORMAL RNASEQ DE ANALYSIS WITH RAW COUNTS
  # create the DGE object with the alredy filtered data
  seqDataTmp <- DGEList(exprDT, group=my_group, genes = rownames(exprDT))
  #calcNormFactors: 
  #it is usual to apply scale normalization to RNA-seq read counts, e.g. TMM normalization 
  seqDataTmp <- try(calcNormFactors(seqDataTmp))
  if(class(seqDataTmp) == "try-error") {
    seqData <- DGEList(exprDT, group=my_group, genes = rownames(exprDT))
    seqData <- calcNormFactors(seqData, method = "none")
    txt <- paste0(toupper(script_name), "> could not compute calcNormFactors with default method, used default = \"none\" \n")
    printAndLog(txt, pipLogFile)
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
} else if(microarray) {
  ## RNASEQ ANALYSIS FOR MICROARRAY
  ## https://support.bioconductor.org/p/67590/
  #fit <- lmFit(y, design)
  #fit <- eBayes(fit, trend=TRUE)
  ## gives (for microarray data) essentially the same effect as:
  #v <- vooma(y, design)
  #fit <- lmFit(v, design)
  #fit <- eBayes(fit)
  ##We generally recommend the former pipeline over the second for microarray data.
  outFile <- paste0(curr_outFold, "/", "mean_variance_vooma.png")
  png(paste0(outFile))
  #voomData <- vooma(exprDT, design = my_design, plot=T) CORRECTED 19.12 - WAS WRONG !!!
  foo <- dev.off()
  vfitData <- lmFit(exprDT, my_design)
  efitData <- eBayes(vfitData, trend=TRUE)
} else if(FPKM){
  ## RNASEQ DE ANALYSIS FOR FPKM DATA
  #https://stat.ethz.ch/pipermail/bioconductor/2013-November/056309.html
  #If FPKM is really all you have, then convert the values to a log2 scale 
  #and do an ordinary limma analysis as you would for microarray data, using 
  #eBayes() with trend=TRUE. 
  exprDT <- exprDT + 0.0001
  exprDT <- log2(exprDT)
  vfitData <- lmFit(exprDT, my_design) 
  efitData <- eBayes(vfitData, trend = TRUE)
} else {
  ## RNASEQ FOR ALREADY NORMALIZED/LOG-TRANSFORMED COUNTS
  vfitData <- lmFit(exprDT, my_design) 
  efitData <- eBayes(vfitData)
}

if(microarray | !applyVoomAndCPM | FPKM){
  DE_topTable <- topTable(efitData, coef=ncol(my_design), number=Inf, sort.by="p")
  DE_topTable$genes <- rownames(DE_topTable)
}else{
  DE_topTable <- topTable(efitData, coef=ncol(voomData$design), number=Inf, sort.by="p") 
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
cat("... end DE\n...prepare output data\n")
#### PREPARE THE QQNORM DATA FOR THE GENES I USED FOR DE ANALYSIS
if(microarray) {
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
if(!microarray & applyVoomAndCPM) stopifnot(all(dim(voomData$E) == dim(exprDT)))
stopifnot(all(rownames(exprDT) == geneList))

DE_rnaseqDT <- exprDT
if(microarray) {
  DE_madnorm_rnaseqDT <- madnorm_exprDT
} else {
  DE_qqnorm_rnaseqDT <- qqnorm_exprDT 
}
DE_geneList <- rna_geneList[rownames(exprDT)]

##### WRITE DATA IN FILE
cat("... write data in files\n")
# the expression data used for the DE analysis (i.e. CAGE seq after filtering minimum cpm count)
save(DE_rnaseqDT, file = paste0(curr_outFold, "/", "DE_rnaseqDT.Rdata"))
cat(paste0("... written: ", paste0(curr_outFold, "/", "DE_rnaseqDT.Rdata"), "\n"))
# the same but qqnorm
if(microarray) {
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

stopifnot(check_pipeline_geneList %in% check_DE_topTable$genes)

#####################################################################################


txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))
