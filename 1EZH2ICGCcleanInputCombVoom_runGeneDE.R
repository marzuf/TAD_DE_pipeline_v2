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


# UPDATE 20.06.2018
# add specification for input data type
stopifnot(normRNAseq)
cat(paste0("> normRNAseq: ", as.character(normRNAseq), "\n"))

stopifnot(!microarray)

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
# UPDATE 20.06.2018 for ICGC
if(normRNAseq) {
  cpm_exprDT <- cpm(rnaseqDT)
  txt <- paste0(toupper(script_name), "> NA in cpm_exprDT: ", 
                sum(is.na(cpm_exprDT)), "/", dim(cpm_exprDT)[1]*dim(cpm_exprDT)[2], " (",
                round((sum(is.na(cpm_exprDT))/(dim(cpm_exprDT)[1]*dim(cpm_exprDT)[2]) * 100),2), "%)\n")
  printAndLog(txt, pipLogFile)
  rowsToKeep <- rowSums(cpm_exprDT, na.rm=T) >= (minCpmRatio * ncol(rnaseqDT))
  
  txt <- paste0(toupper(script_name), "> CPM filter, genes to retain: ", sum(rowsToKeep), "/", nrow(rnaseqDT), "\n")
  printAndLog(txt, pipLogFile)
} else {
  stop("ERROR\n")
}


#############################################################################  UPDATE 29.01.2018 >>> SPLIT EXPRDT TO CREATE 2 VOOMDATA: VOOM DATA ALL

stopifnot(exists("cohort1samples"))
stopifnot(exists("cohort2samples"))

stopifnot(cohort1samples %in% colnames(rnaseqDT))
stopifnot(cohort2samples %in% colnames(rnaseqDT))

# go further only with the ones that passed the filter
exprDT <- rnaseqDT[rowsToKeep, c(cohort1samples, cohort2samples)]
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



#############################################################################  UPDATE 29.01.2018 >>> SPLIT EXPRDT TO CREATE 2 VOOMDATA: VOOM DATA 1

# go further only with the ones that passed the filter
exprDT1 <- rnaseqDT[rowsToKeep, c(cohort1samples)]

stopifnot(all(rownames(exprDT1) == geneList))
stopifnot(nrow(exprDT1) == sum(rowsToKeep))

cat("... prepare design (model matrix)\n")

# for the design:
my_group1 <- unlist(sapply(colnames(exprDT1), function(x) {
  ifelse(x %in% samp1, cond1, ifelse(x %in% samp2, cond2, NA))
}))
stopifnot(!any(is.na(my_group1)))
stopifnot(length(my_group1) == ncol(exprDT1))

# design matrix, cond1 as reference (for voom())
my_group_design1 <- factor(my_group1, levels = c(cond1, cond2))
my_design1 <- model.matrix( ~ my_group_design1)



#############################################################################  UPDATE 29.01.2018 >>> SPLIT EXPRDT TO CREATE 2 VOOMDATA: VOOM DATA 1

# go further only with the ones that passed the filter
exprDT2 <- rnaseqDT[rowsToKeep, c(cohort2samples)]

stopifnot(all(rownames(exprDT2) == geneList))
stopifnot(nrow(exprDT2) == sum(rowsToKeep))

cat("... prepare design (model matrix)\n")

# for the design:
my_group2 <- unlist(sapply(colnames(exprDT2), function(x) {
  ifelse(x %in% samp1, cond1, ifelse(x %in% samp2, cond2, NA))
}))
stopifnot(!any(is.na(my_group2)))
stopifnot(length(my_group2) == ncol(exprDT2))

# design matrix, cond2 as reference (for voom())
my_group_design2 <- factor(my_group2, levels = c(cond1, cond2))
my_design2 <- model.matrix( ~ my_group_design2)

####################################
outFile <- paste0(curr_outFold, "/", "boxplot_MDS_replicates.png")

if(normRNAseq){
  cpm_expr_tmp <- cpm(exprDT, log=T)
  cpm_expr_tmp[is.na(cpm_expr_tmp)] <- 0
}else{
  stop("error\n")
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

if(normRNAseq) {
  ### UPDATE FOR ICGC !! 20.06.2018
  ## NORMAL RNASEQ DE ANALYSIS WITH RAW COUNTS
  # create the DGE object with the alredy filtered data
  
  #############################################################################  UPDATE 29.01.2018 >>> SPLIT EXPRDT TO CREATE 2 VOOMDATA: VOOM DATA ALL
  
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
                      outFile <- paste0(curr_outFold, "/", "mean_variance_voomAllData.png")
                      png(paste0(outFile))
                      voomData <- voom(seqData, design = my_design, plot=T)
                      foo <- dev.off()
  #############################################################################  UPDATE 29.01.2018 >>> SPLIT EXPRDT TO CREATE 2 VOOMDATA: VOOM DATA 1
                      seqDataTmp1 <- DGEList(exprDT1, group=my_group1, genes = rownames(exprDT1))
                      #calcNormFactors:
                      #it is usual to apply scale normalization to RNA-seq read counts, e.g. TMM normalization
                      seqDataTmp1 <- try(calcNormFactors(seqDataTmp1))
                      if(class(seqDataTmp1) == "try-error") {
                        seqData1 <- DGEList(exprDT1, group=my_group1, genes = rownames(exprDT1))
                        seqData1 <- calcNormFactors(seqData1, method = "none")
                        txt <- paste0(toupper(script_name), "> could not compute calcNormFactors with default method, used method = \"none\" \n")
                        printAndLog(txt, pipLogFile)
                        save(exprDT1, file="exprDT1.Rdata") #example: GSE52166_prePf_postPf
                      } else{
                        seqData1 <- DGEList(exprDT1, group=my_group1, genes = rownames(exprDT1))
                        seqData1 <- calcNormFactors(seqData1)
                      }
                      outFile <- paste0(curr_outFold, "/", "mean_variance_voomData1.png")
                      png(paste0(outFile))
                      voomData1 <- voom(seqData1, design = my_design, plot=T)
                      foo <- dev.off()
                      
                      
  #############################################################################  UPDATE 29.01.2018 >>> SPLIT EXPRDT TO CREATE 2 VOOMDATA: VOOM DATA 2                   
                      seqDataTmp2 <- DGEList(exprDT2, group=my_group2, genes = rownames(exprDT2))
                      #calcNormFactors:
                      #it is usual to apply scale normalization to RNA-seq read counts, e.g. TMM normalization
                      seqDataTmp2 <- try(calcNormFactors(seqDataTmp2))
                      if(class(seqDataTmp2) == "try-error") {
                        seqData2 <- DGEList(exprDT2, group=my_group2, genes = rownames(exprDT2))
                        seqData2 <- calcNormFactors(seqData2, method = "none")
                        txt <- paste0(toupper(script_name), "> could not compute calcNormFactors with default method, used method = \"none\" \n")
                        printAndLog(txt, pipLogFile)
                        save(exprDT2, file="exprDT2.Rdata") #example: GSE52166_prePf_postPf
                      } else{
                        seqData2 <- DGEList(exprDT2, group=my_group2, genes = rownames(exprDT2))
                        seqData2 <- calcNormFactors(seqData2)
                      }
                      outFile <- paste0(curr_outFold, "/", "mean_variance_voomData2.png")
                      png(paste0(outFile))
                      voomData2 <- voom(seqData2, design = my_design, plot=T)
                      foo <- dev.off()
                      

  ############################################################################# combine both voom data
                      # length(voomData1@.Data)
                      # # [[1]] -> list of genes
                      # all(voomData1@.Data[[1]] == voomData2@.Data[[1]]) # TRUE
                      # # -> take those from voomData
                      
                      # [[2]] => lib size and norm fact.
                      # voomData1@.Data[[2]]
                      voomDataComb_data2 <- rbind(voomData1@.Data[[2]],voomData2@.Data[[2]])
                      stopifnot(dim(voomDataComb_data2) == dim(voomData@.Data[[2]]))
                      
                      # [[3]] => gene x sample matrix
                      # voomData1@.Data[[3]]
                      voomDataComb_data3 <- cbind(voomData1@.Data[[3]],voomData2@.Data[[3]])
                      stopifnot(dim(voomDataComb_data3) == dim(voomData@.Data[[3]]))
                      
                      # [[4]] => gene x sample matrix
                      voomData1@.Data[[4]]
                      voomDataComb_data4 <- cbind(voomData1@.Data[[4]],voomData2@.Data[[4]])
                      stopifnot(dim(voomDataComb_data4) == dim(voomData@.Data[[4]]))
                      
                      # # [[5]] => contrasts
                      # voomData1@.Data[[5]]
                      # voomData2@.Data[[5]]
                      # # -> take those from voomData
                      
                      voomDataComb <- voomData
                      rm(voomData)
                      voomDataComb[[2]] <- voomDataComb_data2
                      voomDataComb[[3]] <- voomDataComb_data3
                      voomDataComb[[4]] <- voomDataComb_data4

  voomData <- voomDataComb  
                      
  cat(paste0("... written: ", outFile, "\n"))
  vfitData <- lmFit(voomData, voomData$design)
  efitData <- eBayes(vfitData)
} else {
  stop("error\n")
}

cat(paste0("... end DE analysis", "\t", Sys.time(), "\n"))
# printAndLog(txt, pipLogFile)

if(normRNAseq){
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
# UPDATE 20.06.2018 FOR ICGC
if(normRNAseq) {
  cat("... qqnorm the data for other analyses \n")
  qqnorm_exprDT <- t(apply(exprDT, 1, quantNorm))
  stopifnot(all(dim(qqnorm_exprDT) == dim(exprDT)))
  rownames(qqnorm_exprDT) <- rownames(exprDT)
  colnames(qqnorm_exprDT) <- colnames(exprDT)
} else {
  stop("ERROR\n")
}

#### PREPARE DATA TO WRITE IN FILES
## UPDATE 20.06.2018 ICGC DATA
if(normRNAseq) {
  stopifnot(all(dim(voomData$E) == dim(exprDT)))
} else {
  stop("error\n")
}
stopifnot(all(rownames(exprDT) == geneList))

DE_rnaseqDT <- exprDT

if(normRNAseq) {
  DE_qqnorm_rnaseqDT <- qqnorm_exprDT 
} else {
  stop("error\n")
}

DE_geneList <- rna_geneList[rownames(exprDT)]
stopifnot(length(DE_geneList) > 0)

##### WRITE DATA IN FILE
cat("... write data in files\n")
# the expression data used for the DE analysis (i.e. CAGE seq after filtering minimum cpm count)
save(DE_rnaseqDT, file = paste0(curr_outFold, "/", "DE_rnaseqDT.Rdata"))
cat(paste0("... written: ", paste0(curr_outFold, "/", "DE_rnaseqDT.Rdata"), "\n"))
# the same but qqnorm
if(normRNAseq) {
  save(DE_qqnorm_rnaseqDT, file = paste0(curr_outFold, "/", "DE_qqnorm_rnaseqDT.Rdata"))
  cat(paste0("... written: ", paste0(curr_outFold, "/", "DE_qqnorm_rnaseqDT.Rdata"), "\n")) 
} else{
  stop("error\n")
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
