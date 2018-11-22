#################################################################### WITH ALL DATA
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
seqDataTmp <- DGEList(exprDT, group=my_group, genes = rownames(exprDT))
seqData <- DGEList(exprDT, group=my_group, genes = rownames(exprDT))
seqData <- calcNormFactors(seqData)
voomData <- voom(seqData, design = my_design, plot=T)
vfitData <- lmFit(voomData, voomData$design)
efitData <- eBayes(vfitData)

#################################################################### FEW SAMPLES 1

nSamp1 <- 1:20
nSamp2 <- 21:54

# go further only with the ones that passed the filter
exprDT1 <- rnaseqDT[rowsToKeep, nSsamp1]
geneList <- rownames(rnaseqDT)[rowsToKeep]
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
my_group1_design <- factor(my_group1, levels = c(cond1, cond2))
my_design1 <- model.matrix( ~ my_group1_design)
seqData1 <- DGEList(exprDT1, group=my_group1, genes = rownames(exprDT1))
seqData1 <- calcNormFactors(seqData1)
voomData1 <- voom(seqData1, design = my_design1, plot=T)
vfitData1 <- lmFit(voomData1, voomData1$design)
efitData1 <- eBayes(vfitData1)


#################################################################### FEW SAMPLES 2
# go further only with the ones that passed the filter
exprDT2 <- rnaseqDT[rowsToKeep, nSamp2]
geneList <- rownames(rnaseqDT)[rowsToKeep]
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
my_group2_design <- factor(my_group2, levels = c(cond1, cond2))
my_design2 <- model.matrix( ~ my_group2_design)
seqData2 <- DGEList(exprDT2, group=my_group2, genes = rownames(exprDT2))
seqData2 <- calcNormFactors(seqData2)
voomData2 <- voom(seqData2, design = my_design2, plot=T)
vfitData2 <- lmFit(voomData2, voomData2$design)
efitData2 <- eBayes(vfitData2)

#################################################################### FEW SAMPLES 2
length(voomData1@.Data)
# [[1]] -> list of genes
all(voomData1@.Data[[1]] == voomData2@.Data[[1]]) # TRUE
# -> take those from voomData

# [[2]] => lib size and norm fact.
voomData1@.Data[[2]]
voomDataComb_data2 <- rbind(voomData1@.Data[[2]],voomData2@.Data[[2]])
dim(voomDataComb_data2) == dim(voomData@.Data[[2]])

# [[3]] => gene x sample matrix
voomData1@.Data[[3]]
voomDataComb_data3 <- cbind(voomData1@.Data[[3]],voomData2@.Data[[3]])
dim(voomDataComb_data3) == dim(voomData@.Data[[3]])

# [[4]] => gene x sample matrix
voomData1@.Data[[4]]
voomDataComb_data4 <- cbind(voomData1@.Data[[4]],voomData2@.Data[[4]])
dim(voomDataComb_data4) == dim(voomData@.Data[[4]])

# [[5]] => contrasts
voomData1@.Data[[5]]
voomData2@.Data[[5]]
# -> take those from voomData

voomDataComb <- voomData
voomDataComb[[2]] <- voomDataComb_data2
voomDataComb[[3]] <- voomDataComb_data3
voomDataComb[[4]] <- voomDataComb_data4
vfitDataComb <- lmFit(voomDataComb, voomDataComb$design)
efitDataComb <- eBayes(vfitDataComb)


