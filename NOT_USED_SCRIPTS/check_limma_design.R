startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: rna_rnaseqDT.Rdata
# - script0: rna_geneList.Rdata
# - script0: pipeline_geneList.Rdata
################################################################################

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

settingF <- "SETTING_FILES/run_settings_TCGAbrca_lum_bas.R"
# args <- commandArgs(trailingOnly = TRUE)
# stopifnot(length(args) == 1)
# settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline")

script0_name <- "0_prepGeneData"
script_name <- "1_runGeneDE"
stopifnot(file.exists(paste0(pipScriptDir, "/",script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(edgeR, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  
suppressPackageStartupMessages(library(limma, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  

printAndLog <- function(txt, logFile=NULL){
  if(!is.null(logFile)){
    cat(txt, file=logFile, append=T)
  }
  cat(txt)
}

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
curr_outFold <- NULL
# system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- NULL

stopifnot(file.exists(paste0(pipOutFold, "/", script0_name, "/rna_rnaseqDT.Rdata")))
rnaseqDT <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_rnaseqDT.Rdata"))))
rnaseqDT[1:5,1:5]
initRowNbr <- nrow(rnaseqDT)
stopifnot(is.numeric(rnaseqDT[1,1]))

# TAKE ONLY THE GENES FOR WHICH I HAVE POSITIONS
stopifnot(file.exists(paste0(pipOutFold, "/", script0_name,  "/", "rna_geneList.Rdata")))
init_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name,  "/", "rna_geneList.Rdata"))))
# => UPDATE: TAKE ONLY THE GENE LIST PREPARED IN 0_prepGeneData ACCORDING TO CURRENT SETTINGS
stopifnot(file.exists(paste0(pipOutFold, "/", script0_name,  "/", "pipeline_geneList.Rdata")))
rna_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name,  "/", "pipeline_geneList.Rdata"))))

txt <- paste0(toupper(script_name), "> Start with # genes: ", length(rna_geneList), "/", length(init_geneList), "\n")
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

##############################################################################################################################################################
################################################################################################## VERSION DESIGN 1
##############################################################################################################################################################

# for the design:
my_group <- unlist(sapply(colnames(exprDT), function(x) {
  ifelse(x %in% samp1, cond1, ifelse(x %in% samp2, cond2, NA))
}))
stopifnot(!any(is.na(my_group)))
stopifnot(length(my_group) == ncol(exprDT))

# design matrix, cond1 as reference (for voom())
my_group_design <- factor(my_group, levels = c(cond1, cond2))
my_design <- model.matrix( ~ my_group_design)

if(applyVoomAndCPM & !microarray){
  cpm_expr_tmp <- cpm(exprDT, log=T)
  cpm_expr_tmp[is.na(cpm_expr_tmp)] <- 0
}else{
  cpm_expr_tmp <- exprDT
}

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
voomData <- voom(seqData, design = my_design, plot=FALSE)
vfitData <- lmFit(voomData, voomData$design)
efitData <- eBayes(vfitData)

# test2
v2fitData <- lmFit(voomData, voomData$design)
v2fitData <- eBayes(v2fitData)

stopifnot(all.equal(efitData, v2fitData))

if(microarray | !applyVoomAndCPM | FPKM){
  DE_topTable <- topTable(efitData, coef=ncol(my_design), number=Inf, sort.by="p")
  DE_topTable$genes <- rownames(DE_topTable)
}else{
  DE_topTable <- topTable(efitData, coef=ncol(voomData$design), number=Inf, sort.by="p") 
}
stopifnot(all(DE_topTable$genes %in% rownames(exprDT)))
stopifnot(all(DE_topTable$genes %in% names(rna_geneList)))

##############################################################################################################################################################
################################################################################################## VERSION DESIGN 2
##############################################################################################################################################################

data.used <- exprDT
groups <- my_group

# skip from here if not RNAseq
y=DGEList(counts=data.used,group=groups, genes=rownames(data.used))
y = calcNormFactors(y) 

class.names <- c(sort(c(cond1, cond2)))

t=factor(groups)
design=model.matrix(~0+t)
colnames(design) = class.names

# skip if not RNAseq
v=voom(y,design,plot=FALSE)
# end skip

# if not RNAseq, use data.used instead of v
fit = lmFit(v,design)

contrast = paste(class.names,collapse="-")
cont.matrix = makeContrasts(mut_wt=contrast,levels=design)

fit2  <- contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

table = topTable(fit2,coef=1,number=nrow(data.used), sort.by="p")

stopifnot(all.equal(table, DE_topTable))

cat("*** OK **** \n -> no difference found \n")