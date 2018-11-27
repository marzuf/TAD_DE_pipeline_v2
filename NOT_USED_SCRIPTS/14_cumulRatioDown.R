startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script8_name <- "8_runRatioDown"
script_name <- "14_cumulRatioDown"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
#source("run_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(Hmisc, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar

registerDoMC(ifelse(SSHFS,2, nCpu)) # loaded from main_settings.R

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))


#****************************************************************************************************************************

# INPUT DATA
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)


# geneList -> entrezID, names original rownames
# UNCOMMENT ONE OF THE TWO !!!!
#if(useFilterData){
#    ### VERSION 1: script1_name => use only the genes used in DE
#    geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_geneList.Rdata"))))
#    txt <- paste0(toupper(script_name), "> Use filtered rnaseqDT from DE\n")
#    printAndLog(txt, pipLogFile)
#} else {
#    ### VERSION 0: script0_name => use all genes
#    geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_geneList.Rdata"))))
#    txt <- paste0(toupper(script_name), "> Use unfiltered rnaseqDT\n")
#    printAndLog(txt, pipLogFile)
#}

# UPDATE SELECT THE GENES ACCORDING TO THE SETTINGS PREPARED IN 0_PREPGENEDATA
initList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_geneList.Rdata"))))
geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))
txt <- paste0(toupper(script_name), "> Start with # genes: ", length(geneList), "/", length(initList), "\n")
printAndLog(txt, pipLogFile)

stopifnot(!any(duplicated(names(geneList))))

gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% geneList,]
geneNbr <- setNames(as.numeric(table(gene2tadDT$region)), names(table(gene2tadDT$region)))


obs_ratioDown <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/all_obs_ratioDown.Rdata"))))
permut_ratioDown <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/ratioDown_permDT.Rdata"))))

# ensure I used the same set of TADs for the permutation and for the calculation
# (NB: would also be possible to filter the obs_ratioDown, but not the permut_ratioDown)
stopifnot(all(names(obs_ratioDown) %in% rownames(permut_ratioDown)))
stopifnot(all(rownames(permut_ratioDown) %in% names(obs_ratioDown)))

interReg <- intersect(names(obs_ratioDown),rownames(permut_ratioDown) )

### SET OUTPUT
plotType <- "svg"
myHeight <- ifelse(plotType == "png", 480 , 7)
myWidth <- ifelse(plotType == "png", 600, 10)


outFileRank <- paste0(curr_outFold, "/biggerSmallerAtSameRank.", plotType)

outFile1_A <- paste0(curr_outFold, "/ecdf_obs_permut_ratio0_1.", plotType)
outFile1_B <- paste0(curr_outFold, "/ecdf_obs_permut_ratio05_1.", plotType)

outFile2_A <- paste0(curr_outFold, "/cumsum_obs_permut_ratio0_1.", plotType)
outFile2_B <- paste0(curr_outFold, "/cumsum_obs_permut_ratio05_1.", plotType)

outFile3_A <- paste0(curr_outFold, "/sigfit_ecdf_obs_permut_ratio0_1.", plotType)
outFile3_B <- paste0(curr_outFold, "/sigfit_obs_permut_ratio05_1.", plotType)

outFile4_A <- paste0(curr_outFold, "/sigfit_allperm_ecdf_obs_permut_ratio0_1.", plotType)
outFile4_B <- paste0(curr_outFold, "/sigfit_allperm_obs_permut_ratio05_1.", plotType)

###################################################################################################
# get the direction of up/down
###################################################################################################
# retrieve the direction of up/down
DE_topTable <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_topTable.Rdata"))))
DE_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_geneList.Rdata"))))
exprDT <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_rnaseqDT.Rdata"))))
samp1 <- eval(parse(text=load(paste0(setDir, "/", sample1_file))))
samp2 <- eval(parse(text=load(paste0(setDir, "/", sample2_file))))
stopifnot(all(DE_topTable$genes %in% names(DE_geneList)))
stopifnot(!any(duplicated(names(DE_geneList))))
stopifnot(all(colnames(exprDT) %in% c(samp1, samp2)))
stopifnot(all(samp1 %in% colnames(exprDT)))
stopifnot(all(samp2 %in% colnames(exprDT)))
maxDownGene <- DE_topTable$genes[which.min(DE_topTable$logFC)]
stopifnot(maxDownGene %in% rownames(exprDT))
mean_expr1 <- mean(unlist(c(exprDT[maxDownGene, samp1])), na.rm=T)
mean_expr2 <- mean(unlist(c(exprDT[maxDownGene, samp2])), na.rm=T)

if(mean_expr1 > mean_expr2) {
  subtitDir <- paste0("down: ", toupper(cond1), " > ", toupper(cond2))
} else{
  subtitDir <- paste0("down: ", toupper(cond2), " > ", toupper(cond1))
}

############################################################
############################################################ # filter the TADs and sort
############################################################
filter_obs_ratioDown <- sort(obs_ratioDown[interReg], decreasing = T)

filter_permut_ratioDown_unsort <- permut_ratioDown[interReg,]
stopifnot(length(filter_obs_ratioDown) == nrow(filter_permut_ratioDown_unsort))

filter_permut_ratioDown <- apply(filter_permut_ratioDown_unsort, 2, sort, decreasing=T)
rownames(filter_permut_ratioDown) <- NULL
stopifnot(length(filter_obs_ratioDown) == nrow(filter_permut_ratioDown_unsort))

# change so that the ratioDown ranges between 0.5 and 1 (-> e.g. treats 0.1 as 0.9)
filter_obs_ratioDown_half <- abs(filter_obs_ratioDown - 0.5) + 0.5
filter_permut_ratioDown_half <- abs(filter_permut_ratioDown - 0.5) + 0.5


############################################################
############################################################ # rank at the position
############################################################

ii <- cut(filter_obs_ratioDown, breaks = seq(min(filter_obs_ratioDown), max(filter_obs_ratioDown), len = 100),
          include.lowest = TRUE)
colors <- colorRampPalette(c("blue","white", "red"))(99)[ii]

howmanyMoreCoord <- foreach(i = 1:length(filter_obs_ratioDown), .combine='c') %dopar% {
  obsRatio <- filter_obs_ratioDown[i]
  if(obsRatio > 0.5){
    sum(filter_permut_ratioDown[i,] >= obsRatio)
  } else if(obsRatio < 0.5) {
    sum(filter_permut_ratioDown[i,] <= obsRatio)
  } else{
    NA
  }
}

do.call(plotType, list(outFileRank, height=myHeight, width=myWidth))
plot(howmanyMoreCoord ~ c(1:length(filter_obs_ratioDown)), 
     xlab = "position", 
     col=colors,
     ylab="# of perm. at same position with <> ratio down",
     bty="l", pch=16)
mtext(subtitDir, font=3)
foo <- dev.off()
cat(paste0("... written: ", outFileRank, "\n"))

#stop("ok")


############################################################
############################################################ # ecdf
############################################################
do.call(plotType, list(outFile1_A, height=myHeight, width=myWidth))
plot_ecdf(filter_obs_ratioDown, filter_permut_ratioDown)
mtext(subtitDir, font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile1_A, "\n"))

do.call(plotType, list(outFile1_B, height=myHeight, width=myWidth))
plot_ecdf(filter_obs_ratioDown_half, filter_permut_ratioDown_half, halfOnly = T)
mtext(subtitDir, font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile1_B, "\n"))

############################################################
############################################################ # difference to 50
############################################################

do.call(plotType, list(outFile2_A, height=myHeight, width=myWidth))
plot_cumsumDiff05(filter_obs_ratioDown, filter_permut_ratioDown)
mtext(subtitDir, font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile2_A, "\n"))

do.call(plotType, list(outFile2_B, height=myHeight, width=myWidth))
plot_cumsumDiff05(filter_obs_ratioDown_half, filter_permut_ratioDown_half)
mtext(subtitDir, font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile2_B, "\n"))



############################################################
############################################################ # ecdf sigmoid fit
############################################################
do.call(plotType, list(outFile3_A, height=myHeight, width=myWidth))
plot_fittedLogReg(filter_obs_ratioDown, filter_permut_ratioDown)
mtext(subtitDir, font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile3_A, "\n"))

do.call(plotType, list(outFile3_B, height=myHeight, width=myWidth))
plot_fittedLogReg(filter_obs_ratioDown_half, filter_permut_ratioDown_half, halfOnly = TRUE)
mtext(subtitDir, font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile3_B, "\n"))


############################################################
############################################################ # ecdf sigmoid fit, mean of all perm fit
############################################################

do.call(plotType, list(outFile4_A, height=myHeight, width=myWidth))
plot_fittedLogReg_meanPerm(filter_obs_ratioDown, filter_permut_ratioDown)
mtext(subtitDir, font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile4_A, "\n"))

do.call(plotType, list(outFile4_B, height=myHeight, width=myWidth))
plot_fittedLogReg_meanPerm(filter_obs_ratioDown_half, filter_permut_ratioDown_half, halfOnly = TRUE)
mtext(subtitDir, font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile4_B, "\n"))



txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))


