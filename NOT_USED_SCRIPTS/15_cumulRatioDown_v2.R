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
script8_name <- "8c_runAllDown"
script_name <- "15_cumulRatioDown_v2"
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

###########################################################################################################################################
######################################################################################################################## START PLOTTING
###########################################################################################################################################

############################################################
############################################################ # plot the ratioDown and ratioUp
############################################################

par(mfrow=c(1,1))
outFileDownUp <- paste0(curr_outFold, "/cumsum_ratioDownUp_sep.", plotType)
do.call(plotType, list(outFileDownUp, height=myHeight, width=myWidth))
plot_cumsum_down_and_up(filter_obs_ratioDown)
foo <- dev.off()
cat(paste0("... written: ", outFileDownUp, "\n"))


############################################################
############################################################ # plot cumSum with permut for ratioDown and ratioUp separately
############################################################
par(mfrow=c(1,1))
outFileRatioConc <- paste0(curr_outFold, "/cumsum_ratioConcord_with_ks_density_funct.", plotType)
do.call(plotType, list(outFileRatioConc, height=myHeight, width=myWidth*3))
plot_cumsumDiff05_line_mean_perm_ks(observ_vect=filter_obs_ratioDown_half, permut_DT=filter_permut_ratioDown_half,
                                    my_xlab = "regions ranked by decreasing ratioConcordant (= abs(ratioDown-0.5))",
                                    my_ylab = "cumsum(abs(ratioConcordant - 0.5))",
                                    my_main = "Cumsum ratioConcordant departure from 0.5")
foo <- dev.off()
cat(paste0("... written: ", outFileRatioConc, "\n"))


############################################################
############################################################ # ecdf sigmoid fit and ks test
############################################################
################################################################## KS TEST ON ALL REGIONS SIGMOID FIT - take the means of ranked permuts and fit the sig
par(mfrow=c(1,1))
outFileECDFsigMeanPermFit <- paste0(curr_outFold, "/ecdf_sigmoid_meanPerm_fit_with_ks_density_funct.", plotType)
do.call(plotType, list(outFileECDFsigMeanPermFit, height=myHeight, width=myWidth*3))
plot_fittedLogReg_fitMeanPerm_ks(observ_vect = filter_obs_ratioDown, permut_DT = filter_permut_ratioDown)
foo <- dev.off()
cat(paste0("... written: ", outFileECDFsigMeanPermFit, "\n"))


################################################################## KS TEST ON ALL REGIONS SIGMOID FIT - fit all permut then take the mean
par(mfrow=c(1,1))
outFileECDFsigAllFit <- paste0(curr_outFold, "/ecdf_sigmoid_allPermFit_mean_with_ks_density_funct.", plotType)
do.call(plotType, list(outFileECDFsigAllFit, height=myHeight, width=myWidth*3))
plot_fittedLogReg_fitAll_meanPerm_ks(observ_vect = filter_obs_ratioDown, permut_DT = filter_permut_ratioDown)
foo <- dev.off()
cat(paste0("... written: ", outFileECDFsigAllFit, "\n"))


############################################################
############################################################ # plot cumSumRatioConcordant with permut, also for ratioDown and ratioUp separately
############################################################

############## prepare observed and permutated data, only those 05_1
filter_obs_ratioDown_05_1 <- filter_obs_ratioDown[filter_obs_ratioDown >= 0.5]
# filter_obs_ratioDown["chr6_TAD1"]  # 0.625
# filter_obs_ratioDown_05_1["chr6_TAD1"] # 0.625
filter_permut_ratioDown_unsort_05_1 <- permut_ratioDown[names(filter_obs_ratioDown_05_1),]
stopifnot(length(filter_obs_ratioDown_05_1) == nrow(filter_permut_ratioDown_unsort_05_1))
filter_permut_ratioDown_05_1 <- apply(filter_permut_ratioDown_unsort_05_1, 2, sort, decreasing=T)
rownames(filter_permut_ratioDown_05_1) <- NULL
# filter_permut_ratioDown_05_1[1:5,1:5]
stopifnot(length(filter_obs_ratioDown_05_1) == nrow(filter_permut_ratioDown_05_1))

filter_obs_ratioDown_0_05 <- filter_obs_ratioDown[filter_obs_ratioDown < 0.5]
# filter_obs_ratioDown["chr5_TAD26"]  # 0.2
# filter_obs_ratioDown_0_05["chr5_TAD26"] # 0.8
# as it is 0-0.5, rank by increasing
filter_obs_ratioDown_0_05 <- abs(filter_obs_ratioDown_0_05 - 1)
filter_obs_ratioDown_0_05 <- sort(filter_obs_ratioDown_0_05, decreasing=T)

filter_permut_ratioDown_unsort_0_05 <- permut_ratioDown[names(filter_obs_ratioDown_0_05),]
stopifnot(length(filter_obs_ratioDown_0_05) == nrow(filter_permut_ratioDown_unsort_0_05))
filter_permut_ratioDown_unsort_0_05 <- abs(filter_permut_ratioDown_unsort_0_05 - 1)
filter_permut_ratioDown_0_05 <- apply(filter_permut_ratioDown_unsort_0_05, 2, sort, decreasing=T)
stopifnot(length(filter_obs_ratioDown_0_05) == nrow(filter_permut_ratioDown_0_05))
# filter_permut_ratioDown_0_05[1:5,1:5]

### PLOT WITH RATIOCONCORDANT INSTEAD OF RATIODOWN
outFileCumSum <- paste0(curr_outFold, "/cumsum_ratioConcord_ratioDown_ratioUp.", plotType)
do.call(plotType, list(outFileCumSum, height=myHeight, width=myWidth*3))
par(mfrow=c(1,3))
plot_cumsumDiff05(filter_obs_ratioDown_half, filter_permut_ratioDown_half,
                  my_xlab = "regions ranked by decreasing ratioConcordant (= abs(ratioDown-0.5))",
                  my_ylab = "cumsum(abs(ratioConcordant - 0.5))",
                  my_main = "Cumsum ratioConcordant departure from 0.5")
### PLOT ONLY RATIODOWN <0.5
# par(mfrow=c(1,1))
plot_cumsumDiff05(filter_obs_ratioDown_0_05, filter_permut_ratioDown_0_05,
                  my_xlab = "regions with ratioDown < 0.5 ranked by decreasing ratioUp",
                  my_ylab = "cumsum(abs(ratioUp - 0.5))",
                  my_main = "Cumsum ratioUp departure from 0.5 (regions with ratioDown < 0.5 only)")
### PLOT ONLY RATIODOWN >=0.5
# par(mfrow=c(1,1))
plot_cumsumDiff05(filter_obs_ratioDown_05_1, filter_permut_ratioDown_05_1,
                  my_xlab = "regions with ratioDown >= 0.5 ranked by decreasing ratioDown",
                  my_ylab = "cumsum(abs(ratioDown - 0.5))",
                  my_main = "Cumsum ratioDown departure from 0.5 (regions with ratioDown >= 0.5 only)")
foo <- dev.off()
cat(paste0("... written: ", outFileCumSum, "\n"))



txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)
cat(paste0("*** DONE: ", script_name, "\n"))

