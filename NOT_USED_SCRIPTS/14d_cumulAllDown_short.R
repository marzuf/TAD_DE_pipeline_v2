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
script8_name <- "8c_runAllDown"
script_name <- "14d_cumulAllDown_short"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
#source("run_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar

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

### SET OUTPUT
plotType <- "svg"
myHeight <- ifelse(plotType == "png", 480 , 7)
myWidth <- ifelse(plotType == "png", 600, 10)

# "permThresh" the quantile of permutations to take is loaded from main_settings.R

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

# allDown retrieved from main_settings.R


for(curr_ratio_type in allDown) {

    cat(paste0("*** START ", curr_ratio_type, "\n"))
    
    obs_curr_down <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/all_obs_", curr_ratio_type, ".Rdata"))))
    permut_currDown <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/", curr_ratio_type, "_permDT.Rdata"))))
  
    # ensure I used the same set of TADs for the permutation and for the calculation
    # (NB: would also be possible to filter the obs_curr_down, but not the permut_currDown)
    stopifnot(all(names(obs_curr_down) %in% rownames(permut_currDown)))
    stopifnot(all(rownames(permut_currDown) %in% names(obs_curr_down)))

    interReg <- intersect(names(obs_curr_down),rownames(permut_currDown) )

    ############################################################
    ############################################################ # filter the TADs and sort
    ############################################################
    filter_obs_curr_down <- sort(obs_curr_down[interReg], decreasing = T)

    filter_permut_currDown_unsort <- permut_currDown[interReg,]
    stopifnot(length(filter_obs_curr_down) == nrow(filter_permut_currDown_unsort))

    filter_permut_currDown <- apply(filter_permut_currDown_unsort, 2, sort, decreasing=T)
    rownames(filter_permut_currDown) <- NULL
    stopifnot(length(filter_obs_curr_down) == nrow(filter_permut_currDown_unsort))

    # change so that the ratioDown ranges between 0.5 and 1 (-> e.g. treats 0.1 as 0.9)
    filter_obs_curr_down_half <- abs(filter_obs_curr_down - 0.5) + 0.5
    filter_permut_currDown_half <- abs(filter_permut_currDown - 0.5) + 0.5

    ############################################################
    ############################################################ # difference to 50
    ############################################################
    curr_ratio_type_concord <- paste0(curr_ratio_type, "_concord")  
    
    
    ############################### PLOT CUMSUM POLYGON 1) DEPARTURE 0.5, RATIO
    outFile2_A <- paste0(curr_outFold, "/", curr_ratio_type, "departure05_cumsum_obs_permut_ratio0_1.", plotType)
    do.call(plotType, list(outFile2_A, height=myHeight, width=myWidth))
    plot_cumsumDiff05(observ_vect = filter_obs_curr_down, permut_DT = filter_permut_currDown,  
                      my_stat = curr_ratio_type, departureValue = 0.5)
    mtext(subtitDir, font=3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile2_A, "\n"))

    ############################### PLOT CUMSUM POLYGON 2) DEPARTURE 0, RATIO
    outFile2b_A <- paste0(curr_outFold, "/", curr_ratio_type, "departure0_cumsum_obs_permut_ratio0_1.", plotType)
    do.call(plotType, list(outFile2b_A, height=myHeight, width=myWidth))
    plot_cumsumDiff05(observ_vect = filter_obs_curr_down, permut_DT = filter_permut_currDown,  
                      my_stat = curr_ratio_type, departureValue = 0)
    mtext(subtitDir, font=3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile2b_A, "\n"))

    ############################### PLOT CUMSUM POLYGON 3) DEPARTURE 0.5, RATIOCONCORD    
    outFile2_B <- paste0(curr_outFold, "/", curr_ratio_type, "departure05_cumsum_obs_permut_ratio05_1.", plotType)
    do.call(plotType, list(outFile2_B, height=myHeight, width=myWidth))
    plot_cumsumDiff05(filter_obs_curr_down_half, filter_permut_currDown_half, my_stat = curr_ratio_type_concord,
                      departureValue = 0.5)
    mtext(subtitDir, font=3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile2_B, "\n"))

    ############################### PLOT CUMSUM POLYGON 4) DEPARTURE 0.5, RATIOCONCORD    
    outFile2b_B <- paste0(curr_outFold, "/", curr_ratio_type, "departure0_cumsum_obs_permut_ratio05_1.", plotType)
    do.call(plotType, list(outFile2b_B, height=myHeight, width=myWidth))
    plot_cumsumDiff05(filter_obs_curr_down_half, filter_permut_currDown_half, my_stat = curr_ratio_type_concord,
                      departureValue = 0)
    mtext(subtitDir, font=3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile2b_B, "\n"))
    
    ############################### PLOT AUC CUMSUM POLYGON 1) DEPARTURE 0.5, RATIO
    outFile_aucA <- paste0(curr_outFold, "/", curr_ratio_type, "departure05_auc_cumsum_obs_permut_ratio0_1.", plotType)
    do.call(plotType, list(outFile_aucA, height=myHeight, width=myWidth))
    auc_ratio_stat_dep05 <- calc_aucObsPerm_cumsum(observ_vect = filter_obs_curr_down, permut_DT=filter_permut_currDown, 
                               thresh_perm =permThresh, 
                               obsLineCol = "black",
                               permutLineCol = "dodgerblue",
                               polygonPermutCol =  rgb(255/255,102/255,102/255, 0.3),
                               departureValue=0.5, my_stat=curr_ratio_type)
    foo <- dev.off()
    cat(paste0("... written: ", outFile_aucA, "\n"))
    
    
    ############################### PLOT AUC CUMSUM POLYGON 2) DEPARTURE 0, RATIO
    outFile_aucB <- paste0(curr_outFold, "/", curr_ratio_type, "departure0_auc_cumsum_obs_permut_ratio0_1.", plotType)
    do.call(plotType, list(outFile_aucB, height=myHeight, width=myWidth))
    auc_ratio_stat_dep0 <- calc_aucObsPerm_cumsum(observ_vect = filter_obs_curr_down, permut_DT=filter_permut_currDown, 
                                     thresh_perm =permThresh, 
                                     obsLineCol = "black",
                                     permutLineCol = "dodgerblue",
                                     polygonPermutCol =  rgb(255/255,102/255,102/255, 0.3),
                                     departureValue=0, my_stat=curr_ratio_type)
    foo <- dev.off()
    cat(paste0("... written: ", outFile_aucB, "\n"))
    
    ############################### PLOT AUC CUMSUM POLYGON 3) DEPARTURE 0.5, RATIOCONCORD    
    outFile_aucC <- paste0(curr_outFold, "/", curr_ratio_type, "departure05_auc_cumsum_obs_permut_ratio05_1.", plotType)
    do.call(plotType, list(outFile_aucC, height=myHeight, width=myWidth))
    auc_ratio_stat_concord_dep05 <- calc_aucObsPerm_cumsum(observ_vect = filter_obs_curr_down_half, permut_DT=filter_permut_currDown_half, 
                                             thresh_perm =permThresh, 
                                             obsLineCol = "black",
                                             permutLineCol = "dodgerblue",
                                             polygonPermutCol =  rgb(255/255,102/255,102/255, 0.3),
                                             departureValue=0.5, my_stat=curr_ratio_type_concord)
    foo <- dev.off()
    cat(paste0("... written: ", outFile_aucC, "\n"))
    
    ############################### PLOT AUC CUMSUM POLYGON 4) DEPARTURE 0, RATIOCONCORD    
    outFile_aucD <- paste0(curr_outFold, "/", curr_ratio_type, "departure0_auc_cumsum_obs_permut_ratio05_1.", plotType)
    do.call(plotType, list(outFile_aucD, height=myHeight, width=myWidth))
    auc_ratio_stat_concord_dep05 <- calc_aucObsPerm_cumsum(observ_vect = filter_obs_curr_down_half, permut_DT=filter_permut_currDown_half, 
                               thresh_perm =permThresh, 
                               obsLineCol = "black",
                               permutLineCol = "dodgerblue",
                               polygonPermutCol =  rgb(255/255,102/255,102/255, 0.3),
                               departureValue=0, my_stat=curr_ratio_type_concord)
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile_aucD, "\n"))

}


txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))


