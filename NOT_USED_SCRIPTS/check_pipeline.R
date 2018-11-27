SSHFS <- T
setDir <- ifelse(SSHFS, "/media/electron", "")

inDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline")

###########################
# test 2_runWilcoxonTAD
###########################
# GSE90749_adipo
stdWilcox <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/iPS_to_adipo_GSE90749/08.08_wilcoxon_meanTAD/all_ttest_meanTAD_qq.Rdata"))))
pipWilcox <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/GSE90749_adipo/2_runWilcoxonTAD/wilcox_ttest_meanTAD_qq.Rdata"))))
interReg <- intersect(names(stdWilcox), names(pipWilcox))
cat(paste0("> 2_runWilcoxonTAD - std: ", length(stdWilcox), "/", length(interReg), "\n"))
cat(paste0("> 2_runWilcoxonTAD - std: ", length(pipWilcox), "/", length(interReg), "\n"))
stdWilcox_pvalWilcox <- unlist(sapply(stdWilcox, function(x) x[["wilcox_pval"]]))
stdWilcox_pvalWilcox <- stdWilcox_pvalWilcox[interReg]
pipWilcox_pvalWilcox <- unlist(sapply(pipWilcox, function(x) x[["wilcox_pval.wilcox_pval"]]))
pipWilcox_pvalWilcox <- pipWilcox_pvalWilcox[interReg]
head(pipWilcox_pvalWilcox)
head(stdWilcox_pvalWilcox)

# crc_msi_mss
stdWilcox <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/crc_MSI_MSS/28.07_wilcoxon_meanTAD/all_ttest_meanTAD_qq.Rdata"))))
pipWilcox <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/TCGA_crc_MSI_MSS/2_runWilcoxonTAD/wilcox_ttest_meanTAD_qq.Rdata"))))
interReg <- intersect(names(stdWilcox), names(pipWilcox))
cat(paste0("> 2_runWilcoxonTAD - std: ", length(stdWilcox), "/", length(interReg), "\n"))
cat(paste0("> 2_runWilcoxonTAD - std: ", length(pipWilcox), "/", length(interReg), "\n"))
stdWilcox_pvalWilcox <- unlist(sapply(stdWilcox, function(x) x[["wilcox_pval"]]))
stdWilcox_pvalWilcox <- stdWilcox_pvalWilcox[interReg]
pipWilcox_pvalWilcox <- unlist(sapply(pipWilcox, function(x) x[["wilcox_pval.wilcox_pval"]]))
pipWilcox_pvalWilcox <- pipWilcox_pvalWilcox[interReg]
head(pipWilcox_pvalWilcox)
head(stdWilcox_pvalWilcox)

###########################
# test 3_runMeanTADLogFC
###########################
# GSE90749_adipo
stdLogFC <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/iPS_to_adipo_GSE90749/08.08_meanLogFC_TAD/all_meanLogFC_TAD.Rdata"))))
pipLogFC <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/GSE90749_adipo/3_runMeanTADLogFC/all_meanLogFC_TAD.Rdata"))))
interReg <- intersect(names(stdLogFC), names(pipLogFC))
cat(paste0("> 3_runMeanTADLogFC - std: ", length(stdLogFC), "/", length(interReg), "\n"))
cat(paste0("> 3_runMeanTADLogFC - std: ", length(pipLogFC), "/", length(interReg), "\n"))
stdLogFC <- stdLogFC[interReg]
pipLogFC <- pipLogFC[interReg]
head(stdLogFC)
head(pipLogFC)

# crc_msi_mss
stdLogFC <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/crc_MSI_MSS/28.07_meanLogFC_TAD/all_meanLogFC_TAD.Rdata"))))
pipLogFC <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/TCGA_crc_MSI_MSS/3_runMeanTADLogFC/all_meanLogFC_TAD.Rdata"))))
interReg <- intersect(names(stdLogFC), names(pipLogFC))
cat(paste0("> 3_runMeanTADLogFC - std: ", length(stdLogFC), "/", length(interReg), "\n"))
cat(paste0("> 3_runMeanTADLogFC - std: ", length(pipLogFC), "/", length(interReg), "\n"))
stdLogFC <- stdLogFC[interReg]
pipLogFC <- pipLogFC[interReg]
head(stdLogFC)
head(pipLogFC)

###########################
# test 4_runMeanTADCorr
###########################
# GSE90749_adipo
stdCorr <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/iPS_to_adipo_GSE90749/08.08_intraCorr_TAD/all_intracorr_TAD_qq.Rdata"))))
pipCorr <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/GSE90749_adipo/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))))
interReg <- intersect(names(stdCorr), names(pipCorr))
cat(paste0("> 4_runMeanTADCorr - std: ", length(stdCorr), "/", length(interReg), "\n"))
cat(paste0("> 4_runMeanTADCorr - std: ", length(pipCorr), "/", length(interReg), "\n"))
stdCorr <- stdCorr[interReg]
pipCorr <- pipCorr[interReg]
head(stdCorr)
head(pipCorr)

# crc_msi_mss
stdCorr <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/crc_MSI_MSS/28.07_intraCorr_TAD/all_intracorr_TAD_qq.Rdata"))))
pipCorr <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/TCGA_crc_MSI_MSS/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))))
interReg <- intersect(names(stdCorr), names(pipCorr))
cat(paste0("> 4_runMeanTADCorr - std: ", length(stdCorr), "/", length(interReg), "\n"))
cat(paste0("> 4_runMeanTADCorr - std: ", length(pipCorr), "/", length(interReg), "\n"))
stdCorr <- stdCorr[interReg]
pipCorr <- pipCorr[interReg]
head(stdCorr)
head(pipCorr)


###########################
# test 8_runRatioDown
###########################
# GSE90749_adipo
stdRatioDown <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/iPS_to_adipo_GSE90749/08.08_ratioDown/obs_ratioDown.Rdata"))))
pipRatioDown <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/GSE90749_adipo/8_runRatioDown/all_obs_ratioDown.Rdata"))))
interReg <- intersect(names(stdRatioDown), names(pipRatioDown))
cat(paste0("> 8_runRatioDown - std: ", length(stdRatioDown), "/", length(interReg), "\n"))
cat(paste0("> 8_runRatioDown - std: ", length(pipRatioDown), "/", length(interReg), "\n"))
stdRatioDown <- stdRatioDown[interReg]
pipRatioDown <- pipRatioDown[interReg]
head(stdRatioDown)
head(pipRatioDown)

# crc_msi_mss
stdRatioDown <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/crc_MSI_MSS/31.07_ratioDown/obs_ratioDown.Rdata"))))
pipRatioDown <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/TCGA_crc_MSI_MSS/8_runRatioDown/all_obs_ratioDown.Rdata"))))
interReg <- intersect(names(stdRatioDown), names(pipRatioDown))
cat(paste0("> 8_runRatioDown - std: ", length(stdRatioDown), "/", length(interReg), "\n"))
cat(paste0("> 8_runRatioDown - std: ", length(pipRatioDown), "/", length(interReg), "\n"))
stdRatioDown <- stdRatioDown[interReg]
pipRatioDown <- pipRatioDown[interReg]
head(stdRatioDown)
head(pipRatioDown)

