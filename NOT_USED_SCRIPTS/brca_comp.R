gse_data <- eval(parse(text = load("OUTPUT_FOLDER/GSE58135_ERpos_tripleNeg/11_runEmpPvalCombined/emp_pval_combined.Rdata")))
tcga_data <- eval(parse(text = load("OUTPUT_FOLDER/TCGA_brca_lum_bas//11_runEmpPvalCombined/emp_pval_combined.Rdata")))


commonReg <- intersect(names(gse_data), names(tcga_data))

svg("comp_pval_combined_BRCA_GSE58135_TCGA.svg")
plot(log10(gse_data[commonReg]), log10(tcga_data[commonReg]),
     cex = 0.7, pch = 16,
     xlab="emp. p-val combined GSE58135_ERpos_tripleNeg",
     ylab = "emp. p-val combined TCGA_brca_lum_bas")
foo <- dev.off()



random_data <- eval(parse(text = load("OUTPUT_FOLDER/TCGA_stad_msi_gs//11_runEmpPvalCombined/emp_pval_combined.Rdata")))
svg("comp_pval_combined_random_TCGAstad_TCGAbrca.svg")
commonReg2 <- intersect(names(random_data), names(tcga_data))
plot(log10(random_data[commonReg2]), log10(tcga_data[commonReg2]),
     cex = 0.7, pch = 16,
     xlab="emp. p-val combined TCGA_stad_msi_gs",
     ylab = "emp. p-val combined TCGA_brca_lum_bas")
