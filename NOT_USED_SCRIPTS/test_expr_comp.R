SSHFS <- T
setDir <- ifelse(SSHFS, "/media/electron", "")

ucecDT <- eval(parse(text = load("/media/electron/mnt/ed4/marie/other_datasets/TCGA_ucec_cnl_msi/29.08_prepData/msi_cnlow_rnaseqDT.Rdata")))
ucec_cnlow <- eval(parse(text = load("/media/electron/mnt/ed4/marie/other_datasets/TCGA_ucec_cnl_msi/29.08_prepData/cnlow_ID.Rdata")))
ucec_msi <- eval(parse(text = load("/media/electron/mnt/ed4/marie/other_datasets/TCGA_ucec_cnl_msi/29.08_prepData/msi_ID.Rdata")))
rm(msi_ID)
rm(cnlow_ID)

crcDT <- eval(parse(text = load("/media/electron/mnt/ed4/marie/scripts/crc_MSI_MSS/28.07_prepData/msi_mss_rnaseqDT.Rdata")))
crc_msi <- eval(parse(text = load("/media/electron/mnt/ed4/marie/scripts/crc_MSI_MSS/28.07_prepData/msi_ID.Rdata")))
crc_mss <- eval(parse(text = load("/media/electron/mnt/ed4/marie/scripts/crc_MSI_MSS/28.07_prepData/mss_ID.Rdata")))
rm(msi_ID)
rm(mss_ID)

gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/v2_assign/genes2tad/all_genes_positions5.txt") 
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

symbolDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_SYMBOL/final_entrez2syno.txt")
symbDT <- read.delim(paste0(symbolDT_file), header=TRUE, stringsAsFactors = F)
symbDT$entrezID <- as.character(symbDT$entrezID)

tad <- "chrX_TAD28"
tad <- "chr1_TAD108"

tad_genes <- gene2tadDT$entrezID[gene2tadDT$region == tad]
tad_symb <- unlist(sapply(tad_genes, function(x)symbDT$symbol[symbDT$entrezID == x][1]))

ucec_entrez <- rownames(ucecDT)[rownames(ucecDT) %in% tad_genes]
sub_ucec <- ucecDT[ucec_entrez,]
# msi_ucec_mean <- mean(rowMeans(sub_ucec[, ucec_msi], na.rm=T))
# # 193.084
# cnlow_ucec_mean <- mean(rowMeans(sub_ucec[, ucec_cnlow], na.rm=T))
# # 215.0973
mean_msi_ucec <- rowMeans(sub_ucec[, ucec_msi], na.rm=T)
mean_cnlow_ucec <- rowMeans(sub_ucec[, ucec_cnlow], na.rm=T)
mean_genes <- names(mean_msi_ucec)
ucec_DT <- data.frame(entrezID = mean_genes, ucec_symb = tad_symb[mean_genes],
                      mean_msi = mean_msi_ucec[mean_genes], mean_cnlow = mean_cnlow_ucec[mean_genes])

crc_entrez <- rownames(crcDT)[rownames(crcDT) %in% tad_genes]
sub_crc <- crcDT[rownames(crcDT) %in% tad_genes,]
# msi_crc_mean <- mean(rowMeans(sub_crc[, crc_msi], na.rm=T))
# # 132.3892
# mss_crc_mean <- mean(rowMeans(sub_crc[, crc_mss], na.rm=T))
# # 208.5103
mean_msi_crc <- rowMeans(sub_crc[, crc_msi], na.rm=T)
mean_mss_crc <- rowMeans(sub_crc[, crc_mss], na.rm=T)
mean_genes <- names(mean_msi_crc)
crc_DT <- data.frame(entrezID = mean_genes, crc_symb = tad_symb[mean_genes],
                      mean_msi = mean_msi_crc[mean_genes], mean_mss = mean_mss_crc[mean_genes])
