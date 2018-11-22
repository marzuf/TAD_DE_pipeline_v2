curr_data <- "GSE102073"
curr_data <- "GSE90749"

setDir <- "/media/electron"
geneClass1 <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_DI/", curr_data, "_geneAggregExpression.Rdata_notfpkm1"))))
geneClass2 <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_DI/", curr_data, "_geneAggregExpression.Rdata_notfpkm2"))))
geneClass3 <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_DI/", curr_data, "_geneAggregExpression.Rdata_notfpkm3"))))

all.equal(geneClass1, geneClass2)
all.equal(geneClass2, geneClass3)


fpkm_geneClass1 <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_DI/",curr_data, "_geneAggregExpression.Rdata_fpkm1"))))
fpkm_geneClass2 <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_DI/", curr_data, "_geneAggregExpression.Rdata_fpkm2"))))
fpkm_geneClass3 <- eval(parse(text = load(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_DI/", curr_data, "_geneAggregExpression.Rdata_fpkm3"))))

all.equal(fpkm_geneClass1, fpkm_geneClass2)
all.equal(fpkm_geneClass2, fpkm_geneClass3)

class1_notfpkm <- geneClass1$gene[geneClass1$class == 1]
class1_fpkm <- fpkm_geneClass1$gene[fpkm_geneClass1$class == 1]

cat("> class1, notFPKM in FPKM: ", round(sum(class1_notfpkm %in% class1_fpkm)/length(class1_notfpkm)*100, 2), "% /n")
cat("> class1, FPKM in notFPKM: ", round(sum(class1_fpkm %in% class1_notfpkm)/length(class1_fpkm)*100, 2), "% /n")


class2_notfpkm <- geneClass2$gene[geneClass1$class == 2]
class2_fpkm <- fpkm_geneClass1$gene[fpkm_geneClass1$class == 2]

cat("> class2, notFPKM in FPKM: ", round(sum(class2_notfpkm %in% class2_fpkm)/length(class2_notfpkm)*100, 2), "% /n")
cat("> class2, FPKM in notFPKM: ", round(sum(class2_fpkm %in% class2_notfpkm)/length(class2_fpkm)*100, 2), "% /n")


class3_notfpkm <- geneClass1$gene[geneClass1$class == 3]
class3_fpkm <- fpkm_geneClass1$gene[fpkm_geneClass1$class == 3]

cat("> class3, notFPKM in FPKM: ", round(sum(class3_notfpkm %in% class3_fpkm)/length(class3_notfpkm)*100, 2), "% /n")
cat("> class3, FPKM in notFPKM: ", round(sum(class3_fpkm %in% class3_notfpkm)/length(class3_fpkm)*100, 2), "% /n")


class4_notfpkm <- geneClass1$gene[geneClass1$class == 4]
class4_fpkm <- fpkm_geneClass1$gene[fpkm_geneClass1$class == 4]

cat("> class4, notFPKM in FPKM: ", round(sum(class4_notfpkm %in% class4_fpkm)/length(class4_notfpkm)*100, 2), "% /n")
cat("> class4, FPKM in notFPKM: ", round(sum(class4_fpkm %in% class4_notfpkm)/length(class4_fpkm)*100, 2), "% /n")


class5_notfpkm <- geneClass1$gene[geneClass1$class == 5]
class5_fpkm <- fpkm_geneClass1$gene[fpkm_geneClass1$class == 5]

cat("> class5, notFPKM in FPKM: ", round(sum(class5_notfpkm %in% class5_fpkm)/length(class5_notfpkm)*100, 2), "% /n")
cat("> class5, FPKM in notFPKM: ", round(sum(class5_fpkm %in% class5_notfpkm)/length(class5_fpkm)*100, 2), "% /n")
