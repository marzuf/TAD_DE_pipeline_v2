startTime <- Sys.time()

nPerm <- 100
nCpu <- 30


suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
registerDoMC(nCpu)

# TEST1: library(GOSemSim) here -> does not work
# suppressPackageStartupMessages(library(GOSemSim, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 


combineSemSimMethod="BMA"
semSimMetric = "Resnik"
randomGenes <- c("835", "5261","241")

# foo <- foreach(curr_perm = seq_len(nPerm), .combine='c') %dopar% {
#   suppressPackageStartupMessages(library(GOSemSim, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
#   hs <- godata('org.Hs.eg.db', ont="BP", computeIC = TRUE) # works 10 nPerm 2 CPUM 50 nPerm and 30 nCPU
# 
#   random_semSim <- mgeneSim(genes=randomGenes,
#                             # semData=hs, # works 10 nPerm 2 CPUM 50 nPerm and 30 nCPU
#                             semData = godata('org.Hs.eg.db', ont="BP", computeIC = TRUE), # works 10 nPerm 2 CPUM 50 nPerm and 30 nCPU
#                             combine=combineSemSimMethod,
#                             measure=semSimMetric,
#                             verbose=FALSE)
# 
# }
foo <- lapply(c(1:10), function(x) {
  foo <- foreach(curr_perm = seq_len(nPerm), .combine='c') %dopar% {
    suppressPackageStartupMessages(library(GOSemSim, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
    hs <- godata('org.Hs.eg.db', ont="BP", computeIC = TRUE) # works 10 nPerm 2 CPUM 50 nPerm and 30 nCPU
    
    random_semSim <- mgeneSim(genes=randomGenes,
                              # semData=hs, # works 10 nPerm 2 CPUM 50 nPerm and 30 nCPU
                              semData = godata('org.Hs.eg.db', ont="BP", computeIC = TRUE), # works 10 nPerm 2 CPUM 50 nPerm and 30 nCPU
                              combine=combineSemSimMethod,
                              measure=semSimMetric,
                              verbose=FALSE)
    
  }
})

cat(paste0(startTime, "\n", Sys.time(), "\n"))