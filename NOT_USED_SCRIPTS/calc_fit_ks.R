

calc_all_ks_scores <- function(currOutFold) { 

    #****************************************************************************************************************************
    # INPUT DATA
    gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
    gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

    # UPDATE SELECT THE GENES ACCORDING TO THE SETTINGS PREPARED IN 0_PREPGENEDATA
    initList <- eval(parse(text = load(paste0(currOutFold, "/", script0_name, "/rna_geneList.Rdata"))))
    geneList <- eval(parse(text = load(paste0(currOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))
    txt <- paste0(toupper(script_name), "> Start with # genes: ", length(geneList), "/", length(initList), "\n")
    cat(txt)
    stopifnot(!any(duplicated(names(geneList))))
    gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% geneList,]
    geneNbr <- setNames(as.numeric(table(gene2tadDT$region)), names(table(gene2tadDT$region)))
    obs_ratioDown <- eval(parse(text = load(paste0(currOutFold, "/", script8_name, "/all_obs_ratioDown.Rdata"))))
    permut_ratioDown <- eval(parse(text = load(paste0(currOutFold, "/", script8_name, "/ratioDown_permDT.Rdata"))))
    # ensure I used the same set of TADs for the permutation and for the calculation
    # (NB: would also be possible to filter the obs_ratioDown, but not the permut_ratioDown)
    stopifnot(all(names(obs_ratioDown) %in% rownames(permut_ratioDown)))
    stopifnot(all(rownames(permut_ratioDown) %in% names(obs_ratioDown)))
    interReg <- intersect(names(obs_ratioDown),rownames(permut_ratioDown) )

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
    ############################################################ # cumSum with permut for ratioDown - KS pval and D stat
    ############################################################

    cumsumDiff05_line_mean_perm_ks <- function(observ_vect, permut_DT) {
      stopifnot(length(observ_vect) == nrow(permut_DT))
      # observed vect
      observ_vect <- sort(observ_vect, decreasing = T)
      obs_cumsum <- cumsum(abs(observ_vect - 0.5))
      # permut vect - v1 first cumsum than means
      permut_DT_sort <- apply(permut_DT, 2, sort, decreasing=T)
      rownames(permut_DT_sort) <- NULL
      permut_DT_cumsum <- apply(permut_DT_sort, 2, function(x) cumsum(abs(x-0.5)))
      permut_cumsum <- rowMeans(permut_DT_cumsum, na.rm=T)
      x_val <- c(1:length(observ_vect))
      # ks test
      obs_perm_ks_test <- ks.test(obs_cumsum, permut_cumsum)
      return(c(cumsum_permMean_ks_pval = as.numeric(obs_perm_ks_test$p.value), cumsum_permMean_ks_D_stat = as.numeric(obs_perm_ks_test$statistic)))
    }
    cumsum_test_mean <- cumsumDiff05_line_mean_perm_ks(observ_vect=filter_obs_ratioDown_half, permut_DT=filter_permut_ratioDown_half)


    # DO THE SAME BUT INSTEAD OF TAKING THE MEAN VALUE OF EACH ROW FOR THE PERMUT, TAKE THE MAX VALUE
    cumsumDiff05_line_max_perm_ks <- function(observ_vect, permut_DT) {
      stopifnot(length(observ_vect) == nrow(permut_DT))
      # observed vect
      observ_vect <- sort(observ_vect, decreasing = T)
      obs_cumsum <- cumsum(abs(observ_vect - 0.5))
      # permut vect - v1 first cumsum than means
      permut_DT_sort <- apply(permut_DT, 2, sort, decreasing=T)
      rownames(permut_DT_sort) <- NULL
      permut_DT_cumsum <- apply(permut_DT_sort, 2, function(x) cumsum(abs(x-0.5)))
      permut_cumsum <- apply(permut_DT_cumsum, 1, max, na.rm=T)
      x_val <- c(1:length(observ_vect))
      # ks test
      obs_perm_ks_test <- ks.test(obs_cumsum, permut_cumsum)
      return(c(cumsum_permMax_ks_pval = as.numeric(obs_perm_ks_test$p.value), cumsum_permMax_ks_D_stat = as.numeric(obs_perm_ks_test$statistic)))
    }
    cumsum_test_max <- cumsumDiff05_line_max_perm_ks(observ_vect=filter_obs_ratioDown_half, permut_DT=filter_permut_ratioDown_half)



    ############################################################
    ############################################################ # ecdf sigmoid fit and ks test KS TEST ON ALL REGIONS SIGMOID FIT - take the means of ranked permuts and fit the sig
    ############################################################
    fittedLogReg_fitMeanPerm_ks <- function(observ_vect, permut_DT) {
      logisticFunc <- function(x, m, b){
        exp(m*x +b)/(1+exp(m*x+b))
      }
      logisticDeriv <- function(x, m, b) {
        (m*exp(m*x+b))/((1+exp(m*x+b))^2)
      }
      logisticGetX <- function(y, m, b) {
        (log(y/(1-y)) - b)/m
      }
      # prepare the OBSERVED data
      observ_vect <- sort(observ_vect, decreasing = T)
      observ_vect_ecdf <- ecdf(observ_vect)(observ_vect)
      # fit the logistic model
      mod_obs <- glm(observ_vect_ecdf ~ observ_vect, family=binomial(link=logit))
      obs_ipt <- as.numeric(coef(mod_obs)["(Intercept)"])
      obs_slope <- as.numeric(coef(mod_obs)[2])
      obs_der05 <- logisticDeriv(0.5, m=obs_slope, b=obs_ipt)
      obs_y05 <- logisticGetX(0.5, m=obs_slope, b=obs_ipt)
      obs_logFit <- fitted(mod_obs)
      # take the means of the permut and fit a single logistic 
      permut_DT_sort <- apply(permut_DT, 2, function(x){sort(x, decreasing=T)})
      permut_DT_vect <- rowMeans(permut_DT_sort, na.rm=T)  
      stopifnot(length(permut_DT_vect) == nrow(permut_DT))
      permut_DT_vect_ecdf <- ecdf(permut_DT_vect)(permut_DT_vect)
      # fit the logistic model
      mean_perm_obs <- glm(permut_DT_vect_ecdf ~ permut_DT_vect, family=binomial(link=logit))
      mean_perm_ipt <- as.numeric(coef(mean_perm_obs)["(Intercept)"])
      mean_perm_slope <- as.numeric(coef(mean_perm_obs)[2])
      mean_perm_der05 <- logisticDeriv(0.5, m=mean_perm_slope, b=mean_perm_ipt)
      mean_perm_y05 <- logisticGetX(0.5, m=mean_perm_slope, b=mean_perm_ipt)
      # v2: to have the same number of points, fit the observ_vect
      mean_perm_logFit <- logisticFunc(observ_vect, m=mean_perm_slope, b=mean_perm_ipt)
      stopifnot(length(mean_perm_logFit) == length(obs_logFit))
      # KS test between the observed and permut
      vect_ks_test <- ks.test(obs_logFit, mean_perm_logFit)
      return(c(obs_slope_meanPerm = obs_slope, meanPerm_slope = mean_perm_slope, fitMeanPerm_ks_pval = as.numeric(vect_ks_test$p.value), fitMeanPerm_ks_D_stat = as.numeric(vect_ks_test$statistic)))
    }

    fitMeanPerm_test <- fittedLogReg_fitMeanPerm_ks(observ_vect = filter_obs_ratioDown, permut_DT = filter_permut_ratioDown)

    ################################################################## KS TEST ON ALL REGIONS SIGMOID FIT - fit all permut then take the mean
    fittedLogReg_fitAll_meanPerm_ks <- function(observ_vect, permut_DT) {
      logisticFunc <- function(x, m, b){
        exp(m*x +b)/(1+exp(m*x+b))
      }
      logisticDeriv <- function(x, m, b) {
        (m*exp(m*x+b))/((1+exp(m*x+b))^2)
      }
      logisticGetX <- function(y, m, b) {
        (log(y/(1-y)) - b)/m
      }
      # prepare the OBSERVED data
      observ_vect <- sort(observ_vect, decreasing = T)
      observ_vect_ecdf <- ecdf(observ_vect)(observ_vect)
      # fit the logistic model
      mod_obs <- glm(observ_vect_ecdf ~ observ_vect, family=binomial(link=logit))
      obs_ipt <- as.numeric(coef(mod_obs)["(Intercept)"])
      obs_slope <- as.numeric(coef(mod_obs)[2])
      obs_der05 <- logisticDeriv(0.5, m=obs_slope, b=obs_ipt)
      obs_y05 <- logisticGetX(0.5, m=obs_slope, b=obs_ipt)
      obs_logFit <- fitted(mod_obs)
      
      # fit logistic for each permut
      mod_permut_list <- apply(permut_DT, 2, function(i_perm_vec){
        # prepare the OBSERVED data
        i_perm_vec <- sort(i_perm_vec, decreasing = T)
        i_perm_vect_ecdf <- ecdf(i_perm_vec)(i_perm_vec)
        # fit the logistic model
        i_perm_obs <- glm(i_perm_vect_ecdf ~ i_perm_vec, family=binomial(link=logit))
        i_perm_ipt <- as.numeric(coef(i_perm_obs)["(Intercept)"])
        i_perm_slope <- as.numeric(coef(i_perm_obs)[2])
        i_perm_der05 <- logisticDeriv(0.5, m=i_perm_slope, b=i_perm_ipt)
        i_perm_y05 <- logisticGetX(0.5, m=i_perm_slope, b=i_perm_ipt)
        c(perm_ipt = i_perm_ipt, perm_slope= i_perm_slope, 
          perm_der05 = i_perm_der05, perm_y05 = i_perm_y05  )
      })
      stopifnot(ncol(mod_permut_list) == ncol(permut_DT))
      stopifnot(nrow(mod_permut_list) == 4 )
      perm_mean_slope <- mean(mod_permut_list["perm_slope",], na.rm =  T)
      perm_meanLogFit <- logisticFunc(observ_vect, m=perm_mean_slope, b=mean(mod_permut_list["perm_ipt",],na.rm=T))
      stopifnot(length(perm_meanLogFit) == length(obs_logFit))
      # KS test between the observed and permut
      vect_ks_test <- ks.test(obs_logFit, perm_meanLogFit)
      return(c(obs_slope_meanAllPerm = obs_slope, meanAllPerm_slope = perm_mean_slope, fit_meanAllPerm_ks_pval = as.numeric(vect_ks_test$p.value), fit_meanAllPerm_ks_D_stat = as.numeric(vect_ks_test$statistic)))
    }

    fitMeanAllPerm_test <- fittedLogReg_fitAll_meanPerm_ks(observ_vect = filter_obs_ratioDown, permut_DT = filter_permut_ratioDown)

    return(c(cumsum_test_mean, cumsum_test_max, fitMeanPerm_test, fitMeanAllPerm_test))

}
#      return(c(cumsum_ks_pval = as.numeric(obs_perm_ks_test$p.value), cumsum_ks_D_stat = as.numeric(obs_perm_ks_test$statistic)))
#      return(c(obs_slope_meanPerm = obs_slope, meanPerm_slope = mean_perm_slope, fitMeanPerm_ks_pval = as.numeric(vect_ks_test$p.value), fitMeanPerm_ks_D_stat = as.numeric(vect_ks_test$statistic)))
#      return(c(obs_slope_meanAllPerm = obs_slope, meanAllPerm_slope = perm_mean_slope, fit_meanAllPerm_ks_pval = as.numeric(vect_ks_test$p.value), fit_meanAllPerm_ks_D_stat = as.numeric(vect_ks_test$statistic)))

