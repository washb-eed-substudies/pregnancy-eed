# 
# #--------------------------------------------------------------
# # get bootstrap pvalue
# #--------------------------------------------------------------
# #https://stats.stackexchange.com/questions/20701/computing-p-value-using-bootstrap-with-r
# 
# library(bootstrap)
# require(boot)
# 
# ratio <- function(d, w) sum(d$x * w)/sum(d$u * w)
# city.boot <- boot(city, ratio, R = 999, stype = "w", sim = "ordinary")
# boot.ci(city.boot, conf = c(0.90, 0.95),
#         type = c("norm", "basic", "perc", "bca"))
# city.boot$t0
# 
# #The 95% CI for the normal bootstrap is obtained by calculating:
# with(city.boot, 2*t0 - mean(t) + qnorm(c(0.025, 0.975)) %o% sqrt(var(t)[1,1]))
# with(city.boot, pnorm(abs((2*t0 - mean(t) - 1) / sqrt(var(t)[1,1])), lower.tail=F)*2)
# 
# #quantile based
# quantile(city.boot$t, c(0.025, 0.975))
# cvs <- quantile(city.boot$t0 - city.boot$t + 1, c(0.025, 0.975))
# mean(city.boot$t > cvs[1] & city.boot$t < cvs[2])





#--------------------------------------------------------------
# analysis functions
#--------------------------------------------------------------

# data=d
# Y="ln_L_conc_t1"
# X="ln_preg_cort"
# pair = NULL 
# W = NULL #pick_covariates(Y)
# forcedW = NULL
# V = NULL
# id = "clusterid"
# family = "gaussian"
# pval = 0.2
# print = TRUE
# verbose = FALSE
# B = 200


washb_glm_lasso_boot <- function(data, Y, X, pair = NULL, W, forcedW = NULL, V = NULL, id = "clusterid", 
          family = "gaussian", pval = 0.2, print = TRUE, 
          verbose = FALSE, B = 200){
  
  require(tidyverse)
  Subgroups = NULL
  if (!is.null(W)) {
    colnamesW <- W
    W <- data %>% select(all_of(W)) %>% data.frame()
  }

    if(!is.null(W)) {
      glmdat <- data.frame(Y=data[[Y]], X=data[[X]], id=data[[id]], W)
    }else{
      colnamesW=NULL
      glmdat <- data.frame(Y=data[[Y]], X=data[[X]], id=data[[id]])
    }


  n.orig <- dim(glmdat)[1]
  rowdropped <- rep(1, nrow(glmdat))
  rowdropped[which(complete.cases(glmdat))] <- 0
  glmdat <- glmdat[complete.cases(glmdat), ]
  n.sub <- dim(glmdat)[1]
  if (print == TRUE) 
    if (n.orig > n.sub) 
      cat("\n-----------------------------------------\nDropping", 
          n.orig - n.sub, "observations due to missing values in 1 or more variables\n", 
          "Final sample size:", n.sub, "\n-----------------------------------------\n")
  
  # if (!is.null(W)) {
  #   colnamesW <- names(W)
  # }
  # if (!is.null(W)) {
  #   if (!is.null(V)) {
  #     forcedW = c(V, forcedW)
  #   }
  #   if (!is.null(forcedW)) {
  #     screenW <- subset(glmdat, select = colnamesW)
  #     toexclude <- names(screenW) %in% forcedW
  #     if (length(which(toexclude == TRUE)) != length(forcedW)) 
  #       stop("A forcedW variable name is not a variable within the W data frame.")
  #     screenW = screenW[!toexclude]
  #     if (ncol(screenW) == 0) {
  #       screenW <- NULL
  #     }
  #     if (print == TRUE) {
  #       cat("\n-----------------------------------------\nInclude the following adjustment covariates without screening:\n-----------------------------------------\n")
  #       print(forcedW, sep = "\n")
  #     }
  #   }else {
  #     screenW <- subset(glmdat, select = colnamesW)
  #   }
  # }else {
  #   screenW <- NULL
  # }
  # colnamesW <- colnames(screenW)
  

    suppressWarnings(fit <- cowboy_glm(data = glmdat, clusterid = "id", 
                                       Ws = colnamesW, forcedW = forcedW, pair = NULL, family = family, 
                                       B = B, confint.level = 0.95, n.cores = 1))
    # fit$boot.coefs
    # fit$boot.pval
    #NOTE: need to verify the p-value calculation for binary outcomes, so not reporting yet
    
    if (family == "gaussian") {
      modelfit <- data.frame(est = fit$boot.coefs[2], ci.lb = fit$percentile.interval[2, 
                                                                                      1], ci.ub = fit$percentile.interval[2, 2], se = fit$boot.sds[2], pval= fit$boot.pval)
    }else{
      modelfit <- data.frame(RR = exp(fit$boot.coefs[2]), 
                             est = fit$boot.coefs[2], ci.lb = exp(fit$percentile.interval[2, 
                                                                                          1]), ci.ub = exp(fit$percentile.interval[2, 
                                                                                                                                   2]), se = fit$boot.sds[2])
    }
  
  return(list(TR = modelfit, fit = fit))
}



# data = glmdat
# clusterid = "id"
# Ws = colnamesW
# pair = NULL
# confint.level = 0.95
# n.cores = 1
# forcedW = NULL
# with.replacement = T
# family = "gaussian"

cowboy_glm <- function (data, clusterid = "id", Ws, forcedW = NULL, pair = NULL, 
          with.replacement = T, family = "gaussian", B = 200, confint.level = 0.95, 
          n.cores = 8) 
{
  require(rsample)
  set.seed(12345)
  bfull <- paste(c("X", Ws, forcedW, pair), collapse = "+")
  if (!is.null(Ws)) {
    prescreened_Ws <- washb_glmnet_prescreen(Y = data$Y, 
                                             data %>% select(!!(Ws)), family = family)
  }else {
    prescreened_Ws = NULL
  }
  b <- paste(c("X", prescreened_Ws, forcedW, pair), collapse = "+")
  full_model <- as.formula(paste("Y ~ ", bfull, sep = ""))
  model <- as.formula(paste("Y ~ ", b, sep = ""))
  res.or <- glm(full_model, family = family, data = data)
  res.or.screen <- glm(model, family = family, data = data)
  confint.pboundaries = c((1 - confint.level)/2, 1 - (1 - confint.level)/2)
  confint.Zboundaries = qnorm(confint.pboundaries)
  n <- nrow(data)
  p <- length(res.or$coef)
  coefs <- matrix(NA, nrow = B, ncol = p)
  if (with.replacement) {
    f = NA
    D <- data %>% as_tibble() %>% nest(-id)
    bs <- bootstraps(D, times = B)
    for (i in 1:B) {
      set.seed(i)
      dboot <- as.tibble(bs$splits[[i]]) %>% arrange(id) %>% 
        unnest(cols = c(data))
      if (!is.null(Ws)) {
        prescreened_Ws <- washb_glmnet_prescreen(Y = dboot$Y, 
                                                 dboot %>% select(!!(Ws)), family = family)
      }
      else {
        prescreened_Ws = NULL
      }
      b <- paste(c("X", prescreened_Ws, forcedW, pair), 
                 collapse = "+")
      model <- as.formula(paste("Y ~ ", b, sep = ""))
      bootcoef <- tryCatch(coef(glm(model, family = family, 
                                    data = dboot)), error = function(x) rep(as.numeric(NA), 
                                                                            p))
      coefs[i, which(names(res.or$coef) %in% names(bootcoef))] <- bootcoef
    }
  }else{
    cluster <- as.character(data[[clusterid]])
    clusters <- unique(data[[clusterid]])
    Obsno <- split(1:n, cluster)
    f = matrix(clusters, length(clusters), B)
    ff = matrix(f, prod(dim(f)), 1)
    fff = sample(ff)
    f = matrix(fff, length(clusters), B)
    for (i in 1:B) {
      set.seed(i)
      j <- f[, i]
      obs <- unlist(Obsno[j])
      dboot = data[obs, ]
      table(dboot$Y)
      if (!is.null(Ws)) {
        prescreened_Ws <- washb_glmnet_prescreen(Y = dboot$Y, 
                                                 dboot %>% select(!!(Ws)), family = family)
      }
      else {
        prescreened_Ws = NULL
      }
      b <- paste(c("X", prescreened_Ws, forcedW, pair), 
                 collapse = "+")
      model <- as.formula(paste("Y ~ ", b, sep = ""))
      bootcoef <- tryCatch(coef(glm(model, family = family, 
                                    data = dboot)), error = function(x) rep(as.numeric(NA), 
                                                                            p))
      coefs[i, which(names(res.or$coef) %in% names(bootcoef))] <- bootcoef
    }
  }
  
  invalid.samples <- colSums(is.na(coefs))
  names(invalid.samples) <- colnames(coefs) <- names(res.or$coef)
  samples.with.NA.coef <- which(is.na(rowSums(coefs)))
  sdcoefs <- apply(coefs, 2, sd, na.rm = TRUE)
  ci_percentile <- t(apply(coefs, 2, quantile, probs = confint.pboundaries, 
                           na.rm = TRUE))
  ci_parametric <- cbind(res.or$coef + confint.Zboundaries[1] * 
                           sdcoefs, res.or$coef + confint.Zboundaries[2] * sdcoefs)
  ci_BCa <- matrix(NA, 1, 1)
  rownames(ci_percentile) <- dimnames(ci_parametric)[[1]]
  colnames(ci_parametric) <- dimnames(ci_percentile)[[2]]
  
  t0=res.or$coef[2]
  t = coefs[,2]
  b.under.H0 <- t - mean(t)
  boot.pval <- mean(abs(b.under.H0) > abs(t0))
  
  
  result <- list(call = match.call(), model = model, family = family, 
                 B = B, coefficients = coefs, data = data, bootstrap.matrix = f, 
                 subject.vector = clusterid, lm.coefs = res.or$coef, 
                 boot.coefs = colMeans(coefs, na.rm = TRUE), boot.sds = sdcoefs, 
                 boot.pval=boot.pval,ci.level = confint.level, 
                 percentile.interval = ci_percentile, parametric.interval = ci_parametric, 
                 BCa.interval = ci_BCa, samples.with.NA.coef = samples.with.NA.coef, 
                 failed.bootstrap.samples = invalid.samples)
  class(result) <- "clusbootglm"
  return(result)
}



