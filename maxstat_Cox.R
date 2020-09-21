maxstat_Cox <- function(data, formula.cut, formula, seed, ...){
  require(maxstat)
  require(survival)
  
  # cutpoint search using maxstat package  
  cut.var <- labels(terms(formula.cut))
  cut.search <- maxstat.test(formula.cut, data=data, smethod="LogRank",...)

  # create binary variable
  cutpoint <- cut.search$estimate
  data[, paste(cut.var, "bin", sep=".")] <- factor(data[,cut.var]>cutpoint, levels=c(FALSE, TRUE), labels = c("low", "high")) #add dichotomous variable "cut.var.bin"
  
  # update formula for multivariable model
#  terms.updated <- c(labels(terms(formula))[!(labels(terms(formula))==cut.var)], paste0("I(",cut.var, ">",cutpoint,")"))
  terms.updated <- c(labels(terms(formula))[!(labels(terms(formula))==cut.var)], paste(cut.var, "bin", sep=".")) #füge Terme der Formel ein, die nicht die Cutvariable waren und setze sie mit cut.var.bin
  formula <- update.formula(formula, formula(paste("~", paste(terms.updated, collapse="+"))))
  
  # Cox model with dichotomized variable
  mod.multi <- coxph(formula=formula, data=data, model=TRUE, x = TRUE, y= TRUE)
  mod.multi$formula <- formula
  
  # extract estimate and variance of dichotomized variable
  est <- mod.multi$coefficients[grepl(cut.var, names(mod.multi$coefficients))]
  var <- mod.multi$var[grepl(cut.var, names(mod.multi$coefficients)), grepl(cut.var, names(mod.multi$coefficients))]
  
  # compute shrinkage factor (cf. Hollaender et al (2008), page 1703)
  sf <- (est^2 - var)/est^2 
  # adjust effect estimate
  e <- sf*est
  mod.multi$coefficients[grepl(cut.var, names(mod.multi$coefficients))] <- e

  # determine standard error using bootstrap
  B = 100 # number of Bootstrap samples
  eb <- sfb <- NULL
  set.seed(seed)
  for(j in seq(B)) {
    jsam <- sample(seq(nrow(data)), replace=TRUE)
    nbj <- data[jsam,]
    mj <- maxstat.test(formula.cut, data=nbj, ...)$estimate
    nbj[, paste(cut.var, "bin", sep=".")] <- factor(nbj[,cut.var]>mj, levels=c(FALSE, TRUE), labels = c("low", "high"))
    fitj <- coxph(formula, data=nbj)
    eb[j] <- fitj$coef[grepl(cut.var, names(fitj$coef))]   # Estimate
    sfb[j] <- (eb[j]^2-fitj$var[grepl(cut.var, names(fitj$coef)), grepl(cut.var, names(fitj$coef))])/eb[j]^2 # Shrinkage factor
  }
  # compute empirical bootstrap standard deviation
  sdb <- sd(sfb*eb)
  mod.multi$var[grepl(cut.var, names(mod.multi$coef)), grepl(cut.var, names(mod.multi$coef))] <- sdb^2  
  mod.multi
  out <- list(fit = mod.multi, formula = formula, cutpoint = cutpoint)
  out$call <- match.call()
  class(out) <- "maxstat_Cox"
  out
}

predictSurvProb.maxstat_Cox <- function (object, newdata, times, ...) 
{
  predictSurvProb(object[[1]], newdata = newdata, times = times, ...)
}
