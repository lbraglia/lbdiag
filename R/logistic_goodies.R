#' Hosmer Lemeshow test
#'
#' Hosmer Lemeshow test for logistic regression
#' @param mod glm (family = binomial) model
# #' @param y factor with indicator variable
# #' @param yhat predicted probability
#' @param g number of groups
#' @export
hosmerlem <- function(mod,
                      ## y,
                      ## yhat,
                      g = 10) {
    ## y <- as.integer(y) - 1
    y <- mod$y
    yhat <- stats::predict(mod, type = 'response')
    Breaks <- quantile(yhat, probs = seq(0, 1, 1/g)) 
    cutyhat <- cut(yhat, breaks = Breaks, include.lowest=TRUE)
    obs <- xtabs(cbind(1 - y, y) ~ cutyhat)
    expect <- xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
    chisq <- sum((obs - expect)^2/expect)
    df <- g - 2
    P <- 1 - stats::pchisq(q = chisq, df = df)
    return(list('observed' = obs,
                'expected' = expect,
                'HL test' = c('chisq' = chisq,
                              'df' = df,
                              'p.value' = P )))
}


#' Plot logistic model's ROC with AUC and confidence interval
#'
#' Plot logistic model's ROC with AUC and confidence interval
#' @param mod a glm (family = binomial) model
#' @export
logistic_roc <- function(mod){
    db <- data.frame(outcome = mod$y)
    db$pred <- stats::predict(mod, type = 'response')
    ROC <- pROC::roc(outcome ~ pred, data = db)
    pROC::plot.roc(ROC)
    AUC_CI <- pROC::ci.auc(ROC)
    leg <- sprintf("AUC = %.3f (%.3f - %.3f)",
                   AUC_CI[2], AUC_CI[1], AUC_CI[3])
    graphics::legend(0.8, 0.2, legend = leg)
    invisible(ROC)
}
