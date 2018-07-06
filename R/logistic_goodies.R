#' Hosmer Lemeshow test
#'
#' Hosmer Lemeshow test for logistic regression
#' @param mod glm (family = binomial) model
#' @param y factor with indicator variable
#' @param yhat predicted probability
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
    list('observed' = obs,
         'expected' = expect,
         'HL test' = c('chisq' = chisq, 'df' = df, 'p.value' = P ))
}


#' Plot logistic model's ROC with AUC and confidence interval
#'
#' Plot logistic model's ROC with AUC and confidence interval
#' @param mod a glm (family = binomial) model
#' @param plot wheter to plot or not the graph
#' @param add add the roc to a previous existing one
#' @param legend_auc add a legend with the AUC and confidence interval
#' @param ... further params passed to plot.roc
#' @export
logistic_roc <- function(mod, plot = TRUE, add = FALSE,
                         legend_auc = TRUE, ...)
{
    db <- data.frame(outcome = mod$y)
    db$pred <- stats::predict(mod, type = 'response')
    ROC <- pROC::roc(outcome ~ pred, data = db)
    AUC_CI <- pROC::ci.auc(ROC)
    if (plot) {
        pROC::plot.roc(ROC, add = add, mar = par("mar"), ...)
        if (legend_auc) {
            leg <- sprintf("AUC = %.3f (%.3f - %.3f)",
                           AUC_CI[2], AUC_CI[1], AUC_CI[3])
            graphics::legend(0.8, 0.2, legend = leg)
        }
    }
    list("roc" = ROC, "auc" = AUC_CI)
}


#' Test added diagnostic accuracy (estimated with AUC) dued to
#' variable insertion in a logistic model
#'
#' The function add one terms per iteration, estimate the model,
#' calculates the roc and finally compare them with z test using
#' delong variance (by default, but other tests can be asked for,
#' using pROC::roc.test).
#'
#' @param formula a glm full formula
#' @param data data.frame
#' @param perm_tests_sequence permutate model covariate sequence (checking
#'     all possible models/sequences given the covariates specified)
#' @param plot_roc plot the overlapping rocs of the models estimated
#' @param col color of the rocs
#' @param lty lty of the rocs
#' @param legend_style covariates_auc for the right hand formula and
#'     auc, numbers_auc for number between parenthesis and auc, none
#'     for avoiding printing legend
#' @param test_style all_combs for testing all the combinations, AvsC
#'     for comparing all the models to the one with only one
#'     covariate, BcsC to compare the models with the one identical
#'     excepting for the last covariate added
#' @param test_param list of parameter to be passed to pROC::roc.test
#' @source idea in Knottnerus Buntinx (editors), Evidence base of
#'     clinical diagnosis, 2008, Wiley Blackwell, pag 157
#' 
#' @export
added_da <- function(formula, data,
                     perm_test_sequence = FALSE,
                     plot_roc = TRUE,
                     col = 'black', 
                     lty = 1,
                     legend_style = c('covariates_auc', 'numbers_auc', 'none'),
                     test_style = c("all_combs", "BvsC", "AvsC"),
                     test_param = list())
{
    ## todo:
    ## - un domani: altro stile legenda che metta in luce la variabile
    ##   incrementale rispetto al modello precedente
    ##   (tipo knotterus a pag 157)
    legend_style <- match.arg(legend_style)
    test_style <- match.arg(test_style)
    all_vars <- all.vars(formula)
    data <- lbmisc::NA_remove(data[, all_vars])

    vars <- if (perm_test_sequence){
                lapply(combinat::permn(all_vars[-1]),
                       function(x) c(all_vars[1], x))
            } else {
                list(all_vars)
            }

    vars_sequences <-  lapply(vars, function(v){
        paste(v[-1], collapse = '->')
    })
    names(vars) <- paste('Seq:', vars_sequences)

    ## now some utility functions
    
    ## make_formulas creates the incremental sequence of formulas given
    ## a specified char with outcome var1 var2 ... varn: that is
    ## "outcome ~ var1", "outcome ~ var1 + var2", and so on till
    ## "outcome ~ var1 + var2 + ... + varn"
    make_formulas <- function(v){# vars is char with 
        outcome <- v[1]
        covariates <- v[-1]
        indexes <- lapply(seq_along(covariates), seq)
        covariates <- lapply(indexes, function(i){
            paste0(covariates[i], collapse = ' + ')
        })
        formulas <- lapply(covariates, function(x){
            formula(paste0(outcome, '~', x))
        })
        names(formulas) <- covariates
        formulas
    }

    ## the roc (estimates/plot) maker given  a logistic model
    roc_maker <- function(m, a, c, l) logistic_roc(mod = m, 
                                                   plot = plot_roc, 
                                                   add = a,
                                                   legend_auc = FALSE,
                                                   col = c,
                                                   lty = l)

    ## roc_selector from the output of logistic_log
    roc_selector <- function(x) x$roc

    ## the function comparing two rocs
    roc_comparator <- function(x, y, xmod, ymod){
        test_param <- c(list("roc1" = x, "roc2" = y), test_param)
        test <- do.call(pROC::roc.test, test_param)
        test$data.names <- paste('model with ', xmod, 
                                 'vs model with', ymod)
        test
    }

    ## given a roc and a prefix
    legend_maker <- function(roc, prefix) {
        auc <- roc$auc
        sprintf("%s AUC: %.2f (%.2f - %.2f)", 
                prefix, auc[2], auc[1], auc[3])
    }
    
    ## this does the estimates and plotting for each sequence of marker
    analysis <- function(var_seq){
        ## get the specific formulas for each test sequence
        formulas <- make_formulas(v = var_seq)
        models <- lapply(formulas, function(f) {
            glm(formula = f, data = data, family = binomial)
        })
    
        ## add only if it's not the first graph
        add <- c(FALSE, rep(TRUE, length(models) - 1))
        rocs <- Map(roc_maker, models,
                    as.list(add), as.list(col), as.list(lty))

        ## tests comparing models
        if (test_style == 'all_combs'){
            all_combs_indexes <- combn(seq_along(rocs), 2)
            base_index        <- all_combs_indexes[1, ]
            comparison_index  <- all_combs_indexes[2, ]
        } else if (test_style == "AvsC"){
            base_index       <- 1
            comparison_index <- -1
        } else if (test_style == "BvsC") {
            base_index       <- -length(rocs)
            comparison_index <- -1
        }
        else {
            stop("here test_style should be one of AvsC or BvsC")
        }
        
        test_base       <- lapply(rocs[base_index], roc_selector)
        test_comparison <- lapply(rocs[comparison_index], roc_selector)
        tests <- Map(roc_comparator, test_base, test_comparison,
                     as.list(names(test_base)),
                     as.list(names(test_comparison)))
        names(tests) <- paste(names(test_base), 'vs', names(test_comparison))

        ## plotting
        if (plot_roc && legend_style != 'none'){
            if (legend_style == 'covariates_auc'){
                prefixes <- paste(names(rocs), ' - ')
            } else if (legend_style == 'numbers_auc'){
                prefixes <- paste0('(', seq_along(rocs), ')')
            } else 
                stop("legend_style here should be one ",
                     "of covariates_auc or numbers_auc")
            legends <- unlist(Map(legend_maker, rocs, as.list(prefixes)))
            ## legend right alignment here, see ?legend example
            longest_legend <- legends[which.max(nchar(legends))]
                                        # browser()
            temp <- legend("bottomright", 
                           legend = rep(" ", length(rocs)),
                           ## qui metto il reverse perché l'asse x è invertito!
                           text.width = - strwidth(longest_legend),
                           bty = 'n', #do not draw the box
                           col = col, 
                           lty = lty, 
                           xjust = 1, yjust = 1)
            text(temp$rect$left + temp$rect$w, temp$text$y,
                 legends, 
                 pos = 2)
        }
        invisible(list("rocs" = rocs, "tests" = tests))

    }

    rval <- lapply(vars, analysis)

    ## if only one sequence tested, simplify by default
    if (length(rval) == 1) rval[[1]] else rval

}
