#' A quick ROC plot.
#' 
#' 
#' A quick ROC plot. It's a simple wrapper around \code{pROC::roc},
#' which adds Youden's criterion optimal cutoff and conditional
#' positioning for AUC (+CI) estimates.
#' 
#' 
#' @param test Quantitative marker
#' @param refstd Reference standard (binary variable)
#' @param roc_params list of parameters passed to pROC::roc
#' @param plot Logical. Plot ROC curve?
#' @param plot_params list of parameters passed to pROC::plot.roc
#' @param plot_auc add auc estimates to the plot
#' @param auc_cex cex for auc estimates in the plot
#' @return The function plot the graph and return a list with ROC statistics
#' @export
quick_roc <-
    function(test = NULL,
             refstd = NULL,
             roc_params = list(direction = c("auto", "<", ">")),
             plot = TRUE,
             plot_params = list(
                 print.thres = TRUE,
                 mar = par("mar"),
                 print.thres.pattern = "Thresh: %.2f\nSp: %.2f\nSe: %.2f",
                 main = ""),
             plot_auc = TRUE,
             auc_cex = 1)
{

    ## Estimates
    roc_params <- c(list(response = refstd, predictor = test), roc_params)
    my.roc <- do.call(pROC::roc, roc_params)
    uc.roc <- unclass(my.roc)
    uc.roc$youden <- with(uc.roc, sensitivities + specificities - 1)
    uc.roc$best.thresh <- uc.roc$thresholds[which.max(uc.roc$youden)] 
    my.ci.auc <- as.numeric(as.character(pROC::ci.auc(my.roc)))
    uc.roc$auc.ci <- my.ci.auc[c(1, 3)]

    ## Plotting
    if (plot) {
        ## Roc plotting
        plot_params <- c(list(x = my.roc), plot_params)
        do.call(pROC::plot.roc, plot_params)
        ## AUC estimates with sensible positioning
        if (plot_auc){
            y_auc_base <- 0.15
            x_auc_base <- 0.30
            y_ci_base  <- 0.1
            x_ci_base  <- 0.30
            if (my.ci.auc[2] > 0.5) {
                y_auc <- y_auc_base
                x_auc <- x_auc_base
                y_ci  <- y_ci_base 
                x_ci  <- x_ci_base 
            } else {
                y_auc <- 1 - y_auc_base
                x_auc <- 1 - x_auc_base
                y_ci  <- 1 - y_ci_base 
                x_ci  <- 1 - x_ci_base 
            }
            text(x = x_auc, y = y_auc, "AUC: ", pos = 2, cex = auc_cex)
            text(x = x_ci,  y = y_ci, "95% CI: ", pos = 2 , cex = auc_cex)
            text(x = x_auc, y = y_auc, sprintf("%.2f", my.ci.auc[2]),
                 pos = 4, cex = auc_cex)
            text(x = x_ci,  y = y_ci,
                 sprintf("%.2f - %.2f", my.ci.auc[1], my.ci.auc[3]),
                 cex = auc_cex, pos=4)
        }
    }
    uc.roc
}

#' ROC graph with boxplot below
#'
#' ROC graph with boxplot below
#' 
#' @param test Marker
#' @param refstd Reference standard (binary variable)
#' @param direction Charachter. Direction passed to \code{pROC::roc}
#' @param layout_heights layout heights for roc/boxplot image proportions
#' @param oma oma parameter
#' @param mar mar parameter
#' @export
roc_with_boxplot <- function(test, refstd,
                             direction = c("auto", "<", ">"),
                             layout_heights = c(3, 1),
                             oma = c(0, 0, 0, 0),
                             mar = c(4, 6, 0, 0) + 0.1
                             )
{
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    direction <- match.arg(direction)
    layout(mat = matrix(1:2), heights = layout_heights)
    par(oma = oma, mar = mar)
    roc <- quick_roc(test = test, refstd = refstd,
                     roc_params = list(direction = direction))
    cutoff <- roc$best.thresh
    bp <- boxplot(test ~ refstd,
                  ylab = "",
                  xlab = comment(test),
                  horizontal = TRUE,
                  las = 1)
    abline(v = cutoff, col = 'red')
    invisible(list("roc" = roc, "bp" = bp, "cutoff" = cutoff))
}
