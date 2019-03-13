#' Calculate diagnostic accuracy measures for several cutoffs of a quantitative
#' marker.
#' 
#' Calculate diagnostic accuracy measures for several cutoffs of a quantitative
#' marker.
#' 
#' @param cutoffs Cutoffs considered
#' @param direction > or <: if > a test value greater than the cutoff
#'     will be interpreted as a positive test, in < a test value below
#'     the threshold
#' @param test Test
#' @param refstd Reference standard
#' @param digits Rounding digits
#' @param ... parameters passed to da
#' @return A data.frame for diagnostic accuracy studies.
#' @export
cutoff_da <- function(cutoffs = NULL,
                      direction = c(">", "<"),
                      test = NULL,
                      refstd = NULL,
                      digits = 4,
                      ## parameters passed to da
                      ...
                      ){
  res <- list()
  direction <- match.arg(direction)
  for (tres in cutoffs)  {
      index <- which(cutoffs %in% tres)
      test_dummy <- if (direction == ">") (test > tres) else (test < tres)
      tmp <- da(test = test_dummy, refstd = refstd, digits = digits, ...)
      tmp <- tmp[["stats"]]
      tmp$thresh <- tres
      res[[index]] <- tmp[, c("thresh", "Est", "Low.CI", "Up.CI")]
  }

  do.call("rbind", res)
  
  ## Esempio chesi procalcitonina
  ## pct.gram <- roc(gram ~ pct, 
  ## direction="<", 
  ## data = db
  ## )
  ## ## Per verifiche a manina
  ## listRoc <- unclass(pct.gram)
  ## thresh.db <- da.at.cutoff(cutoffs=listRoc$thresholds,
  ## test=db$pct,
  ## refstd=db$gram)
  ## th.spl <- split(thresh.db, thresh.db$stat )
  ## ppv <- th.spl$PPV
  ## npv <- th.spl$NPV

}
