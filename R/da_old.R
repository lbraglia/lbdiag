#' Calculate diagnostic accuracy measures for binary measures (test, reference
#' standard)
#' 
#' Calculate diagnostic accuracy measures for binary measures (test,
#' reference standard) and confidence intervals (two sided binomial
#' for all the indexes except PPV and NPV, calculated via bdpv a-la
#' Mercaldo)
#' 
#' @param test test variable (dichotomic factor or logical)
#' @param refstd reference standard (dichotomic factor or logical)
#' @param alpha type I error for two sided confidence interval
#' @param digits rounding digits
#' @param positive_first logical: display positive test and reference
#'     standard in first row column
#' @param ppv_npv_prev prevalence adopted for ppv and npv confidence
#'     interval (if NULL it's estimated from sample)
#' @param ppv_npv_force_unadj logical: force unadjusted ppv npv estimates (and
#'     confidence interval)
#' @examples
#' 
#' ## CASS Example (Pepe pag 22)
#' db <- dadb(tn = 327, fn = 208, fp = 115, tp = 815)
#' with(db, da(test = test, refstd = refstd))
#' 
#' ## Alzheimer Example (Mercaldo 2007)
#' db <- dadb(tn = 288, fn = 178 , fp = 87, tp = 240)
#' with(db, da(test = test, refstd = refstd, ppv_npv_prev = .03))
#' 
da_old <- function(test                = NULL,
                   refstd              = NULL,
                   alpha               = 0.05,
                   digits              = 4L,
                   positive_first      = TRUE,
                   ppv_npv_prev        = NULL,
                   ppv_npv_force_unadj = FALSE)
{

    ## Test input
    test <- binary_preproc(test, "test")
    refstd <- binary_preproc(refstd, "refstd")

    ## Test parametri rimanenti
    ok.input <-
        (is.logical(positive_first)) &
        (is.null(ppv_npv_prev) || is.numeric(ppv_npv_prev))
    
    ## testa l'immissione
    if (!ok.input)
        stop("OR positive_first not logical\n",
             "OR 'ppv_npv_prev'  not (null|numeric)")

    
    ## ##############   TABLE ##########################
    
    ## Gestione delle label della tabella per la stampa
    tab <- table(test, refstd)
    tn <- tab[1, 1]
    tp <- tab[2, 2]
    fn <- tab[1, 2]
    fp <- tab[2, 1]
    
    ## positive first handling
    ## ora che ho estratto i dati da una forma tipica dell'elaborazione 
    ## li metto in forma tipica per la stampa, se desiderato
    dn <- list("test" = rev(levels(test)), "refstd" = rev(levels(refstd)))
    tab.positive_first <- matrix(c(tab[2,2], tab[2,1], 
                                   tab[1,2], tab[1,1]),
                                 ncol = 2,
                                 byrow = TRUE,
                                 dimnames = dn)
    class(tab.positive_first) <- class(tab)

    ## Table for results depends on positive_first
    tab.for.results <- if (positive_first) {
        tab.positive_first
    } else {
        tab
    }

    ## ##############   STATS   ##########################
	
    n.diseased <- sum(tab[,2])
    n.non.diseased <- sum(tab[,1])
    n.positive <- sum(tab[2,])
    n.negative <- sum(tab[1,])
    n <- sum(tab)
    
    ## per il calcolo degli intervalli di confidenza di ppv npv
    z <- qnorm(1 - alpha/2)

    ## prevalence
    prevalence <- n.diseased / n
    prevalence.ci <- (binom.test( n.diseased, n, conf.level = 1 - alpha )$conf.int)[1:2]

    ### Sensitivity
    sensitivity <- tp/n.diseased
    sensitivity.ci <- (binom.test( tp, n.diseased, conf.level = 1 - alpha )$conf.int)[1:2]
    
    ## Specificity
    specificity <- tn/n.non.diseased
    specificity.ci <- (binom.test( tn, n.non.diseased, conf.level = 1 - alpha )$conf.int)[1:2]

    ## Accuracy
    accuracy    <- (tn+tp)/ n
    accuracy.ci <- (binom.test( tn+tp, n, conf.level = 1 - alpha )$conf.int)[1:2]
    
    
    ## Odd ratio diagnostico
    ## or <- (sensitivity/(1-sensitivity))/ ((1-specificity)/specificity)
    or <- (tp * tn)/ (fp*fn)
    or.se <- sqrt(1/tp + 1/tn + 1/fp + 1/fn)
    or.ci <- exp( c(log(or) - z*or.se, log(or) + z*or.se) )
    
    ## Youden index
    youden <- sensitivity + specificity - 1
    
    
    ## Predictive value: for more stuff check out library(bdpv)
    ## -------------------------
    ## Intervalli di confidenza di npv npp adottando la formula
    ## di mercaldo (pag 108 di zhou e stat in medicine 2007) che
    ## permette di specificare una prevalenza diversa da quella
    ## riscontrata nel campione.  Tutte e due le implementazioni
    ## del paper disponibili.  Prima la standard logit, qui
    ## sotto, poi la adjusted logit

    if( is.null(ppv_npv_prev) ) {
        ## Se l'utente non ha specificato un valore di prevalenza, 
        ## prendere quello del campione
        prev <- prevalence
    } else if (is.numeric(ppv_npv_prev) & (ppv_npv_prev>=0 &
                                           ppv_npv_prev<=1 )){
        ## Se invece l'ha specificato, posto che sia ragionevo
        prev <- ppv_npv_prev
    } else {
        stop("ppv_npv_prev must be NULL or numeric [0,1]")		
    }	
	
    if (all(tab>0) | ( any(tab==0) & (ppv_npv_force_unadj==TRUE)) ) {
        ## PPV NPV standard logit estimates (mercaldo)
        
        ## PPV standard logit
        ppv <- (sensitivity*prev)/((sensitivity*prev) + (1-specificity)*(1-prev) )
        logit.ppv <- log((sensitivity * prev) /
                         ((1-specificity)*(1-prev)))
        var.logit.ppv <- ((1-sensitivity)/sensitivity)*(1/n.diseased) +
            (specificity/(1-specificity))*(1/n.non.diseased)

        ppv.wald.low <- exp(logit.ppv - z*sqrt(var.logit.ppv))/
            (1 + exp(logit.ppv - z*sqrt(var.logit.ppv)))
        ppv.wald.up <- exp(logit.ppv + z*sqrt(var.logit.ppv))/
            (1 + exp(logit.ppv + z*sqrt(var.logit.ppv)))
        ppv.ci <- c(ppv.wald.low, ppv.wald.up)

        ## NPV standard logit
        npv <- 	(specificity*(1- prev))/
            (  (1-sensitivity)*(prev)  +  specificity*(1- prev) )
        
        logit.npv <- log((specificity * (1-prev)) /
                         ((1-sensitivity)*(prev)))
        var.logit.npv <- (sensitivity/(1-sensitivity))*(1/n.diseased) +
            ((1-specificity)/specificity)*(1/n.non.diseased)

        npv.wald.low <- exp(logit.npv - z*sqrt(var.logit.npv))/
            (1 + exp(logit.npv - z*sqrt(var.logit.npv)))
        npv.wald.up <- exp(logit.npv + z*sqrt(var.logit.npv))/
            (1 + exp(logit.npv + z*sqrt(var.logit.npv)))
        npv.ci <- c(npv.wald.low, npv.wald.up)
        
    } else {
	
        ## PPV NPV adjusted logit estimates (mercaldo)
        xmat <- unname(as.matrix(tab.positive_first))
        class(xmat) <- "matrix"
        pv.tmp <- bdpv::BDtest(xmat = xmat, 
                         pr=prev, 
                         conf.level = 1 - alpha)[["PPVNPVDAT"]]
        ppv.row <- row.names(pv.tmp) %in% "PPV"
        npv.row <- row.names(pv.tmp) %in% "NPV"
        ppv <- pv.tmp[ppv.row,1]
        ppv.ci <-  unlist(unname(c(pv.tmp[ppv.row, 3:4])))
        npv <- pv.tmp[npv.row,1]
        npv.ci <- unlist(unname(c(pv.tmp[npv.row, 3:4])))
	
    }
	
    ## Statistiche complessive
    stat <- as.data.frame( rbind(
        c(prevalence, prevalence.ci),
        c(sensitivity, sensitivity.ci),
        c(specificity, specificity.ci),
        c(accuracy, accuracy.ci),
        c(ppv, ppv.ci),
        c(npv, npv.ci),
        c(or, or.ci),
        c(youden, rep(NA,2))
	))
    
    stat <- round(stat,digits)
    names(stat) <- c("Est","Low.CI", "Up.CI")
    rownames(stat) <- c("Prevalence","Sensitivity",
                        "Specificity","Accuracy",
                        "PPV","NPV","OR","Youden")
    
    ## #####################   OUTPUT   ##########################
		
    ## return results
    list("table"= addmargins(tab.for.results), 
         "stats"= stat)
	
}
