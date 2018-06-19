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
#' @export
da <- function(test                = NULL,
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


#' compare diagnostic accuracy measures for 2 binary tests evaluated
#' in a paired setting (each test evaluated each patients) having a
#' reference tandard
#'
#' @param test1 oldest test available, logical value with FALSE = non
#'     diseased and TRUE diseased or factor with two levels
#' @param test2 new test, same coding as test1
#' @param test1_lab test 1 label/name
#' @param test2_lab test 2 label/name
#' @param refstd reference standard diagnosis same coding as test1
#' @param alpha Type I error
#' @param boot_R bootstrap repetition
#' @param boot_parallel use parallel 
#' @param boot_ncpus number of cpus dedicated to bootstrapping
#' @examples
#' 
#' db  <- dadb(tp = 70, fp = 20, fn = 30, tn = 80)
#' db2 <- dadb(tp = 50, fp = 45, fn = 50, tn = 55)
#' db$test2 <- db2$test
#' ## test2 is worst than test
#' acc <- da_compare(test1 = db$test, test2 = db$test2, refstd = db$refstd)
#' 
#' @export
da_compare <- function(test1 = NULL, test2 = NULL, refstd = NULL,
                       test1_lab     = 'test1',
                       test2_lab     = 'test2',
                       alpha         = 0.05,
                       boot_R        = 10000,
                       boot_parallel = 'multicore',
                       boot_ncpus    = 8L)
{
    
    test1  <- binary_preproc(test1,  "test1")
    test2  <- binary_preproc(test2,  "test2")
    refstd <- binary_preproc(refstd, "refstd")

    db <- lbmisc::NA_remove(data.frame(test1, test2, refstd))
    nunits <- nrow(db)

    ## single tests performances
    da_test1 <- da(test = db$test1, refstd = db$refstd, alpha = alpha)
    da_test2 <- da(test = db$test2, refstd = db$refstd, alpha = alpha)
    
    group <- function(t, r, postfix){

        tposlev <- levels(t)[2]
        rposlev <- levels(r)[2]

        tpos <- t == tposlev
        tneg <- !tpos
        rpos <- r == rposlev
        rneg <- !rpos

        group <-
            1L * (tpos & rpos) + # true positive
            2L * (tpos & rneg) + # false positive
            3L * (tneg & rpos) + # false negative
            4L * (tneg & rneg)   # true negative

        if (!all(group %in% 1:4))
            stop("strangely some units were not",
                 "classified as TP, FP, FN or TN")

        ## a quale gruppo appartiene il paziente
        group <- c("tp", "fp", "fn", "tn")[group]
        ## il test ha portato ad una diagnosi corretta
        correct <- factor(group %in% c("tp", "tn"),
                          levels = c(TRUE, FALSE),
                          labels = c('Correct', 'Uncorrect'))
        rval <- data.frame(group = group, correct = correct)
        names(rval) <- paste(names(rval), postfix, sep = "_")
        rval
    }

    db <- cbind(db, group(t = db$test1, r = db$refstd, postfix = "1"))
    db <- cbind(db, group(t = db$test2, r = db$refstd, postfix = "2"))
    
    positive     <- db$refstd == levels(db$refstd)[2]
    negative     <- db$refstd == levels(db$refstd)[1]
    diseased     <- c("tp", "fn")
    not_diseased <- c("tn", "fp")

    ## db at this point is like this  
    
    ##   test1 test2 refstd group_1 correct_1 group_2 correct_2
    ## 1     -     - absent      tn   Correct      tn   Correct
    ## 2     -     - absent      tn   Correct      tn   Correct
    ## 3     -     - absent      tn   Correct      tn   Correct
    ## 4     -     - absent      tn   Correct      tn   Correct

    ## ----------------------------------------
    ## tables for bonett_price (Sens, Spec, Acc
    ## ----------------------------------------
    sens_tab <- with(db[positive, ], {
        ## per mantenere costante la struttura indipendentemente dai
        ## dati specifico i livelli
        table(factor(group_1, levels = diseased),
              factor(group_2, levels = diseased))
    })

    spec_tab <- with(db[negative, ], {
        table(factor(group_1, levels = not_diseased),
              factor(group_2, levels = not_diseased))
    })

    acc_tab <- with(db, table(correct_1, correct_2))

    tabs <- list("Sensitivity" = sens_tab,
                 "Specificity" = spec_tab,
                 "Accuracy"    = acc_tab)
    cis <- lapply(tabs, function(x)
        ## qui li devo fare reversed per avere test2 - test1 (invece
        ## che test1 - test2)
        bonett_price(n12 = x[2,1],
                     n21 = x[1,2],
                     N = sum(x),
                     alpha = alpha))
    cis <- do.call(rbind, cis)

    mcnemar <- lapply(tabs,function(t) stats::mcnemar.test(x = t)$p.value)
    mcnemar <- unlist(mcnemar)
                      
    ## -----------------------
    ## PPV and NPV (bootstrap
    ## -----------------------
    ## http://www.statmethods.net/advstats/bootstrapping.html
    ## https://en.wikipedia.org/wiki/Bootstrapping_(statistics)
    
    pv <- function(x){# x is a vector of tp, fp, fn, tn
        ppv <- sum(x %in% 'tp')/sum(x %in% c('tp', 'fp'))
        npv <- sum(x %in% 'tn')/sum(x %in% c('tn', 'fn'))
        c("ppv" = ppv, "npv" = npv)
    }

    pv2 <- pv(db$group_2)
    pv1 <- pv(db$group_1)
    pv_diff <- pv2 - pv1

    boot_pv <- function(data, indexes){
        d <- data[indexes, ]
        pv2 <- pv(d$group_2)
        pv1 <- pv(d$group_1)
        pv2 - pv1
    }
    
    pv_boot <- boot::boot(data = db,
                          statistic = boot_pv,
                          R = boot_R,
                          parallel = boot_parallel,
                          ncpus = boot_ncpus)

    extract <- function(x) x$bca[1, 4:5]
    ppv_diff_ci <- extract(boot::boot.ci(pv_boot, type = 'bca', index = 1)) ## ppv
    npv_diff_ci <- extract(boot::boot.ci(pv_boot, type = 'bca', index = 2)) ## npv

    ## add estimate and set names
    ppv_diff_ci <- setNames(c(pv_diff["ppv"], ppv_diff_ci), c('Difference', 'Lower.CI', 'Upper.CI'))
    npv_diff_ci <- setNames(c(pv_diff["npv"], npv_diff_ci), c('Difference', 'Lower.CI', 'Upper.CI'))
    pvs <- do.call(rbind, list('PPV' = ppv_diff_ci, 'NPV' = npv_diff_ci))
    cis <- rbind(cis, pvs)
                              
    ## -----------------------
    ## Return
    ## -----------------------
    rval <- list("test1" = da_test1, "test2" = da_test2,
                 "diffs" = cis, "mcnemar" = mcnemar)
    names(rval)[1:2] <- c(test1_lab, test2_lab)
    rval
}


#' Confidence interval for the difference between two
#' paired proportions in a 2x2 table
#'
#' This function computes the Wald confidence interval with
#' Bonett-Price Laplace adjustment confidence interval for the
#' difference between two paired proportions. It gives the estimate as
#' well.
#' 
#' @param n12 corresponding cell frequency in 2x2 setting
#' @param n21 corresponding cell frequency in 2x2 setting
#' @param N table total frequency
#' @param alpha type I error
#'
#' @references
#'
#' Fagerland MW, Lydersen S, Laake P. Recommended tests and confidence
#'     intervals for paired binomial proportions.  Stat Med. 2014 Jul
#'     20;33(16):2850-75. doi: 10.1002/sim.6148. Epub 2014 Mar 20.
#' 
#' @examples
#' 
#' # table V fagerland etal, recommended tests and confindence intervals
#' # for paired binomial proportions  ,stat in med 2014
#' 
#' bonett_price(n12 = 1, n21 = 7, N =21)
#' @export
bonett_price <- function(n12, n21, N, alpha = 0.05){
    p12 <- (n12 + 1)/(N+2)
    p21 <- (n21 + 1)/(N+2)
    pdiff <- p12 - p21
    z <- qnorm(1 - alpha/2)
    se <- sqrt( (p12 + p21 - (pdiff)^2) / (N+2))
    est <- (n12 - n21) / N
    lower <- max(pdiff - z * se, -1L)
    upper <- min(pdiff + z * se, +1L)
    c('Difference' = est, 'Low.CI' = lower, 'Up.CI' = upper)
}
