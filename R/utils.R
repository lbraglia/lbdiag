binary_preproc <- function(x, varname){
    if (!((is.factor(x) && nlevels(x) == 2L) || is.logical(x))){
        msg <- sprintf("%s must be a logical or a factor with 2 levels", varname)
        stop(msg)
    }
    if (is.logical(x)) factor(x, levels = c(FALSE, TRUE)) else x
}
