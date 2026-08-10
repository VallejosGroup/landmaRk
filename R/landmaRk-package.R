#' landmaRk
## usethis namespace: start
#' @import dplyr
#' @importFrom foreach `%dopar%`
#' @importFrom methods new is show
#' @importFrom prodlim Hist
#' @importFrom stats as.formula formula reformulate sigma predict model.matrix
#' @importFrom survival Surv
#' @importFrom utils capture.output
#' @importFrom utils head
## usethis namespace: end
#' @docType package
#' @name landmaRk
#' @keywords internal
"_PACKAGE"

# `fgwt` is not a global variable; it is looked up as a column of
# `finegray_data` by .fit_finegray_survival()'s call to survival::coxph()
# (as in ?survival::finegray's own example).
utils::globalVariables("fgwt")
