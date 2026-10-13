#' Dose calibration of anti-epileptic drugs data
#'
#' @description The SANAD (Standard and New Anti-epileptic Drugs) study (Marson
#'   et al., 2007) is a randomized control trial of standard and new
#'   anti-epileptic drugs, comparing effects on longer term clinical outcomes.
#'   The data consists of longitudinal measurements of calibrated dose for the
#'   groups randomized to a standard drug (CBZ) and a new drug (LTG). The
#'   objective of the analysis is to investigate the effect of drug titration on
#'   the relative effects of LTG and CBZ on treatment failure (withdrawal of the
#'   randomized drug). There are several baseline covariates available, and also
#'   data on the time to withdrawal from randomized drug. This version of the
#'   data has been modified from the joineR package.
#'
#' @usage data(epileptic)
#'
#' @format This is a data frame in the unbalanced format, that is, with one row
#'   per observation. The data consists of columns for patient identifier, time
#'   of measurement, calibrated dose, baseline covariates, and survival data.
#'   The column names are identified as follows:
#'
#'   \describe{
#'
#'   \item{\code{id}}{integer: patient identifier.}
#'
#'   \item{\code{time}}{integer: timing of clinic visit at which dose recorded
#'   (days).}
#'
#'   \item{\code{with.time}}{integer: time of drug withdrawal/maximum follow up
#'   time (days).}
#'
#'
#'   \item{\code{with.status}}{censoring indicator (\code{1 = }withdrawal of
#'   randomized drug and \code{0 = }not withdrawn from randomized drug/lost to
#'   follow up).}
#'
#'   \item{\code{dose}}{numeric: calibrated dose.}
#'
#'   \item{\code{treat}}{factor: randomized treatment (CBZ or LTG).}
#'
#'   \item{\code{age}}{numeric: age of patient at randomization (years).}
#'
#'   \item{\code{gender}}{factor: gender of patient. \code{F = }female, \code{M
#'   = }male.}
#'
#'   \item{\code{learn.dis}}{factor: learning disability.}
#'
#'   }
#' @keywords datasets
#' @source SANAD Trial Group, University of Liverpool. Data modified from the
#'   joineR package.
#' @docType data
#' @seealso [joineR::epileptic], which this dataset was modified from.
#' @references
#'
#' Marson AG, Appleton R, Baker GA, et al. A randomised controlled trial
#' examining longer-term outcomes of standard versus new antiepileptic drugs.
#' The SANAD Trial. \emph{Health Tech Assess}. 2007; \strong{11(37)}.
#'
#' Marson AG, Al-Kharusi AM, Alwaidh M, et al. The SANAD study of effectiveness
#' of carbamazepine, gabapentin, lamotrigine, oxcarbazepine, or topiramate for
#' treatment of partial epilepsy: an unblinded randomised controlled trial.
#' \emph{Lancet}. 2007; \strong{365}: 2007-2013.
#'
#' Williamson PR, Kolamunnage-Dona R, Philipson P, Marson AG. Joint modelling of
#' longitudinal and competing risks data. \emph{Stats Med.} 2008;
#' \strong{27(30)}: 6426-6438.
"epileptic"

#' Advanced ovarian cancer data
#'
#' @description Data from two clinical trials comparing combination
#'   chemotherapies for advanced ovarian cancer (Neijt et al. 1984, 1987).
#'   The data were originally
#'   distributed in the archived \pkg{dynpred} package (as \code{ova}).
#'
#'
#' @usage data(ova)
#'
#' @format A data frame with 358 rows, one per patient:
#'
#'   \describe{
#'
#'   \item{\code{tyears}}{numeric: time in years until death or last
#'   follow-up.}
#'
#'   \item{\code{d}}{survival status; \code{1 = }death, \code{0 = }censored.}
#'
#'   \item{\code{Karn}}{numeric: Karnofsky score, recoded from 6
#'   (\eqn{\le 60}) to 10 (100).}
#'
#'   \item{\code{Broders}}{factor: Broders' histological grade, with levels
#'   \code{"1"}, \code{"2"}, \code{"3"}, \code{"4"}, \code{"unknown"}.}
#'
#'   \item{\code{FIGO}}{factor: FIGO stage, with levels \code{"III"} and
#'   \code{"IV"}.}
#'
#'   \item{\code{Ascites}}{factor: presence of ascites, with levels
#'   \code{"absent"}, \code{"present"}, \code{"unknown"}.}
#'
#'   \item{\code{Diam}}{factor: diameter of the residual tumour, with levels
#'   \code{"micr."}, \code{"<1cm"}, \code{"1-2cm"}, \code{"2-5cm"},
#'   \code{">5cm"}.}
#'
#'   \item{\code{id}}{integer: patient identifier.}
#'
#'   }
#' @keywords datasets
#' @source The archived \pkg{dynpred} package (GPL (>= 2)), companion to
#'   van Houwelingen & Putter (2012).
#' @docType data
#' @references
#'
#' van Houwelingen HC, Putter H (2012). \emph{Dynamic Prediction in Clinical
#' Survival Analysis}. Chapman & Hall.
#'
#' Neijt JP, ten Bokkel Huinink WW, van der Burg MEL, et al. Randomised trial
#' comparing two combination chemotherapy regimens (Hexa-CAF vs CHAP-5) in
#' advanced ovarian carcinoma. \emph{Lancet}. 1984; \strong{2}: 594-600.
#'
#' Neijt JP, ten Bokkel Huinink WW, van der Burg MEL, et al. Randomized trial
#' comparing two combination chemotherapy regimens (CHAP-5 vs CP) in advanced
#' ovarian carcinoma. \emph{J Clin Oncol}. 1987; \strong{5}: 1157-1168.
"ova"

#' Chronic myeloid leukaemia data: static covariates and survival outcome
#'
#' @description Survival outcome and baseline (static) risk factors for 190
#'   patients from the Benelux Chronic Myeloid Leukaemia (CML) trial
#'   (Kluin-Nelemans et al. 1998).
#'   The data were originally distributed in the archived \pkg{dynpred}
#'   package (as \code{wbc1}).
#'
#' @usage data(cml)
#'
#' @format A data frame with 190 rows, one per patient:
#'
#'   \describe{
#'
#'   \item{\code{id}}{integer: patient identifier, shared with
#'   \code{\link{cml_wbc}}.}
#'
#'   \item{\code{tyears}}{numeric: time in years until death or last
#'   follow-up.}
#'
#'   \item{\code{d}}{survival status; \code{1 = }death, \code{0 = }censored.}
#'
#'   \item{\code{sokal}}{numeric: Sokal score, a clinical risk index based on
#'   spleen size, percentage of circulating blasts, platelet count and age at
#'   diagnosis.}
#'
#'   \item{\code{age}}{numeric: age at diagnosis (years).}
#'
#'   }
#' @keywords datasets
#' @source The archived \pkg{dynpred} package (GPL (>= 2)), companion to
#'   van Houwelingen & Putter (2012).
#' @docType data
#' @seealso \code{\link{cml_wbc}} for the companion longitudinal biomarker
#'   data.
#' @references
#'
#' van Houwelingen HC, Putter H (2012). \emph{Dynamic Prediction in Clinical
#' Survival Analysis}. Chapman & Hall.
#'
#' Kluin-Nelemans JC, Delannoy A, Louwagie A, le Cessie S, Hermans J, van der
#' Burgh JF, Hagemeijer AM, van den Berghe H, Benelux CML Study Group (1998).
#' Randomized study on hydroxyurea alone versus hydroxyurea combined with
#' low-dose interferon-alpha 2b for chronic myeloid leukemia. \emph{Blood}.
#' 1998; \strong{91}: 2713-2721.
"cml"

#' Chronic myeloid leukaemia data: longitudinal white blood cell count
#'
#' @description Repeated white blood cell count (WBC) measurements for the
#'   190 patients of \code{\link{cml}}, recorded every one to two months
#'   during the first stage of the Benelux Chronic Myeloid Leukaemia (CML)
#'   trial (Kluin-Nelemans et al. 1998). The data were originally distributed in the archived
#'   \pkg{dynpred} package (as \code{wbc2}).
#'
#' @usage data(cml_wbc)
#'
#' @format A data frame with 5600 rows, in long format (one row per
#'   measurement):
#'
#'   \describe{
#'
#'   \item{\code{id}}{integer: patient identifier, shared with
#'   \code{\link{cml}}.}
#'
#'   \item{\code{time}}{numeric: time in years from diagnosis at which the
#'   measurement was taken.}
#'
#'   \item{\code{lwbc}}{numeric: WBC, transformed and centered as
#'   \eqn{{}^{10}\log(\mathrm{WBC}) - 0.95}.}
#'
#'   }
#' @keywords datasets
#' @source The archived \pkg{dynpred} package (GPL (>= 2)), companion to
#'   van Houwelingen & Putter (2012).
#' @docType data
#' @seealso \code{\link{cml}} for the companion static covariates and
#'   survival outcome.
#' @references
#'
#' van Houwelingen HC, Putter H (2012). \emph{Dynamic Prediction in Clinical
#' Survival Analysis}. Chapman & Hall.
#'
#' Kluin-Nelemans JC, Delannoy A, Louwagie A, le Cessie S, Hermans J, van der
#' Burgh JF, Hagemeijer AM, van den Berghe H, Benelux CML Study Group (1998).
#' Randomized study on hydroxyurea alone versus hydroxyurea combined with
#' low-dose interferon-alpha 2b for chronic myeloid leukemia. \emph{Blood}.
#' 1998; \strong{91}: 2713-2721.
"cml_wbc"
