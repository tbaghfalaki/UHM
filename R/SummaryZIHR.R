#' Summary of ZIHR
#'
#' @description
#' Computing a summary of the outputs of the ZIHR function
#'
#' @details
#' It provides a summary of the output of the ZIHR function, including parameter estimations.
#'
#'
#' @param object 	an object inheriting from class ZIHR
#'
#'
#' @return
#' Estimation list of posterior summary includes estimation, standard deviation, lower and upper bounds for 95% credible intervals, and Rhat (when n.chain > 1).
#' DIC  deviance information criterion
#' LPML Log Pseudo Marginal Likelihood (LPML) criterion
#'
#'
#'
#' @author Taban Baghfalaki \email{t.baghfalaki@gmail.com}, Mojtaba Ganjali \email{m-ganjali@sbu.ac.ir}
#'
#'
#' @example inst/exampleZIHR.R
#'
#' @md
#' @seealso \code{\link{ZIHR}}
#' @export


SummaryZIHR <- function(object) {
  VV <- object
  Estimation <- VV$Estimation
  DIC <- VV$DIC
  LPML <- VV$LPML
  list(Estimation = Estimation, DIC = DIC, LPML = LPML)
}
