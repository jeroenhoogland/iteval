#' Calibration model for predicted ITEs
#'
#' Fits a calibration model based on predicted risk under both
#' treatments and an observed binary outcome. Note that the ITE predictions
#' that are evaluated should be for those that were treated (ind.B).
#'
#' @param y0hat Predicted outcome under control treatment.
#' @param y1hat Predicted outcome under treatment.
#' @param y Observed outcome of interest.
#' @param ind.B Index for individuals in the treated group (should be a vector
#' of length y with the indices that belong to treated individuals).
#'
#' @return Calibration model for predicted ITEs
#' @export
#'
#' @examples
#' # generate some  data
#' n <- 250
#' x <- rnorm(n) # covariate
#' trt <- rbinom(n, 1, .5) # treatment
#' beta <- c(0, 1, 1, -.75) # coefficients
#' lp <- cbind(1, x, trt, x*trt) %*% beta # linear predictor
#' y <- rbinom(n, 1, plogis(lp)) # simulated y
#'
#' # model
#' mod <- glm(y ~ x * trt, family = "binomial")
#'
#' # predict ite's
#' control.data <- cbind(1, x, 0, 0)
#' treated.data <- cbind(1, x, 1, x)
#' y0hat <- plogis(control.data %*% coef(mod))
#' y1hat <- plogis(treated.data %*% coef(mod))
#'
#' ind.B <- which(trt == 1)
#' cal(y0hat,y1hat,y,ind.B) # apparent
cal <- function(y0hat,y1hat,y,ind.B=NULL){
  if(is.null(ind.B)){
    ind.B <- 1:length(y1hat)
  }
  stats::glm(y[ind.B] ~ I(stats::qlogis(y1hat[ind.B]) - stats::qlogis(y0hat[ind.B])),
      family="quasibinomial", offset=stats::qlogis(y0hat[ind.B]))
}
