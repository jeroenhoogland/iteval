#' Calibration model for predicted ITEs
#'
#' Fits a calibration model based on predicted risk under both
#' treatments and an observed binary outcome. Note that the ITE predictions
#' that are evaluated should be for those that were treated (ind.B).
#'
#' @param y0hat Predicted outcome under control treatment.
#' @param y1hat Predicted outcome under treatment.
#' @param y Observed outcome of interest.
#' @param trt Vector with treatment allocations (1 for treated, 0 for control).
#' @param p0 Estimate of the potential outcome under the control condition based on independent
#'  data (or known truth in a simulation setting).
#' @param p1 Estimate of the potential outcome under the treated condition based on independent
#' data (or known truth in a simulation setting).
#' @param type Default type="both", indicating that the ITE predictions are evaluated in the
#' control arm based on an offset (the estimated outcome under the treated condition)
#' estimated in the treated group and given as p1, and in the treated arm based on an offset
#' (the estimated outcome under the control condition) as estimated in the control group and
#' given as p0. Types "control" and "treated" only evaluate in this names arm.
#' @return Calibration model for predicted ITEs
#' @export
#'
#' @examples
#' # generate some  data
#' set.seed(123)
#' n <- 500
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
#' # truth probabilities for current sim
#' p0 <- plogis(control.data %*% c(0, 1, 1, -.75))
#' p1 <- plogis(treated.data %*% c(0, 1, 1, -.75))
#'
#' cal(y0hat,y1hat,y,trt,p0,p1,type="both")
cal <- function(y0hat,y1hat,y,trt,p0=NULL,p1=NULL,type="both"){
  if(!type %in% c("treated", "control", "both")){
    stop("Supply type")
  }
  Z <- ifelse(trt==1,1,-1)
  deltahatlp <- stats::qlogis(y1hat) - stats::qlogis(y0hat)
  deltahatlpstar <- ifelse(trt==1, deltahatlp, -deltahatlp)
  if(is.null(p0) & is.null(p1)){
    if(type=="treated"){
      ind.B <- which(trt == 1)
      mod <- stats::glm(y[ind.B] ~ deltahatlp[ind.B],
                 family="quasibinomial", offset= stats::qlogis(y0hat[ind.B]))
    } else if(type=="control"){
      ind.A <- which(trt == 0)
      minusOne <- rep(-1, length(ind.A))
      mod <- stats::glm(y[ind.A] ~ 0 + minusOne + I(-deltahatlp[ind.A]),
                 family="quasibinomial", offset= stats::qlogis(y1hat[ind.A]))
    } else if(type=="both"){
      os <- ifelse(trt==1, stats::qlogis(y0hat), stats::qlogis(y1hat))
      mod <- stats::glm(y ~ 0 + Z + deltahatlpstar,
                 family="quasibinomial", offset=os)
    }
  } else if(!is.null(p0) & !is.null(p1)){
    if(type=="treated"){
      ind.B <- which(trt == 1)
      mod <- stats::glm(p1[ind.B] ~ deltahatlp[ind.B],
                 family="quasibinomial", offset= stats::qlogis(p0[ind.B]))
    } else if(type=="control"){
      ind.A <- which(trt == 0)
      minusOne <- rep(-1, length(ind.A))
      mod <- stats::glm(p0[ind.A] ~ 0 + minusOne + I(-deltahatlp[ind.A]),
                 family="quasibinomial", offset= stats::qlogis(p1[ind.A]))
    } else if(type=="both"){
      os <- ifelse(trt==1, stats::qlogis(p0), stats::qlogis(p1))
      p01 <- ifelse(trt==1,p1,p0)
      mod <- stats::glm(p01 ~ 0 + Z + deltahatlpstar,
                 family="quasibinomial", offset=os)
    } else {
      stop("Supply both p0 and p1 (for estimand under simulation) or neither (empirical)")
    }
  }
  return(mod)
}
