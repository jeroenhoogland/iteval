#' Model-based concordance for benefit (mbcb)
#'
#' @param y0hat Predicted outcome under control treatment.
#' @param y1hat Predicted outcome under treatment.
#' @param y0hat.updated Predicted outcome under control treatment based on an
#' updated model (e.g., based on validation data)
#' @param y1hat.updated Predicted outcome under treatment based on an
#' updated model (e.g., based on validation data)
#'
#' @return The model-based concordance statistic for predicted ITE
#' @export
#'
#' @examples
#' #' # generate some  data
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
#' # mbcb
#' mbcb(y0hat,y1hat) # apparent mbcb
#'
#' # generate external data
#' n <- 500
#' x <- rnorm(n) # covariate
#' trt <- rbinom(n, 1, .5) # treatment
#' beta <- c(0, .8, .8, -.5) # somewhat different population
#' lp <- cbind(1, x, trt, x*trt) %*% beta # linear predictor
#' y <- rbinom(n, 1, plogis(lp)) # simulated y
#'
#' # external model
#' mod.ext <- glm(y ~ x * trt, family = "binomial")
#'
#' # predict ite's
#' control.data <- cbind(1, x, 0, 0)
#' treated.data <- cbind(1, x, 1, x)
#' y0hat <- plogis(control.data %*% coef(mod))
#' y1hat <- plogis(treated.data %*% coef(mod))
#' y0hat.updated <- plogis(control.data %*% coef(mod.ext))
#' y1hat.updated <- plogis(treated.data %*% coef(mod.ext))
#'
#' # external mbcb
#' mbcb(y0hat,y1hat,y0hat.updated,y1hat.updated)
#'
#' # small simulation
#' \dontrun{
#' B <- 50
#' dev <- ext <- rep(NA, B)
#' for(b in 1:B){
#'   n <- 250
#'   x <- rnorm(n) # covariate
#'   trt <- rbinom(n, 1, .5) # treatment
#'   beta <- c(0, 1, 1, -.75) # coefficients
#'   lp <- cbind(1, x, trt, x*trt) %*% beta # linear predictor
#'   y <- rbinom(n, 1, plogis(lp)) # simulated y
#'
#'   # model
#'   mod <- glm(y ~ x * trt, family = "binomial")
#'
#'   # predict ite's
#'   control.data <- cbind(1, x, 0, 0)
#'   treated.data <- cbind(1, x, 1, x)
#'   y0hat <- plogis(control.data %*% coef(mod))
#'   y1hat <- plogis(treated.data %*% coef(mod))
#'
#'   # mbcb
#'   dev[b] <- mbcb(y0hat,y1hat) # apparent mbcb
#'
#'   # generate external data
#'   n <- 500
#'   x <- rnorm(n) # covariate
#'   trt <- rbinom(n, 1, .5) # treatment
#'   beta <- c(0, .8, .8, -.5) # somewhat different population
#'   lp <- cbind(1, x, trt, x*trt) %*% beta # linear predictor
#'   y <- rbinom(n, 1, plogis(lp)) # simulated y
#'
#'   # external model
#'   mod.ext <- glm(y ~ x * trt, family = "binomial")
#'
#'   # predict ite's
#'   control.data <- cbind(1, x, 0, 0)
#'   treated.data <- cbind(1, x, 1, x)
#'   y0hat <- plogis(control.data %*% coef(mod))
#'   y1hat <- plogis(treated.data %*% coef(mod))
#'   y0hat.updated <- plogis(control.data %*% coef(mod.ext))
#'   y1hat.updated <- plogis(treated.data %*% coef(mod.ext))
#'
#'   # external mbcb
#'   ext[b] <- mbcb(y0hat,y1hat,y0hat.updated,y1hat.updated)
#' }
#' mean(dev)
#' mean(ext)
#' }
mbcb <- function(y0hat,y1hat,y0hat.updated=NULL,y1hat.updated=NULL){
  deltahat <- y1hat - y0hat
  if(is.null(y0hat.updated) | is.null(y1hat.updated)){
    y0hat.updated <- y0hat
    y1hat.updated <- y1hat
  }
  args <- list(
    p_0 = y0hat.updated,
    p_1 = y1hat.updated,
    tau_hat = deltahat,
    sample_weight = rep(1, length(deltahat))
  )
  mbcb <- do.call("mbcb_cpp", args)
  return(mbcb[["C Index"]])
}
