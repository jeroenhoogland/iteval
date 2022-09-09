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
  n <- length(deltahat)
  snum <- sden <- rep(NA, n)
  g <- rbind(c(0,1,0,0),
             c(0,0,1,0),
             c(0,1,1,0),
             c(1,1,1,0),
             c(0,1,1,1))
  j=1

  if(is.null(y0hat.updated) | is.null(y1hat.updated)){
    y0hat.updated <- y0hat
    y1hat.updated <- y1hat
  }

  for(j in 1:n){
    jprime <- (1:n)[-j]

    deltaless <- ifelse(deltahat[j] < deltahat[jprime], 1, 0)
    deltaeq <- ifelse(deltahat[j] == deltahat[jprime], 1, 0)

    ggj <- ifelse(g[ ,1], y1hat.updated[j], 1-y1hat.updated[j]) * ifelse(g[ ,2], y0hat.updated[j], 1-y0hat.updated[j])
    ggjprime <- (outer(g[ ,3],  y1hat.updated[-j], FUN = function(X,Y) X*Y) + outer(1-g[ ,3],  1-y1hat.updated[-j], FUN = function(X,Y) X*Y)) *
      (outer(g[ ,4],  y0hat.updated[-j], FUN = function(X,Y) X*Y) + outer(1-g[ ,4],  1-y0hat.updated[-j], FUN = function(X,Y) X*Y))
    oless <- apply(ggjprime, 2, function(x) sum(x*ggj))

    snum[j] <- sum(deltaless * oless) + 1/2 * sum(deltaeq * oless)
    sden[j] <- sum(oless)
  }
  sum(snum) / sum(sden) # mbcben
}
