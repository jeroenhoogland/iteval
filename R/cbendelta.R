#' Title
#'
#' This function derives the c-for-benefit based on 1:1 control:treated matches,
#' with matching based on predicted individualized treatment effect. The main
#' reference is <DOI:10.1016/j.jclinepi.2017.10.021>.
#'
#' @param deltahat Predicted individualized treatment effects
#' @param y Observed outcome of interest
#' @param ind.A Index for individuals in the control group (should be a vector
#' of length y with the indices that belong to control individuals).
#' @param ind.B Index for individuals in the treated group (should be a vector
#' of length y with the indices that belong to treated individuals).
#' @param nresample In case of unequal  numbers in the treatment groups, the
#' concordance estimate is based on repeated matching of the smallest group to
#' equal sized subsamples of the larger group. nresample denotes the number of
#' replications of this procedure.
#' @param get.all get the c-for-ben for all resamplings.
#'
#' @return The concordance estimate(s)
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
#' deltahat <- plogis(treated.data %*% coef(mod)) -
#' plogis(control.data %*% coef(mod))
#'
#' # cbendelta
#' table(trt) # check treatment group sizes
#' ind.A <- which(trt == 0)
#' ind.B <- which(trt == 1)
#' res <- cbendelta(deltahat, y, ind.A, ind.B, nresample = 100, get.all=TRUE)
#' mean(res) # cbendelta
#' hist(res); abline(v=mean(res), col="red")

cbendelta <- function(deltahat, y, ind.A, ind.B, nresample=1, get.all=FALSE){
  order.A <- order(deltahat[ind.A])
  ind.A.ord <- ind.A[order.A]
  order.B <- order(deltahat[ind.B])
  ind.B.ord <- ind.B[order.B]
  if(length(ind.A.ord) == length(ind.B.ord)){
    delta.avg <- (deltahat[ind.A.ord] + deltahat[ind.B.ord]) / 2
    obs.ben <- y[ind.B.ord] - y[ind.A.ord]
    return(Hmisc::rcorr.cens(delta.avg, obs.ben)["C Index"])
  } else {
    cstat <- rep(NA, nresample)
    for(r in 1:nresample){
      if(length(ind.A.ord) < length(ind.B.ord)){
        omit <- abs(length(ind.A.ord) - length(ind.B.ord))
        omit <- sample(ind.B.ord, omit, replace = FALSE)
        ind.B.ord.new <- ind.B.ord[-which(ind.B.ord %in% omit)]
        delta.avg <- (deltahat[ind.A.ord] + deltahat[ind.B.ord.new]) / 2
        obs.ben <- y[ind.B.ord.new] - y[ind.A.ord]
        cstat[r] <- Hmisc::rcorr.cens(delta.avg, obs.ben)["C Index"]
      }
      if(length(ind.A.ord) > length(ind.B.ord)){
        omit <- abs(length(ind.A.ord) - length(ind.B.ord))
        omit <- sample(ind.A.ord, omit, replace = FALSE)
        ind.A.ord.new <- ind.A.ord[-which(ind.A.ord %in% omit)]
        delta.avg <- (deltahat[ind.A.ord.new] + deltahat[ind.B.ord]) / 2
        obs.ben <- y[ind.B.ord] - y[ind.A.ord.new]
        cstat[r] <- Hmisc::rcorr.cens(delta.avg, obs.ben)["C Index"]
      }
    }
    if(!get.all){
      return(c("C Index"=mean(cstat)))
    } else {
      return(c("C Index"=cstat))
    }
  }
}
