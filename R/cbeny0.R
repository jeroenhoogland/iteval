#' Matching-based (on y0hat) concordance for benefit
#'
#' @param y0hat Predicted outcome under control treatment.
#' @param y1hat Predicted outcome under treatment.
#' @param y Observed outcome of interest.
#' @param ind.A Index for individuals in the control group (should be a vector
#' of length y with the indices that belong to control individuals).
#' @param ind.B Index for individuals in the treated group (should be a vector
#' of length y with the indices that belong to treated individuals).
#' @param match Currently set to "y0hat" (i.e., match on predicted outcome under
#' control treatment).
#' @param nresample In case of unequal  numbers in the treatment groups, the
#' concordance estimate is based on repeated matching of the smallest group to
#' equal sized subsamples of the larger group. nresample denotes the number of
#' replications of this procedure.
#' @param get.all get the c-for-ben for all resamplings.
#'
#' @return The concordance estimate(s) for predicted ITE
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
#' # cbeny0
#' table(trt) # check treatment group sizes
#' ind.A <- which(trt == 0)
#' ind.B <- which(trt == 1)
#' res <- cbeny0(y0hat, y1hat, y, ind.A, ind.B, "y0hat",
#' nresample = 250, get.all=TRUE)
#' mean(res) # cbeny0
#' hist(res); abline(v=mean(res), col="red")

cbeny0 <- function(y0hat, y1hat, y,
                     ind.A, ind.B,
                     match="y0hat",
                     nresample=1,
                     get.all=FALSE){

  deltahat <- y1hat - y0hat

  order.A <- order(y0hat[ind.A])
  ind.A.ord <- ind.A[order.A]
  order.B <- order(y0hat[ind.B])
  ind.B.ord <- ind.B[order.B]

  if(length(ind.A.ord) == length(ind.B.ord)){
    delta.ij <- y1hat[ind.B.ord] - y0hat[ind.A.ord]
    obs.ben.ij <- y[ind.B.ord] - y[ind.A.ord]
    return(Hmisc::rcorr.cens(delta.ij, obs.ben.ij)["C Index"])
  } else {
    cstat <- rep(NA, nresample)
    for(r in 1:nresample){
      if(length(ind.A.ord) < length(ind.B.ord)){
        omit <- abs(length(ind.A.ord) - length(ind.B.ord))
        omit <- sample(ind.B.ord, omit, replace = FALSE)
        ind.B.ord.new <- ind.B.ord[-which(ind.B.ord %in% omit)]
        delta.ij <- y1hat[ind.B.ord.new] - y0hat[ind.A.ord]
        obs.ben.ij <- y[ind.B.ord.new] - y[ind.A.ord]
        cstat[r] <- Hmisc::rcorr.cens(delta.ij, obs.ben.ij)["C Index"]
      }
      if(length(ind.A.ord) > length(ind.B.ord)){
        omit <- abs(length(ind.A.ord) - length(ind.B.ord))
        omit <- sample(ind.A.ord, omit, replace = FALSE)
        ind.A.ord.new <- ind.A.ord[-which(ind.A.ord %in% omit)]
        delta.ij <- y1hat[ind.B.ord] - y0hat[ind.A.ord.new]
        obs.ben.ij <- y[ind.B.ord] - y[ind.A.ord.new]
        cstat[r] <- Hmisc::rcorr.cens(delta.ij, obs.ben.ij)["C Index"]
      }
    }
    if(!get.all){
      return(c("C Index"=mean(cstat)))
    } else {
      return(c("C Index"=cstat))
    }
  }
}
