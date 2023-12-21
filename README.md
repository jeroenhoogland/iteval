
<!-- README.md is generated from README.Rmd -->

# iteval

<!-- badges: start -->
<!-- badges: end -->

The goal of iteval is to ease implementation of the discrimination and
calibration measures for predicted individualized treatment effect as
described by Hoogland et al. \[yet to insert DOI\].

## Installation

You can install the development version of **iteval** from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jeroenhoogland/iteval")
```

and load the package using

``` r
library(iteval)
```

## Example

This is a basic example which illustrates usage of the main functions.
For this purpose, we first simulate some data from a data generating
mechanism with a single covariate `x`, a treatment indicator `trt`, an
interaction between the two and observed binary outcome `y` and fit a
basic logistic regression model.

``` r
# generate some  data
set.seed(123) # for replicability
n <- 500
x <- rnorm(n) # covariate
trt <- rbinom(n, 1, .5) # treatment
beta <- c(0, 1, 1, -.75) # coefficients
lp <- cbind(1, x, trt, x*trt) %*% beta # linear predictor
y <- rbinom(n, 1, plogis(lp)) # simulated y

# model
mod <- glm(y ~ x * trt, family = "binomial")
```

Next, we obtain predicted of y under control treatment (`y0hat`),
treated condition (`y1hat`) and predicted ite’s (`deltahat`).

``` r
# Predicited ite's
control.data <- cbind(1, x, 0, 0) 
treated.data <- cbind(1, x, 1, x)

# prediction under the control condition
y0hat <- plogis(control.data %*% coef(mod))

# prediction under the treated condition
y1hat <- plogis(treated.data %*% coef(mod))

# and their difference
deltahat <- y1hat - y0hat

# true outcome probabilities given data generating mechanism
p0 <- plogis(control.data %*% c(0, 1, 1, -.75))
p1 <- plogis(treated.data %*% c(0, 1, 1, -.75))
```

Subsequently, we apply the performance measures in the simulated data.
NB.: note that these estimates will be optimistic since this same data
set was used to fit the prediction model; in practice, indepdenet data
or an internal validation procedure is required.

### Measures of discriminative performance

Each of the measures need to be aware of the treatment assignments.
`ind.A` contains the indices for the controls and `ind.B` contains the
indices for the treated.

``` r
ind.A <- which(trt == 0)
ind.B <- which(trt == 1)
```

All is in place now. The treatment assignments are as follows:

``` r
table(trt)
#> trt
#>   0   1 
#> 251 249
```

Hence, 1:1 matching excluded some individuals and repeated subsampling
of the larger group is used more of the data (`nresample=100`).

``` r
# cbendelta
res <- cbendelta(deltahat, y, ind.A, ind.B, nresample = 250, get.all=TRUE)
mean(res) # cbendelta
#> [1] 0.6359874
hist(res, xlab="cbendelta", main=""); abline(v=mean(res), col="red")
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

``` r

# cbeny0
res <- cbeny0(y0hat, y1hat, y, ind.A, ind.B, "y0hat", nresample = 250, get.all=TRUE)
mean(res) # cbendelta
#> [1] 0.6424884
hist(res, xlab="cbeny0", main=""); abline(v=mean(res), col="red")
```

<img src="man/figures/README-unnamed-chunk-6-2.png" width="100%" />

``` r

# mbcb (does not depend on matching)
mbcb(y0hat, y1hat)
#> [1] 0.6318088

# 'true' mbcb for this model (possible since sim setting)
mbcb(y0hat, y1hat, y0hat.updated = p0, y1hat.updated = p1)
#> [1] 0.6266713
```

### Calibration performance

``` r
# Example 'true' calibration (possible since sim setting)
# otherwise, cal() takes estimates p0 and p1 based on 
# independent data.
(cal <- cal(y0hat,y1hat,y,trt,p0,p1,type="both"))
#> 
#> Call:  stats::glm(formula = p01 ~ 0 + Z + deltahatlpstar, family = "quasibinomial", 
#>     offset = os)
#> 
#> Coefficients:
#>              Z  deltahatlpstar  
#>        0.09199         0.90527  
#> 
#> Degrees of Freedom: 500 Total (i.e. Null);  498 Residual
#> Null Deviance:       166.4 
#> Residual Deviance: 1.083e-14     AIC: NA
```
