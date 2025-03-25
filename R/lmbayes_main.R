
##' @description Conjugate update in Bayesian linear model
##' @title Conjugate update in Bayesian linear model
##' @param formula. A formula
##' @param data. A dataframe
##' @param prior named list m (Mean in prior distribution of beta), C (Variance in prior distribution of beta), sd (Conditional variance of y given beta is sd I)
##' @param update.sd Should conditional variance be updated (emperical Bayes)
##' @return An object of class 'blm'
##' @author Søren Højsgaard
##'
##' @examples
##' form <- dist ~ speed + I(speed^2)
##' mm <- lm(form, data=cars)
##' coef(summary(mm))
##' sigma(mm)
##'
##' bb <- lmb(form, data=cars, prior=list(m=c(1, 1, 1), C=diag(10,3), sd=2))
##' bb
##'
##' coef(bb)
##' sigma(bb)
##'
##' if (require(ggplot2)){
##' pp <- qplot(speed, dist, data=cars)
##' pp + geom_line(aes(x=speed, y=predict(bb)), col="red") +
##'      geom_line(aes(x=speed, y=predict(mm)), col="blue")
##' }

#' @export
lm2prior <- function(lm1, scale=1){
    prior <- list(m=coef(lm1), C=diag(diag(vcov(lm1))) * scale, sd=sigma(lm1)*scale)
    prior
}


#' Linear Model with Empirical Bayes Shrinkage
#'
#' Fits a linear regression model with optional shrinkage toward a prior mean, using empirical Bayes logic.
#' If no prior is provided, the function defaults to standard OLS via \code{lm()}.
#'
#' @param formula. A formula specifying the linear model (e.g., \code{y ~ x + I(x^2)}).
#' @param data. A data frame containing the variables in the model.
#' @param prior A list containing prior information with components:
#'   \describe{
#'     \item{\code{m}}{Numeric vector of prior means for the regression coefficients.}
#'     \item{\code{C}}{Prior covariance matrix (not currently used in simple shrinkage).}
#'     \item{\code{sd}}{(Optional) Prior residual standard deviation (can be used for scaling/shrinkage).}
#'   }
#'   If \code{prior} is \code{NULL}, the function fits a standard OLS model.
#' @param alpha A numeric value in \code{[0, 1]} controlling shrinkage:
#'   \describe{
#'     \item{\code{alpha = 0}}{Full reliance on data (pure OLS).}
#'     \item{\code{alpha = 1}}{Full reliance on prior (ignore data).}
#'     \item{\code{0 < alpha < 1}}{Convex combination of prior mean and OLS estimate.}
#'   }
#' @param update.sd Logical. If \code{TRUE}, estimate the residual standard deviation from the data.
#'
#' @return An object of class \code{\"lmb_class\"}, which is a list with the following components:
#' \describe{
#'   \item{\code{call}}{The matched function call.}
#'   \item{\code{formula}}{The model formula used.}
#'   \item{\code{prior}}{The prior provided (or \code{NULL} if none).}
#'   \item{\code{posterior}}{A list with posterior components:}
#'     \describe{
#'       \item{\code{m}}{Posterior mean vector of regression coefficients.}
#'       \item{\code{C}}{Posterior covariance matrix (may be placeholder if not computed).}
#'       \item{\code{sd}}{Posterior residual standard deviation (from OLS if \code{update.sd = TRUE}).}
#'     }
#'   \item{\code{X}}{The design matrix used in the model.}
#'   \item{\code{data}}{The data used to fit the model.}
#' }
#'
#' @examples
#' ## Example: Bayesian shrinkage using the built-in 'cars' dataset
#' data(cars)
#' 
#' # A pretend prior: intercept = 0, slope = 3, with some vague prior sd
#' prior <- list(
#'  m = c(0, 3),              # prior mean: intercept = 0, slope = 3
#'  C = diag(c(100, 10)),     # vague prior covariance (not used in this shrinkage)
#'  sd = 15                   # prior residual sd (not used here)
#')
#'
#' # OLS
#' fit0 <- lmb(dist ~ speed, data = cars)
#'
#' # Bayesian
#' fit1 <- lmb(dist ~ speed, data = cars, prior = prior)
#'
#' # Shrinkage combination of prior and OLS
#' fit2 <- lmb(dist ~ speed, data = cars, prior = prior, alpha = 0.4)
#'
#' # Posterior coefficients
#' fit0$posterior ## Nothing
#' fit1$posterior
#' fit2$posterior

#' @export
lmb <- function(formula., data., prior=NULL, alpha=0, update.sd=TRUE){

    cl  <- match.call()
    mf  <- model.frame(formula., data=data.)
    X   <- model.matrix(mf, data=data.)
    nms <- colnames(X) 
    ols <- lm(formula., data=data.)
    sd.ols   <- sigma(ols)
    beta.ols <- coef(ols)
    
    if (is.null(prior)) {
        ## message("No prior supplied — fitting OLS model via lm().")
        return(ols)
    }
    
    if (alpha > 0) {
        m  <- prior$m
        m2 <- alpha * m + (1-alpha) * beta.ols
        posterior <- list(m=m2)        
    } else {            
        m  <- prior$m
        C  <- prior$C 
        sd <- prior$sd 

        yy <- response(ols)
        if (any(is.na(yy)))
            stop("missing values in response\n")
        
        if (!is.matrix(C)) {
            if (length(C) == 1) {
                C <- diag(C, length(m))            
            }
        }
        
        if (update.sd){
            R <- diag(sd.ols^2, length(yy))
        } else {
            R <- diag(sd^2, length(yy))
        }
        
        XCXt <- X %*% C %*% t(X)
        Cyy <- XCXt + R
        Cyb <- X %*% C 
        Cby <- t(Cyb)    
        Ey  <- as.numeric(X %*% m)    
        m2  <- m + Cby %*% solve(Cyy, (yy - Ey)) ## posterior
        C2  <- C - Cby %*% solve(Cyy, Cyb)       ## posterior
        m2  <- as.numeric(m2)
        names(m2) <- nms
        
        if (update.sd){
            X  <- model.matrix(ols)
            pp <- as.numeric(X %*% m2)
            rr <- yy - pp
            sd2 <- sd(rr)
        } else {
            sd2 <- sd
        }
        posterior <- list(m=m2,  C=C2,  sd=sd2)        
    }
    
    out <- list(call=cl, formula=formula.,
                prior     = prior,
                posterior = posterior,
                X=X, data=data.)
    
    class(out) <- "lmb_class"
    out
}

#' Convert lmb_class object to lm object
#' @export
#' @param object A lmb_class object
as_lm <- function(object){
    if (!inherits(object, "lmb_class"))
        stop("'object' must be lmb_class\n")

    out <- lm(object$formula, object$data)
    return(out)
}



    ## if (is.numeric(sd) && length(sd)==1)
        
    ## else
    ##     stop("sd is not a number")

#' 
#' @importFrom stats "coef" "fitted" "formula" "lm" "model.frame" sigma simulate predict
#'               "model.matrix" "quantile" "resid" "update" "vcov"


#' @export
coef.lmb_class <- function(object, ...){
    object$posterior$m
}

#' @export
vcov.lmb_class <- function(object, ...){
    object$posterior$C
}

#' @export
sigma.lmb_class <- function(object, ...){
    object$posterior$sd
}

#' @export
model.matrix.lmb_class <- function(object, ...){
    object$X
}

#' @export
print.lmb_class <- function(x, ...){
    cat("Bayesian linear model object\n")
    print(formula(x))
    print(as.numeric(coef(x)))
    invisible(x)
}

#' @importFrom broom tidy
#' @export
tidy.lmb_class <- function(x, ...){
    data.frame(estimate=coef(x),std.error=diag(vcov(x)))
}

get_model_matrix <- function(object, newdata=NULL){
    if (is.null(newdata))
        model.matrix(object)
    else model.matrix(object$formula, data=newdata)    
}

#' @export
predict.lmb_class <- function(object, newdata=NULL, ...){
    as.numeric(get_model_matrix(object, newdata) %*% coef(object))
}

#' @export
residuals.lmb_class <- function(object, ...){
  yy <- object$data[,deparse(object$formula[[2]])]
  pp <- predict(object)
  rr <- yy - pp
  rr
}



#' @importFrom MASS mvrnorm
#' @export
simulate.lmb_class <- function(object, nsim=1, seed=NULL, newdata=NULL, ...){
    if (nsim <=1){
        stop("Must have nsim > 1\n")
    }
    bb <- MASS::mvrnorm(n=nsim, mu=coef(object), Sigma=vcov(object))
    bb <- t(bb)
    mm <- get_model_matrix(object, newdata)
    pp <- mm %*% bb
    er <- sigma(object)^2 * diag(1, nrow(mm))
    ee <- MASS::mvrnorm(n=nsim, mu=rep(0, nrow(mm)), Sigma=er)
    ee <- t(ee)
    ## print(ee)    
    pp + ee
}

#' Prediction interval
#' @param object A 'lmb_class' object
#' @param level Level for interval
#' @param nsim Number of simulations to base the interval on
#' @param newdata New dataframe
#' @param ... Additional arguments, currently unused.
#' @author Søren Højsgaard

#' @export
predint <- function(object, level=0.95, nsim=100, newdata=NULL, ...){

    ss <- simulate(object, nsim=nsim, newdata=newdata)

    ## ee <- MASS::mvrnorm(n=nsim, mu=rep(0,nrow(sigma(object))), Sigma=sigma(object))
    ## print(dim(ss)); print(dim(ee))
     
    pp <- c((1-level) / 2, 1-(1-level)/2)
    out <- t(sapply(split(ss, row(ss)), function(r) quantile(r, pp)))
    out <- data.frame(cbind(predict(object, newdata), out))
    names(out) <- c("pred", "lwr", "upr")
    out
}
