
##' @description Conjugate update in Bayesian linear model
##' @title Conjugate update in Bayesian linear model
##' @param formula A formula
##' @param data A dataframe
##' @param m Mean in prior distribution of beta 
##' @param C Variance in prior distribution of beta
##' @param phi Conditional variance of y given beta is phi I.
##' @return An object of class 'blm'
##' @author Søren Højsgaard
##'
##' @examples
##' form <- dist ~ speed + I(speed^2)
##' mm <- lm(form, data=cars)
##' coef(summary(mm))
##' sigma(mm)
##'
##' bb <- bayes_lm(form, data=cars, m=c(1, 1, 1), C=diag(10,3), phi=2)
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
bayes_lm <- function(formula, data, m, C, phi){
    mf <- model.frame(formula, data=data)
    X = model.matrix(mf, data=data)
    
    dummy <- lm(update(formula, ~1), data=data)
    yy <- unname(fitted(dummy) + resid(dummy))

    if (any(is.na(yy))) stop("missing values in response\n")

    if (!is.matrix(C))
        if (length(C) == 1) C <- diag(C, length(m))

    if (is.numeric(phi) && length(phi)==1)
        W <- diag(phi, length(yy))
    else
        stop("phi is not a number")
        
    Cyy = X %*% C %*% t(X) + W
    Cyb <- X %*% C 
    Cby <- t(Cyb)
    
    Ey <- as.numeric(X %*% m)
    
    m2 <- as.numeric(m + Cby %*% solve(Cyy, (yy - Ey)))
    C2 <- C - Cby %*% solve(Cyy, Cyb)
    
    out <- list(formula=formula, m=m2, C=C2, phi=phi, W=W, X=X)
    class(out) <- "blm"
    out
}

#' @importFrom stats "coef" "fitted" "formula" "lm" "model.frame" sigma simulate predict
#'               "model.matrix" "quantile" "resid" "update" "vcov"


#' @export
coef.blm <- function(object, ...){
    object$m
}

#' @export
vcov.blm <- function(object, ...){
    object$C
}

#' @export
sigma.blm <- function(object, ...){
    object$phi
}

#' @export
model.matrix.blm <- function(object, ...){
    object$X
}

#' @export
print.blm <- function(x, ...){
    cat("Bayesian linear model object\n")
    print(formula(x))
    print(coef(x))
    invisible(x)
}


get_model_matrix <- function(object, newdata){
    if (missing(newdata)) model.matrix(object)
    else model.matrix(object$formula, data=newdata)    
}

#' @export
predict.blm <- function(object, newdata, ...){
    as.numeric(get_model_matrix(object, newdata) %*% coef(object))
}

#' @importFrom MASS mvrnorm
#' @export
simulate.blm <- function(object, nsim=1, seed=NULL, newdata, ...){
    bb <- MASS::mvrnorm(n=nsim, mu=coef(object), Sigma=vcov(object))
    if (!is.matrix(bb)) bb <- matrix(bb, nrow=1)
    get_model_matrix(object, newdata) %*% t(bb)    
}

#' Prediction interval
#' @param object A 'blm' object
#' @param level Level for interval
#' @param nsim Number of simulations to base the interval on
#' @param newdata New dataframe
#' @param ... Additional arguments, currently unused.
#' @author Søren Højsgaard

#' @export
predint <- function(object, level=0.95, nsim=100, newdata, ...){

    ss <- simulate(object, nsim=100, newdata=newdata)

    ## ee <- MASS::mvrnorm(n=nsim, mu=rep(0,nrow(sigma(object))), Sigma=sigma(object))
    ## print(dim(ss)); print(dim(ee))
     
    pp <- c((1-level) / 2, 1-(1-level)/2)
    out <- t(sapply(split(ss, row(ss)), function(r) quantile(r, pp)))
    out <- data.frame(cbind(predict(object, newdata), out))
    names(out) <- c("pred", "lwr", "upr")
    out
}


