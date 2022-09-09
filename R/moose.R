#' moose: mean squared out-of-sample error projection
#'
#' This function projects the mean squared out-of-sample error for a linear regression
#' 
#' @param reg an lm object containing the regression to project out-of-sample
#' @param dataset a data.frame containing new cases for out-of-sample projection
#' @keywords generalization
#' @export
#' @examples
#' # set the seed for reproducibility of the example
#' set.seed(04251978)
#' # randomly generate 100 observations of data
#' mydata <- data.frame(x1=rnorm(100),x2=rnorm(100),x3=rnorm(100))
#' # true outcome variable is y = x1 + x2 + x3 + e
#' y <- mydata$x1 + mydata$x2 + mydata$x3 + rnorm(100)
#' # regression with the first 25 observations from the dataset
#' reg <- lm(y ~ x1 + x2 + x3,data=cbind(y,mydata)[1:25,])
#' # using the predictor values from the first 25 observations,
#' # project the out-of-sample error we can expect in the case of
#' # "non-stochastic" predictors whose values are the same in the
#' # test sample as in the training sample.
#' # note that mydata does not include the outcome variable.
#' same.predictor.values.error <- moose(reg,mydata[1:25,])
#' # by comparison, the in-sample R-squared value observed
#' # in training is:
#' summary(reg)$r.squared
#' # using the predictor values from the next 75 obsevervations,
#' # project the out-of-sample error we can expect in the case
#' # of stochastic predictors whose values potentially differ
#' # from those used in training.
#' new.predictor.values.error <- moose(reg,mydata[26:100,])
#' # by comparison, the actual mse and out-of-sample R-squared value
#' # obtained from observations 26-100 of this random sample are:
#' mse <- mean((y[26:100]-predict(reg,mydata[26:100,]))^2)
#' mse
#' m.total.sqs <- mean((y[26:100]-mean(y[26:100]))^2)
#' r2o <- 1-mse/m.total.sqs
#' r2o
moose <- function(reg,dataset){

  if(!("lm" %in% class(reg))){
    stop("First input must be an lm object.\n")
  }
  if(!("data.frame" %in% class(dataset))){
    stop("Second input must be a data.frame.\n")
  }
  if("weights" %in% names(reg)){
    stop("Weights not supported.\n")
  }

  x.locs <- match(names(reg$model)[2:ncol(reg$model)],names(dataset))

  if(any(is.na(x.locs))){
    missing.vars <- paste(names(reg$model)[1+which(is.na(x.locs))],collapse=", ")
    stop(paste("Predictors not found in out-of-sample data:",missing.vars,"\n"))
  }

  xmatrix <- as.matrix(reg$model)
  Xo <- data.frame(dataset)[,c(x.locs)]
  Xo <- Xo[complete.cases(Xo),]

  if(nrow(Xo)==0){
    stop("No complete cases found in out-of-sample data.\n")
  }

  if(names(reg$coefficients)[1]=="(Intercept)"){
    xmatrix[,1] <- 1
    Xo <- cbind(rep(1,nrow(Xo)),Xo)
  } else {
    xmatrix <- xmatrix[,2:dim(xmatrix)[2]]
  }

	xpxinv.xp <- chol2inv(qr.R(reg$qr)) %*% t(xmatrix)
  Ho <- as.matrix(Xo) %*% xpxinv.xp
  Ho2 <- Ho^2
  hat <- diag(xmatrix %*% xpxinv.xp)
  e <- as.numeric(reg$residuals)
  sigma2 <- as.matrix(e^2/(1-hat))
  mse <- mean(as.numeric((Ho+Ho2) %*% sigma2))
  hato <- as.numeric(rowSums(Ho2))

  insample.squares <- as.matrix((reg$model[,1]-mean(reg$model[,1]))^2)
  tss <- mean(as.numeric(Ho %*% insample.squares))
  R2o <- 1-mse/tss

  output <- list(mse=mse,R2o=R2o,hat=hato)

  cat(paste("Complete out-of-sample cases:     ",nrow(Xo),"\n"))
  cat(paste("Projected mean squared oos error: ",mse,"\n"))
  cat(paste("Projected out-of-sample R-squared:",R2o,"\n"))

  return(output)

}

