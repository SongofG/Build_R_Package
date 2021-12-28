#' Lin and Stiver's MLE Based Test Under Heteroscedasticity
#' 
#' Performs Lin and Stiver's MLE Based Test Under Heteroscedasticity on a partially matched sample in a form of R vector.
#' This test is under the assumption of unequal variance of the two input vectors. By F-test, if the two vectors are considered to have the same variance, a warning message will be shown after the result is printed.
#' 
#' @param x a (non-empty) numeric vector of data values with some missing value(NA).
#' @param y a (non-empty) numeric vector of data values with some missing value(NA).
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' 
#' @return Results of test including p-value will be printed.
#' 
#' @examples 
#' # Generating Toy Examples
#' set.seed(123)
#' x <- rnorm(20)
#' x[sample(1:20, 3)] <- NA # Deliverately generating some missing values
#' y <- (rnorm(20) + 1)/3
#' y[sample(which(!is.na(x)), 4)] <- NA
#' lin.stiver.test(x, y)
#' lin.stiver.test(x, y, alternative = "greater")
#' lin.stiver.test(x, y, alternative = "less")
#' 
#' @export

lin.stiver.test <- function(x, y, alternative="two.sided"){
  if(is.null(x) | is.null(y)){
    stop("Both of the input vectors should not be NULL")
  }
  else if(length(x) != length(y)){
    stop("Two input vectors should have the same length")
  }
  else if(is.double(x) == F | is.double(y) == F){
    stop("Only double vectors are accepted")
  }
  else if(length(x)==length(which(is.na(x))) | length(y)==length(which(is.na(y)))){
    stop("The vector is consisted of all NAs.")
  }
  else if(length(x)==1 | length(y)==1){
    stop("Test cannot be conducted due to the small size of sample")
  }
  else if(alternative != "two.sided" & alternative != "greater" & alternative != "less"){
    stop("Argument 'alternative' should be one of 'two.sided', 'greater', or 'less'.")
  }
  
  n1 <- length(which(!is.na(x)&!is.na(y))) # The number of complete matched pairs
  n2 <- length(which(!is.na(x)&is.na(y))) # The number of samples that are not NA in x but NA in y
  n3 <- length(which(is.na(x)&!is.na(y))) # The number of samples that are NA in x but not in y
  
  # n1, n2, n3 case control.
  if(n1 == 0){ # No paired data
    warning("Since there are no paired data, the test will automatically become two sample t-test.")
    return(t.test(x,y,alternative = alternative, var.equal = equal))
  }
  else if(n2==0 & n3==0){ # all paired data
    warning("Since there is no missing value from all the paired data in the sample, the test will automatically become paired two sample t-test.")
    return(t.test(x,y,alternative = alternative, paired = T, var.equal = equal))
  }
  else if((n2==0 & n3!=0) | (n2!=0 & n3==0)){
    warning("Since there are missing values from only one input vector, the test will automatically become two sample t-test")
    return(t.test(x,y,alternative = alternative, var.equal = equal))
  }
  
  STN1 <- cov(x[which(!is.na(x)&!is.na(y))], y[which(!is.na(x)&!is.na(y))])
  ST1 <- sd(x[which(!is.na(x)&!is.na(y))])
  SN1 <- sd(y[which(!is.na(x)&!is.na(y))])
  r <- STN1/(ST1*SN1)
  f <- (n1*(n1+n3+n2*(STN1/ST1^2)))/((n1+n2)*(n1+n3) - n2*n3*r^2)
  g <- (n1*(n1+n2+n3*(STN1/SN1^2)))/((n1+n2)*(n1+n3) - n2*n3*r^2)
  V_1 <- ((f^2/n1 + (1-f)^2/n2)*ST1^2*(n1-1) + (g^2/n1 + (1-g)^2/n3)*SN1^2*(n1-1) - 2*f*g*STN1*(n1-1)/n1)/(n1-1)
  
  Z_LS <- (f*(mean(x[which(!is.na(x) & !is.na(y))]) - mean(x[which(!is.na(x)&is.na(y))]))
           - g*(mean(y[which(!is.na(x) & !is.na(y))]) - mean(y[which(!is.na(y)&is.na(x))]))
           + mean(x[which(!is.na(x)&is.na(y))]) - mean(y[which(!is.na(y)&is.na(x))]))/sqrt(V_1)
  
  if(alternative == "greater"){
    p <- pt(Z_LS, n1, lower.tail = F)
    message <- "alternative hypothesis: true difference in means is greater than 0"
  }
  else if(alternative == "less"){
    p <- pt(Z_LS, n1)
    message <- "alternative hypothesis: true difference in means is less than 0"
  }
  else{
    p <- 2*pt(abs(Z_LS), n1, lower.tail = F)
    message <- "alternative hypothesis: true difference in means is not equal to 0"
  }
  
  cat("       ", "Lin and Stiver's MLE Based Test Under Heteroscedasticity\n\n", "p-value =", p, "\n", message, "\n", "number of matched:", n1, "\n", "number of partially matched pairs:", n2+n3, "\n", "Test Statistic:", Z_LS, "\n\n")
  if(var.test(x,y)$p.value <= 0.05){
    warning("F-test of the two vectors present a p-value less than 0.05. This test is under the equal variance assumption.")
  }
}

