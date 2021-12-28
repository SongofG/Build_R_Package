#' Looney and Jones' Corrected Z-test
#' 
#' Performs Looney and Jones's Corrected Z-test on a partially matched sample in a form of R vector.
#' Depending on the number of missing values in the two input vectors, the test can become Two Sample Z-test or Paired Sample Z-test.
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
#' corrected.z.test(x, y)
#' corrected.z.test(x, y, alternative = "greater")
#' corrected.z.test(x, y, alternative = "less")
#' 
#' @export
corrected.z.test <- function(x, y, alternative = "two.sided"){
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
  
  T_star_bar <- mean(x[which(!is.na(x))])
  N_star_bar <- mean(y[which(!is.na(y))])
  
  ST_star <- sd(x[which(!is.na(x))])
  SN_star <- sd(y[which(!is.na(y))])
  STN1 <- cov(x[which(!is.na(x)&!is.na(y))], y[which(!is.na(x)&!is.na(y))]) # Sample Covariance of completely paired samples
  
  Z_corr <- (T_star_bar - N_star_bar)/sqrt(ST_star^2/(n1+n2) + SN_star^2/(n1+n3) - (2*n1*STN1)/((n1+n2)*(n1+n3)))
  if(alternative == "greater"){
    p <- pnorm(Z_corr, lower.tail = F)
    message <- "alternative hypothesis: true difference in means is greater than 0"
  }
  else if(alternative == "less"){
    p <- pnorm(Z_corr)
    message <- "alternative hypothesis: true difference in means is less than 0"
  }
  else{
    p <- 2*pnorm(abs(Z_corr), lower.tail = F)
    message <- "alternative hypothesis: true difference in means is not equal to 0"
  }
  
  if(n1 == 0){ # Two Sample Z-test
    cat("       ", "Two Sample Z-test\n\n", "p-value =", p, "\n", message, "\n", "number of matched:", n1, "\n", "number of partially matched pairs:", n2+n3, "\n", "Test Statistic:", Z_corr, "\n\n")
  }
  else if(n2==0 & n3==0){ # Paired sample Z-test
    cat("       ", "Paired sample Z-test\n\n", "p-value =", p, "\n", message, "\n", "number of matched:", n1, "\n", "number of partially matched pairs:", n2+n3, "\n", "Test Statistic:", Z_corr, "\n\n")
  }
  else{ # The same Corrected Z-test but warning message if the condition is met
    cat("       ", "Looney and Zones' Corrected Z-test\n\n", "p-value =", p, "\n", message, "\n", "number of matched:", n1, "\n", "number of partially matched pairs:", n2+n3, "\n", "Test Statistic:", Z_corr, "\n\n")
    if((n2==0&n3!=0) | (n2!=0&n3==0)){
      warning("There is only one input vector with missing values")
    }
  }
}
