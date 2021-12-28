#' Kim et al.'s Modified T-test
#' 
#' Performs Kim et al.'s Modified T-test on a partially matched sample in a form of R vector.
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
#' modified.t.test(x, y)
#' modified.t.test(x, y, alternative = "greater")
#' modified.t.test(x, y, alternative = "less")
#' 
#' @export
modified.t.test <- function(x, y, alternative="two.sided"){
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
  
  equal = F
  if(var.test(x,y)$p.value > 0.05){
    equal = T
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
  
  D_bar <- mean(x[which(!is.na(x)&!is.na(y))]-y[which(!is.na(x)&!is.na(y))])
  T_bar <- mean(x[which(!is.na(x)&is.na(y))])
  N_bar <- mean(y[which(!is.na(y)&is.na(x))])
  
  SD <- sd(x[which(!is.na(x)&!is.na(y))]-y[which(!is.na(x)&!is.na(y))])
  ST <- sd(x[which(!is.na(x)&is.na(y))])
  SN <- sd(y[which(!is.na(y)&is.na(x))])
  
  nH <- 2/(1/n2 + 1/n3)
  
  t3 <- (n1*D_bar + nH*(T_bar - N_bar))/sqrt(n1*SD^2 + nH^2*(SN^2/n3 + ST^2/n2))
  
  p <- NULL
  
  if(alternative == "greater"){
    p <- pnorm(t3, lower.tail = F)
    message <- "alternative hypothesis: true difference in means is greater than 0"
  }
  else if(alternative == "less"){
    p <- pnorm(t3)
    message <- "alternative hypothesis: true difference in means is less than 0"
  }
  else{
    p <- 2*pnorm(abs(t3), lower.tail = F)
    message <- "alternative hypothesis: true difference in means is not equal to 0"
  }
  
  cat("       ", "Kim et al.'s Modified T-test\n\n", "p-value =", p, "\n", message, "\n", "number of matched:", n1, "\n", "number of partially matched pairs:", n2+n3, "\n", "Test Statistic:", t3, "\n\n")
  
}
