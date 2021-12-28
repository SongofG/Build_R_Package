#' Liptak's Weighted Z-test
#' 
#' Performs Liptak's Weighted Z-test on a partially matched sample in a form of R vector.
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
#' weighted.z.test(x, y)
#' weighted.z.test(x, y, alternative = "greater")
#' weighted.z.test(x, y, alternative = "less")
#' @export
weighted.z.test <- function(x, y, alternative="two.sided"){
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
  
  w1 <- sqrt(2*n1) # weight of the complete matched pairs
  w2 <- sqrt(n2+n3) # weight of the partially matched pairs
  
  if (alternative == "greater"){
    p1 <- t.test(x[which(!is.na(x)&!is.na(y))], y[which(!is.na(x)&!is.na(y))], alternative = "greater", paired = T, var.equal=equal)$p.value
    p2 <- t.test(x[which(!is.na(x)&is.na(y))], y[which(is.na(x)&!is.na(y))], alternative = "greater", var.equal=equal)$p.value
    Z1 <- qnorm(1-p1)
    Z2 <- qnorm(1-p2)
    p <- pnorm((w1*Z1 + w2*Z2)/sqrt(w1^2 + w2^2), lower.tail = F)
    message <- "alternative hypothesis: true difference in means is greater than 0"
  }
  else if (alternative == "less"){
    p1 <- t.test(x[which(!is.na(x)&!is.na(y))], y[which(!is.na(x)&!is.na(y))], alternative = "less", paired = T, var.equal=equal)$p.value
    p2 <- t.test(x[which(!is.na(x)&is.na(y))], y[which(is.na(x)&!is.na(y))], alternative = "less", var.equal=equal)$p.value
    Z1 <- qnorm(1-p1)
    Z2 <- qnorm(1-p2)
    p <- pnorm((w1*Z1 + w2*Z2)/sqrt(w1^2 + w2^2))
    message <- "alternative hypothesis: true difference in means is less than 0"
  }
  else{
    p1 <- t.test(x[which(!is.na(x)&!is.na(y))], y[which(!is.na(x)&!is.na(y))], alternative = "greater", paired = T, var.equal=equal)$p.value
    p2 <- t.test(x[which(!is.na(x)&is.na(y))], y[which(is.na(x)&!is.na(y))], alternative = "greater", var.equal=equal)$p.value
    Z1 <- qnorm(1-p1)
    Z2 <- qnorm(1-p2)
    p <- pnorm((w1*Z1 + w2*Z2)/sqrt(w1^2 + w2^2))
    if (p < 0.5){
      p <- 2*p
    }
    else{
      p <- 2*(1-p)
    }
    message <- "alternative hypothesis: true difference in means is not equal to 0"
  }
  
  cat("       ", "Liptak's Weighted Z-test\n\n", "p-value =", p, "\n", message, "\n", "number of matched:", n1, "\n", "number of partially matched pairs:", n2+n3, "\n", "weight of matched pairs:", w1, "\n", "weight of partially matched pairs", w2, "\n", "Z score of matched pairs:", Z1, "\n", "Z score of partially matched pairs:", Z2, "\n")
}
