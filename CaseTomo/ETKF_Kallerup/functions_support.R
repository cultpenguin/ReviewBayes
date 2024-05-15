loc_entkf <- function(real_obs, r_matrix, y_matrix, x_matrix, x_avg, y_avg, n_e){
  # This function is following the algorithm in the paper, adapted to implemention in R
  r_inv <- solve(r_matrix)
  p_tilda_mx <- solve(t(y_matrix)%*%r_inv%*%y_matrix + (n_e-1)*diag(n_e))
  k_mx <- (x_matrix)%*%p_tilda_mx%*%t(y_matrix)%*%r_inv
  x_a <- x_avg + k_mx%*%(real_obs-y_avg)
  p_tilda_mx_not_inv <- t(y_matrix)%*%r_inv%*%y_matrix + (n_e-1)*diag(n_e)
  temp_svd <- svd(p_tilda_mx_not_inv)
  mU <- temp_svd$u
  mV <- temp_svd$v
  mS <- diag(temp_svd$d)
  vG <- 1/sqrt(mS)
  
  
  # double diag() to get it into a matrix format
  mG <- mV %*% diag(diag(vG)) %*% t(mV)
  
  x_a_mx <- x_matrix%*%(sqrt(n_e-1)*mG)
  
  #print(c("dim output:letkf",dim(x_a), dim(x_a_mx)))
  
  x_a_vec <- matrix(rep(x_a,ncol(x_a_mx)), nrow = nrow(x_a), ncol = ncol(x_a_mx))

  return((x_a_vec+x_a_mx))
  
}

# Function which is finding the intercept and the gradient corresponding to each of 
# the pair of senders and recievers that have their trace recorded
find.line <- function(p) {
  # (x,y) - sender
  p1 <- c(p[1],p[2])
  # (x,y) - reciever
  p2 <- c(p[3],p[4])
  x <- c(p1[1], p2[1])
  y <- c(p1[2], p2[2])
  fit <- as.numeric(coef(lm(y ~ x)))
  # b - intercept, m - gradient
  return(c(b=fit[1], m=fit[2]))
}