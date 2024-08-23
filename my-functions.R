#### Functions needed for simulations.R and calibration-application.R


##################################################
#### Define functions for equilibrium
##################################################


Bert_foc <- function(p,own_down,alpha,delta,cost,p_W,sumFOC = FALSE){
  # first create ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  ownership <- t(sapply(own_down, own_fun_down) )
  
  # then calculate foc
  shares <- (exp(delta - alpha*p))/(1+sum(exp(delta - alpha*p)))
  m <- p - cost - p_W
  ownd <- -alpha*shares*(1-shares)
  crossd <- alpha*shares%*%t(shares)
  dd <- crossd
  diag(dd) <- ownd
  omega <- (ownership * t(dd))
  
  foc <- omega %*% m + shares
  
  if (sumFOC == FALSE) {
    return(foc)
  } else {
    out <- sum(foc^2)
    return(out)
  }
  
}


Bert_foc_vert <- function(p,own_down,own_up,alpha,delta,
                          c_R,p_W,c_W, sumFOC = FALSE){
  
  # first create ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  own_R <- t(sapply(own_down, own_fun_down) )
  
  own_fun_up <- function(x) {as.numeric(x == own_up)}
  own_W <- t(sapply(own_up, own_fun_up) )
  
  VI_D <- as.numeric(own_down == own_up)
  k <- which(VI_D == 1)
  
  # VI_id identifies integrated firm. Should be singleton or empty:
  VI_id <- unique(own_down[which(own_down %in% own_up)])
  if (length(VI_id) > 1) { stop('Function cannot handle multiple integrated firms')}
  if (length(VI_id) == 0) { VI_id <- ""}
  
  VI_U_w <- as.numeric(VI_id == own_up) - VI_D
  VI_U_r <- as.numeric(VI_id == own_down) - VI_D
  
  # calculate Bertrand FOCs
  shares <- (exp(delta - alpha*p))/(1+sum(exp(delta - alpha*p)))
  m <- (1-VI_D) * (p - p_W - c_R) + (VI_D) * (p - c_W - c_R) #EDM Effect
  ownd <- -alpha*shares*(1-shares)
  crossd <- alpha*shares%*%t(shares)
  dd <- crossd
  diag(dd) <- ownd
  out <- (own_R * t(dd)) %*% m + shares
  
  UPP_effect <- (VI_D+VI_U_r) * (dd %*% (VI_U_w*(p_W - c_W)) )
  
  foc <- out + UPP_effect
  if (sumFOC == FALSE) {
    return(foc)
  } else {
    out <- sum(foc^2)
    return(out)
  }
}



Barg_foc <- function(p_W,own_down,own_up,alpha,delta,c_W,c_R,lambda,p_R,
                     returnGFT = FALSE){
  J <- length(p_W)
  
  # construct ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  own_R <- t(sapply(own_down, own_fun_down) )
  
  own_fun_up <- function(x) {as.numeric(x == own_up)}
  own_W <- t(sapply(own_up, own_fun_up) )
  
  # calculate shares
  shares <- (exp(delta - alpha*p_R))/(1+sum(exp(delta - alpha*p_R)))
  
  # counterfactual shares
  shares_tilde <- vector("list",J) 
  
  denom_tilde <- (1-diag(J)) %*% exp(delta - alpha*p_R)
  for (j in (1:J)) {
    shares_tilde[[j]] <- (exp(delta - alpha*p_R))/(1+denom_tilde[j])
    shares_tilde[[j]][j] <- 0
  }
  shares_tilde <- matrix(unlist(shares_tilde), ncol = J, byrow = FALSE)
  
  pi_w <- own_W %*% ((p_W - c_W)*shares)
  
  pi_w_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_W_tilde <- own_W
    own_W_tilde[j,j] <- 0
    temp <- own_W_tilde %*% ((p_W - c_W)*shares_tilde[,j])
    pi_w_tilde[[j]] <- temp[j]
  }
  
  pi_r <- own_R %*% ((p_R - p_W - c_R)*shares)
  
  pi_r_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_R_tilde <- own_R
    own_R_tilde[j,j] <- 0
    temp <- own_R_tilde %*% ((p_R - p_W - c_R)*shares_tilde[,j])
    pi_r_tilde[[j]] <- temp[j]
  }
  
  foc <- lambda*(1)*(pi_w-pi_w_tilde) - (1-lambda)*(pi_r - pi_r_tilde)*(1)
  if (returnGFT == FALSE) {
    return(sum(foc^2))
  } else {
    return(list("r_gft" = (pi_r - pi_r_tilde),"w_gft" = (pi_w-pi_w_tilde)))
  }
  
}


Barg_foc_vert <- function(p_W,own_down,own_up,alpha,delta,c_W,c_R,lambda,p_R,sumFOC = TRUE){
  
  J <- length(p_W)
  
  # first create ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  own_R <- t(sapply(own_down, own_fun_down) )
  
  own_fun_up <- function(x) {as.numeric(x == own_up)}
  own_W <- t(sapply(own_up, own_fun_up) )
  
  VI_D <- as.numeric(own_down == own_up)
  k <- which(VI_D == 1)
  
  # VI_id identifies integrated firm. Should be singleton or empty:
  VI_id <- unique(own_down[which(own_down %in% own_up)])
  if (length(VI_id) > 1) { stop('Function cannot handle multiple integrated firms')}
  if (length(VI_id) == 0) { VI_id <- ""}
  
  VI_U_w <- as.numeric(VI_id == own_up) - VI_D
  VI_U_r <- as.numeric(VI_id == own_down) - VI_D
  
  # calculate shares and disagreement shares
  shares <- (exp(delta - alpha*p_R))/(1+sum(exp(delta - alpha*p_R)))
  
  shares_tilde <- vector("list",J)
  
  denom_tilde <- (1-diag(J)) %*% exp(delta - alpha*p_R)
  for (j in (1:J)) {
    shares_tilde[[j]] <- (exp(delta - alpha*p_R))/(1+denom_tilde[j])
    shares_tilde[[j]][j] <- 0
  }
  shares_tilde <- matrix(unlist(shares_tilde), ncol = J, byrow = FALSE)
  
  # define margins, depends on whether vertically integrated
  margin_up <-   (1-VI_D)*(p_W - c_W)       + VI_D*(p_R - c_W - c_R)
  margin_down <- (1-VI_D)*(p_R - p_W - c_R) + VI_D*(p_R - c_W - c_R)
  
  # specify payoffs and disagreement payoffs
  pi_w <- own_W %*% (margin_up*shares)
  
  pi_w_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_W_tilde <- own_W
    own_W_tilde[j,j] <- 0
    temp <- own_W_tilde %*% (margin_up*shares_tilde[,j])
    pi_w_tilde[[j]] <- temp[j]
  }
  
  pi_r <- own_R %*% (margin_down*shares)
  
  pi_r_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_R_tilde <- own_R
    own_R_tilde[j,j] <- 0
    temp <- own_R_tilde %*% (margin_down*shares_tilde[,j])
    pi_r_tilde[[j]] <- temp[j]
  }
  
  # RRC Effect 
  RRC_effect2 <- t(shares_tilde - replicate(J,shares)) %*% (VI_U_r * margin_down)
  RRC_effect <- VI_U_w * RRC_effect2
  
  # Recapture Leverage Effect
  Recap_effect <- t(shares_tilde - replicate(J,shares)) %*% (VI_U_w * margin_up)
  Recap_effect <- VI_U_r * Recap_effect
  
  foc <- lambda*(1)*(pi_w - pi_w_tilde - RRC_effect) - 
    (1-lambda)*(pi_r - pi_r_tilde - Recap_effect)*(1)
  
  foc[VI_D] <- 0  # set foc value for integrated goods to zero
  
  if (sumFOC == FALSE) {
    return(foc)
  } else {
    out <- sum(foc^2)
    return(out)
  }
}


##################################################
#### Functions that allow for GNL
##################################################


share_calc <- function(p,delta,alpha){
  out <- (exp(delta - alpha*p))/(1+sum(exp(delta - alpha*p)))
  return(out)
}

share_calc_gnl <- function(p,delta,alpha,a_jk=NA,B=NA,mu=NA,
                           returnLogsum=FALSE){
  # a is a J-by-K matrix of allocation parameters
  # B is a J-by-K matrix of indicators designating nests
  # mu is a vector length K of nesting parameters
  
  # If no GNL parameters, treat as standard logit. One nest. mu=1.
  J <- length(p)
  
  if (any(is.na(B))) {
    K <- 1
    B <- matrix(1, ncol = 1, nrow = J)
    a_jk <- B
    mu <- rep(1,K)
  }
  
  J <- dim(B)[1]
  K <- dim(B)[2]
  
  V_j <- delta - alpha * p
  
  temp1 <- a_jk * matrix(rep(exp(V_j),K), ncol = K, nrow = J)
  temp2 <- temp1 ^ matrix(rep((1/mu),J), ncol = K, nrow = J, byrow = TRUE)
  P_k_num <- colSums(temp2) ^ mu
  #P_k_denom <- sum(P_k_num)  # NO outside option
  P_k_denom <- 1+sum(P_k_num)  # WITH outside option
  P_k <- P_k_num / P_k_denom   
  
  P_j_Bk <- temp2 / matrix(rep(colSums(temp2),J), ncol = K, nrow = J, byrow = TRUE)
  # this should be a JxK matrix with columns that sum to 1, unless empty nest
  # fix for empty nests:
  P_j_Bk[is.na(P_j_Bk)] <- 0
  
  P_j <- rowSums(P_j_Bk * matrix(rep((P_k),J), ncol = K, nrow = J, byrow = TRUE))
  # P_j should be a vector of length J that sums to one.
  
  if (returnLogsum == FALSE) {
    out <- P_j
    return(out)
  }
  if (returnLogsum == TRUE) {
    out <- log(P_k_denom)
    return(out)
  }
  
}


Bert_foc_vert_gnl <- function(p,own_down,own_up,alpha,delta,c_R,p_W,c_W,
                              a_jk=NA, B=NA, mu=NA, sumFOC = FALSE){
  
  # If no GNL parameters, treat as standard logit. One nest. mu=1.
  J <- length(p)
  if (any(is.na(B))) {
    K <- 1
    B <- matrix(1, ncol = 1, nrow = J)
    a_jk <- B
    mu <- rep(1,K)
  }
  
  # first create ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  own_R <- t(sapply(own_down, own_fun_down) )
  
  own_fun_up <- function(x) {as.numeric(x == own_up)}
  own_W <- t(sapply(own_up, own_fun_up) )
  
  VI_D <- as.numeric(own_down == own_up)
  k <- which(VI_D == 1)
  
  # VI_id identifies integrated firm. Should be singleton or empty:
  VI_id <- unique(own_down[which(own_down %in% own_up)])
  if (length(VI_id) > 1) { stop('Function cannot handle multiple integrated firms')}
  if (length(VI_id) == 0) { VI_id <- ""}
  
  VI_U_w <- as.numeric(VI_id == own_up) - VI_D
  VI_U_r <- as.numeric(VI_id == own_down) - VI_D
  
  # calculate Bertrand FOCs
  shares <- share_calc_gnl(p = p, alpha = alpha, delta = delta,
                           a_jk = a_jk, B=B, mu = mu)
  m <- (1-VI_D) * (p - p_W - c_R) + (VI_D) * (p - c_W - c_R) #EDM Effect
  dd <- jacobian(share_calc_gnl, x = p, delta = delta, alpha = alpha,
                 a_jk=a_jk,B=B,mu=mu)
  omega <- (own_R * t(dd))
  
  # if all costs aviailable:
  foc1 <- omega %*% m + shares
  
  UPP_effect <- (VI_D+VI_U_r) * (dd %*% (VI_U_w*(p_W - c_W)) )
  
  foc <- foc1 + UPP_effect
  if (sumFOC == FALSE) {
    return(foc)
  } else {
    out <- sum(foc^2)
    return(out)
  }
}



##################################################
#### Functions for calibrating
##################################################

Bert_foc2 <- function(p,own_down,alpha,delta,cost,p_W,sumFOC = FALSE){
  # first create ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  ownership <- t(sapply(own_down, own_fun_down) )
  
  # then calculate foc
  shares <- (exp(delta - alpha*p))/(1+sum(exp(delta - alpha*p)))
  m <- p - cost - p_W
  exl_na <- which(!is.na(m))
  
  ownd <- -alpha*shares*(1-shares)
  crossd <- alpha*shares%*%t(shares)
  dd <- crossd
  diag(dd) <- ownd
  omega <- (ownership * t(dd))
  
  # This to deal with possible missing/NA costs
  margin_tilde <- m[exl_na]
  omega_tilde <- omega[exl_na,exl_na]
  
  foc <- omega_tilde %*% margin_tilde + shares[exl_na]
  
  if (sumFOC == FALSE) {
    return(foc)
  } else {
    out <- sum(foc^2)
    return(out)
  }
}

Bert_foc_calibrate_alpha2 <- function(param,own_down,price,shares,cost,p_W){
  
  J <- length(price)
  alpha <- param[1]
  
  s_0 <- 1 - sum(shares)
  delta_j <- log(shares) - log(s_0) + alpha*price
  
  x0 <- price
  out1 <- BBoptim(f = Bert_foc2, par = x0, 
                  own_down = own_down, alpha= alpha, 
                  delta = delta_j, cost = cost, 
                  p_W = p_W, sumFOC = TRUE)
  
  p1 <- out1$par
  share1 <- (exp(delta_j - alpha*p1))/(1+sum(exp(delta_j - alpha*p1)))
  
  pdiff <- price - p1
  
  objfxn <- c(pdiff) %*% diag(J) %*% c(pdiff)
  return(objfxn)
}



Bert_foc_calibrate_costs <- function(param,own_down,price,shares,alpha,delta,p_W){
  
  J <- length(price)
  cost <- param
  
  s_0 <- 1 - sum(shares)
  delta_j <- log(shares) - log(s_0) + alpha*price
  
  x0 <- price
  out1 <- BBoptim(f = Bert_foc2, par = x0, 
                  own_down = own_down, alpha= alpha, 
                  delta = delta_j, cost = cost, 
                  p_W = p_W, sumFOC = TRUE)
  
  p1 <- out1$par
  share1 <- (exp(delta_j - alpha*p1))/(1+sum(exp(delta_j - alpha*p1)))
  
  pdiff <- price - p1
  
  objfxn <- c(pdiff) %*% diag(J) %*% c(pdiff)
  return(objfxn)
}



#### Next functions are for calibrating upstream costs in simultaneous model
#### with GNL demand.

Barg_foc_gnl_cal_c_w <- function(c_W,lambda,p_W,own_down,own_up,alpha,delta,c_R,p_R,
                                 a_jk=NA, B=NA, mu=NA, sumFOC = FALSE){
  
  # If no GNL parameters, treat as standard logit. One nest. mu=1.
  J <- length(p_W)
  if (any(is.na(B))) {
    K <- 1
    B <- matrix(1, ncol = 1, nrow = J)
    a_jk <- B
    mu <- rep(1,K)
  }
  
  # construct ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  own_R <- t(sapply(own_down, own_fun_down) )
  
  own_fun_up <- function(x) {as.numeric(x == own_up)}
  own_W <- t(sapply(own_up, own_fun_up) )
  
  # calculate shares
  shares <- share_calc_gnl(p = p_R, alpha = alpha, delta = delta,
                           a_jk = a_jk, B=B, mu = mu)
  
  # counterfactual shares
  shares_tilde <- vector("list",J) 
  
  for (i in (1:J)) {
    delta_cf <- delta
    delta_cf[i] <- -Inf
    cf_sharei <- share_calc_gnl(p = p_R, alpha = alpha, delta = delta_cf,
                                a_jk = a_jk, B=B, mu = mu)
    shares_tilde[[i]] <- cf_sharei
  }
  
  shares_tilde <- matrix(unlist(shares_tilde), ncol = J, byrow = FALSE)
  
  pi_w <- own_W %*% ((p_W - c_W)*shares)
  
  pi_w_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_W_tilde <- own_W
    own_W_tilde[j,j] <- 0
    temp <- own_W_tilde %*% ((p_W - c_W)*shares_tilde[,j])
    pi_w_tilde[[j]] <- temp[j]
  }
  
  pi_r <- own_R %*% ((p_R - p_W - c_R)*shares)
  
  pi_r_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_R_tilde <- own_R
    own_R_tilde[j,j] <- 0
    temp <- own_R_tilde %*% ((p_R - p_W - c_R)*shares_tilde[,j])
    pi_r_tilde[[j]] <- temp[j]
  }
  
  foc <- lambda*(1)*(pi_w-pi_w_tilde) - (1-lambda)*(pi_r - pi_r_tilde)*(1)
  if (sumFOC == TRUE) {
    return(sum(foc^2))
  } else {
    return(foc)
  }
  
}


Barg_NP_seq3_vert <- function(w_start,product_max,p_W,own_down,own_up,alpha,delta,
                              c_W,c_R,lambda,p_R0,sigma,showAll = FALSE,
                              maxJointProfits = FALSE){
  
  p_W[product_max] <- w_start
  J <- length(p_W)
  
  # construct ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  own_R <- t(sapply(own_down, own_fun_down) )
  
  own_fun_up <- function(x) {as.numeric(x == own_up)}
  own_W <- t(sapply(own_up, own_fun_up) )
  
  # indicator of which goods are integrated
  VI_D <- as.numeric(own_down == own_up)
  VI_D_idx <- which(own_down == own_up)
  
  VI_U_w <- as.numeric(own_up %in% own_down) - VI_D
  VI_U_r <- as.numeric(own_down %in% own_up) - VI_D
  # matrix form
  own_R_up <- t(sapply(own_down, 
                       function(i) as.numeric(own_up %in% i) ))
  own_W_down <- t(sapply(own_up, 
                         function(i) as.numeric(own_down %in% i) ))
  # remove integrated goods from partner profits
  own_R_up[,VI_D_idx] <- 0
  own_W_down[,VI_D_idx] <- 0
  
  # Given the inputs, calculate optimal retail prices
  outtest <- BBoptim(f = Bert_foc_vert, par = p_R0, 
                     own_down = own_down, own_up = own_up, 
                     alpha= alpha, 
                     delta = delta, c_R = c_R,
                     p_W = p_W, c_W = c_W, sumFOC = TRUE,
                     control = list(trace=FALSE),
                     quiet = TRUE)
  
  p_R <- outtest$par
  
  # calculate shares
  shares <- (exp(delta - alpha*p_R))/(1+sum(exp(delta - alpha*p_R)))
  
  # counterfactual shares
  shares_tilde <- matrix(data = 0, nrow = J, ncol = J)
  p_R_tilde <- matrix(data = 0, nrow = J, ncol = J)
  
  for (j in (1:J)) {
    delta_tilde <- delta
    delta_tilde[j] <- -Inf
    
    outtest_tilde <- BBoptim(f = Bert_foc_vert, par = p_R0, 
                             own_down = own_down, own_up = own_up, 
                             alpha= alpha, 
                             delta = delta_tilde, c_R = c_R,
                             p_W = p_W, c_W = c_W, sumFOC = TRUE,
                             control = list(trace=FALSE),
                             quiet = TRUE)
    
    p_R_tilde[,j] <- outtest_tilde$par
    
    shares_tilde[,j] <- share_calc(p=p_R_tilde[,j], delta = delta_tilde,
                                       alpha = alpha)
  }

  
  # define margins, depends on whether vertically integrated
  margin_up <-   (1-VI_D)*(p_W - c_W)       + VI_D*(p_R - c_W - c_R)
  margin_down <- (1-VI_D)*(p_R - p_W - c_R) + VI_D*(p_R - c_W - c_R)
  
  # specify payoffs and disagreement payoffs
  margin_up_tilde <- matrix(data = 0, nrow = J, ncol = J)
  margin_down_tilde <- matrix(data = 0, nrow = J, ncol = J)
  for (j in (1:J)) {
    margin_up_tilde[,j] <-   (1-VI_D)*(p_W - c_W)       + VI_D*(p_R_tilde[,j] - c_W - c_R)
    margin_down_tilde[,j] <- (1-VI_D)*(p_R_tilde[,j] - p_W - c_R) + VI_D*(p_R_tilde[,j] - c_W - c_R)
  }
  
  pi_w <- own_W %*% (margin_up*shares) + own_W_down %*% (margin_down*shares)
  
  pi_w_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_W_tilde <- own_W
    own_W_tilde[j,j] <- 0
    temp <- own_W_tilde %*% (margin_up_tilde[,j]*shares_tilde[,j]) + 
      own_R %*% (VI_U_r*margin_down_tilde[,j]*shares_tilde[,j])
    pi_w_tilde[[j]] <- temp[j]
  }
  
  pi_r <- own_R %*% (margin_down*shares) + own_R_up %*% (margin_up*shares)
  
  pi_r_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_R_tilde <- own_R
    own_R_tilde[j,j] <- 0
    temp <- own_R_tilde %*% (margin_down_tilde[,j]*shares_tilde[,j]) + 
      own_W %*% (VI_U_w*margin_up_tilde[,j]*shares_tilde[,j])
    pi_r_tilde[[j]] <- temp[j]
  }
  
  
  NP <- (pi_r - pi_r_tilde)^lambda * (pi_w-pi_w_tilde)^(1-lambda)
  L <- lambda^lambda * (1-lambda)^(1-lambda)
  NP_tpt <- L * (pi_r - pi_r_tilde + pi_w - pi_w_tilde)

  # Option that maximizes joint profits
  if (maxJointProfits == TRUE) {
    NP_tpt <- (pi_r + pi_w)
  }
  
  if (showAll == TRUE) {
    out <- (1-sigma) * NP + sigma * NP_tpt
    F_j <- (1-lambda) * (pi_r - pi_r_tilde) - lambda * (pi_w - pi_w_tilde) 
    return(list("foc" = out,
                "r_gft" = (pi_r - pi_r_tilde),
                "w_gft" = (pi_w - pi_w_tilde),
                "F_j" = F_j) )
  } else {
    out <- -(1-sigma) * NP - sigma * NP_tpt
    out <- out[product_max]
    return(out)
  }
  
}  




Barg_NP_seq3_vert_gnl <- function(w_start,product_max,p_W,own_down,own_up,alpha,delta,
                                  c_W,c_R,lambda,p_R0,sigma,
                                  a_jk,B,mu,showAll = FALSE){
  
  p_W[product_max] <- w_start
  J <- length(p_W)
  
  # construct ownership matrices
  own_fun_down <- function(x) {as.numeric(x == own_down)}
  own_R <- t(sapply(own_down, own_fun_down) )
  
  own_fun_up <- function(x) {as.numeric(x == own_up)}
  own_W <- t(sapply(own_up, own_fun_up) )
  
  # indicator of which goods are integrated
  VI_D <- as.numeric(own_down == own_up)
  VI_D_idx <- which(own_down == own_up)
  
  # VI_U_w = all non-integrated goods sold by upstream divisions of integrated firms. Calculate by finding all upstream goods that are owned by an integrated firm, and subtract off the integrated goods.
  # VI_U_r = all non-integrated goods sold by downstream divisions of integrated firms. Calculate by finding all downstream goods that are sold by an integrated firm, and subtract off the integrated goods.
  VI_U_w <- as.numeric(own_up %in% own_down) - VI_D
  VI_U_r <- as.numeric(own_down %in% own_up) - VI_D
  # matrix form
  own_R_up <- t(sapply(own_down, 
                       function(i) as.numeric(own_up %in% i) ))
  own_W_down <- t(sapply(own_up, 
                         function(i) as.numeric(own_down %in% i) ))
  # remove integrated goods from partner profits
  own_R_up[,VI_D_idx] <- 0
  own_W_down[,VI_D_idx] <- 0
  
  # Given the inputs, calculate optimal retail prices
  outtest <- BBoptim(f = Bert_foc_vert_gnl, par = p_R0, 
                     own_down = own_down, own_up = own_up, 
                     alpha= alpha, 
                     delta = delta, c_R = c_R_vec,
                     p_W = p_W, c_W = c_W_vec, 
                     a_jk=a_jk, B=B, mu=mu, sumFOC = TRUE,
                     control = list(trace=FALSE),
                     quiet = TRUE)
  
  p_R <- outtest$par
  
  # calculate shares
  shares <- share_calc_gnl(p=p_R,delta,alpha,a_jk,B,mu)
  
  # counterfactual shares
  shares_tilde <- matrix(data = 0, nrow = J, ncol = J)
  p_R_tilde <- matrix(data = 0, nrow = J, ncol = J)
  
  for (j in (1:J)) {
    delta_tilde <- delta
    delta_tilde[j] <- -Inf
    
    outtest_tilde <- BBoptim(f = Bert_foc_vert, par = p_R0, 
                             own_down = own_down, own_up = own_up, 
                             alpha= alpha, 
                             delta = delta_tilde, c_R = c_R,
                             p_W = p_W, c_W = c_W, sumFOC = TRUE,
                             control = list(trace=FALSE),
                             quiet = TRUE)
    
    p_R_tilde[,j] <- outtest_tilde$par
    
    shares_tilde[,j] <- share_calc_gnl(p=p_R_tilde[,j], delta = delta_tilde,
                                   alpha = alpha, a_jk=a_jk, B=B, mu=mu)
  }
  
  
  # define margins, depends on whether vertically integrated
  margin_up <-   (1-VI_D)*(p_W - c_W)       + VI_D*(p_R - c_W - c_R)
  margin_down <- (1-VI_D)*(p_R - p_W - c_R) + VI_D*(p_R - c_W - c_R)
  
  # specify payoffs and disagreement payoffs
  margin_up_tilde <- matrix(data = 0, nrow = J, ncol = J)
  margin_down_tilde <- matrix(data = 0, nrow = J, ncol = J)
  for (j in (1:J)) {
    margin_up_tilde[,j] <-   (1-VI_D)*(p_W - c_W)       + VI_D*(p_R_tilde[,j] - c_W - c_R)
    margin_down_tilde[,j] <- (1-VI_D)*(p_R_tilde[,j] - p_W - c_R) + VI_D*(p_R_tilde[,j] - c_W - c_R)
  }
  
  pi_w <- own_W %*% (margin_up*shares) + own_W_down %*% (margin_down*shares)
  
  pi_w_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_W_tilde <- own_W
    own_W_tilde[j,j] <- 0
    temp <- own_W_tilde %*% (margin_up_tilde[,j]*shares_tilde[,j]) + 
      own_R %*% (VI_U_r*margin_down_tilde[,j]*shares_tilde[,j])
    pi_w_tilde[[j]] <- temp[j]
  }
  
  pi_r <- own_R %*% (margin_down*shares) + own_R_up %*% (margin_up*shares)
  
  pi_r_tilde <- vector("numeric",J)
  for (j in (1:J)) {
    own_R_tilde <- own_R
    own_R_tilde[j,j] <- 0
    temp <- own_R_tilde %*% (margin_down_tilde[,j]*shares_tilde[,j]) + 
      own_W %*% (VI_U_w*margin_up_tilde[,j]*shares_tilde[,j])
    pi_r_tilde[[j]] <- temp[j]
  }
  
  
  NP <- (pi_r - pi_r_tilde)^lambda * (pi_w-pi_w_tilde)^(1-lambda)
  L <- lambda^lambda * (1-lambda)^(1-lambda)
  NP_tpt <- L * (pi_r - pi_r_tilde + pi_w - pi_w_tilde)
  
  if (showAll == TRUE) {
    out <- (1-sigma) * NP + sigma * NP_tpt
    return(out)
  } else {
    out <- -(1-sigma) * NP - sigma * NP_tpt
    out <- out[product_max]
    return(out)
  }
  
}  



#### Calibration functions



##### cw2 calibrates costs assuming symmetric costs

Barg_NP_seq3_vert_cal_cw2 <- function(c_w_val2,price_w,own_down,
                                      own_up, alpha,delta,
                                      lambda,c_R,price_r,sigma,
                                      setTol = 0.01,
                                      setMaxIter = 500,
                                      showAll = FALSE){
  
  c_w_val <- rep(c_w_val2, length(delta))
  error <- rep(1,J)
  tol <- setTol
  
  p_W0 <- price_w
  p_R0 <- price_r
  iter <- 1
  
  while (max(error) > tol & iter < setMaxIter) {
    for (x in 1:J) {
      
      w_start <- p_W0[x]

      checktest <- optimize(f = Barg_NP_seq3_vert,
                            product_max = x, p_W = p_W0,
                            own_down = own_down, own_up = own_up, 
                            alpha= alpha, delta = delta, 
                            c_W = c_w_val, c_R = c_R, lambda = lambda, 
                            p_R0 = price_r, sigma = sigma, showAll = FALSE,
                            lower = 0, upper = 5)
      
      
      pw_test <- checktest$minimum
      
      error[x] <- abs(pw_test - p_W0[x])
      
      p_W0[x] <- pw_test
      
      
      # recover p_R at these p_W and update r_R0
      outtest_r <- BBoptim(f = Bert_foc_vert, par = p_R0, 
                           own_down = own_down, own_up = own_up,
                           alpha= alpha, 
                           delta = delta, c_R = c_R,
                           p_W = p_W0, c_W = c_w_val, sumFOC = TRUE,
                           control = list(trace=FALSE),
                           quiet = TRUE)
      
      p_R0 <- outtest_r$par
      
      iter <- iter + 1
      if (iter > setMaxIter) {warning("Max iterations reached")}
    }
  }
  
  price_w2 <- p_W0
  price_r2 <- p_R0
  
  shares2 <- (exp(delta1 - alpha1*price_r2))/(1+sum(exp(delta1 - alpha1*price_r2)))
  
  if (showAll == TRUE) {
    return(list("price_w" = price_w2,"price_r" = price_r2,"shares" = as.numeric(shares2)) )
  } else {
    out <- sum((price_w2 - price_w)^2)
    return(out)
  }
  
}  





