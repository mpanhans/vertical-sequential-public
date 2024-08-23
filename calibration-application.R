##################################################
## Filename: calibration-application.R
## Project: Vertical Sequential Model
## Author: mtp
## Date: 4.11.2024
## Task: Empirical application with calibration
##################################################

## This script runs in about 10 minutes on a laptop with an AMD Ryzen 5 7640U
## processor and 16 GB of RAM.

rm(list=ls())  # clear workspace

library(BB)
library(rootSolve)
library(Matrix)
library(ggplot2)
library(beepr)
library(numDeriv)  # for jacobian() function
library(kableExtra)
library(dplyr)

options(knitr.kable.NA = '-')


##################################################
#### Load necessary functions
##################################################

source("my-functions.R")


##################################################
#### Define observable market outcomes
##################################################


#### observables
shares_obs <- c(.20,.20,.20,.20)

price_r <- c(10,10,10,10)
price_w <- c(4,4,4,4)

lambda <- c(0.5,0.5,0.5,0.5)

c_r_NA <- c(3, 3, NA, NA)

J <- length(shares_obs)

#### ownership
own_down_pre <- paste0("R",rep(c(1,2),each=2))
own_up_pre <- paste0("W",rep(c(1,2),2))
own_down_pre
own_up_pre

own_up_post <- own_up_pre
own_down_post <- own_down_pre
own_down_post[own_down_post == "R1"] <- "W1"


##################################################
#### Simultaneous model (Sheu and Taragin)
##################################################

#### calibrate demand
x00b <- 1.1
Bert_foc_calibrate_alpha2(param = x00b, 
                          own_down = own_down_pre, price = price_r, 
                          shares = shares_obs, cost  = c_r_NA, p_W = price_w )


out1a <- BBoptim(f = Bert_foc_calibrate_alpha2, par = x00b, 
                 own_down = own_down_pre, price = price_r, 
                 shares = shares_obs, cost  = c_r_NA, p_W = price_w)

alpha1 <- out1a$par
delta1 <- log(shares_obs) - log(1-sum(shares_obs)) + alpha1*price_r

alpha1
delta1

#### calibrate missing downstream costs

x00b <- rep(1.5,J)
Bert_foc_calibrate_costs(param = x00b, 
                         own_down = own_down_pre, price = price_r, 
                         shares = shares_obs, alpha = alpha1, delta = delta1,
                         p_W = price_w )


out1b <- BBoptim(f = Bert_foc_calibrate_costs, par = x00b, 
                 own_down = own_down_pre, price = price_r, 
                 shares = shares_obs, alpha = alpha1, delta = delta1,
                 p_W = price_w)

out1b$par
c_r_1 <- out1b$par

#### calibrate c_w

find_c_w <- multiroot(f = Barg_foc_gnl_cal_c_w, start = price_w*.8,
                      lambda=lambda,
                      p_W = price_w, own_down = own_down_pre,
                      own_up = own_up_pre,
                      alpha = alpha1, delta = delta1, 
                      c_R = c_r_1, p_R = price_r, sumFOC = FALSE)

c_w1 <- find_c_w$root


#### check outcome

tol <- .0001
error <- 1

p_W0 <- price_w * 1.1
p_R0 <- price_r * 1.1

while (error > tol) {
  
  out1 <- BBoptim(f = Bert_foc, par = p_R0, 
                  own_down = own_down_pre, alpha= alpha1, 
                  delta = delta1, cost = c_r_1,
                  p_W = p_W0, sumFOC = TRUE)
  
  p_R1 <- out1$par
  
  out2 <- BBoptim(par = as.numeric(p_W0), fn = Barg_foc, 
                  own_down = own_down_pre, own_up = own_up_pre, 
                  alpha= alpha1, delta = delta1, 
                  c_W = c_w1, c_R = c_r_1, lambda = lambda, 
                  p_R = p_R1)
  
  p_W1 <- out2$par
  
  error <- max(abs(c(p_W1-p_W0,p_R1-p_R0)))
  print(error)
  
  p_W0 <- p_W1
  p_R0 <- p_R1
}

price_r1 <- p_R1
price_w1 <- p_W1
print(price_r1)
print(price_w1)

shares1 <- (exp(delta1 - alpha1*price_r1))/(1+sum(exp(delta1 - alpha1*price_r1)))
as.numeric(shares1)



##################################################
#### Sequential model with linear pricing
##################################################

#### Assume no lump sum payment
#### Demand and cost parameters stay the same as above

#### Calibrate c_w. Use symmetry to get answer more quickly.
check2b <- optimize(f = Barg_NP_seq3_vert_cal_cw2, price_w=price_w,
                own_down=own_down_pre, own_up=own_up_pre,
                alpha=alpha1,delta=delta1,
                lambda=lambda, c_R=c_r_1,
                price_r=price_r,
                sigma=0, showAll = FALSE,
                setTol = 0.002,
                lower = 0, upper = 5.0)

print(check2b$minimum)

c_w_2 <- rep(check2b$minimum,J)  # If running code all the way through, use this


##################################################
#### Sequential model with two part tariff
##################################################


#### Calibrate c_w. Use symmetry to get answer more quickly.
check3b <- optimize(f = Barg_NP_seq3_vert_cal_cw2, price_w=price_w,
                    own_down=own_down_pre, own_up=own_up_pre,
                    alpha=alpha1,delta=delta1,
                    lambda=lambda, c_R=c_r_1,
                    price_r=price_r,
                    sigma=1.0, showAll = FALSE,
                    lower = 0, upper = 5.0)

print(check3b$minimum)

c_w_3 <- rep(check3b$minimum,J)  # If running code all the way through, use this


##################################################
#### Make table with merger simulations of all models
##################################################

alpha <- alpha1
delta <- delta1
c_R_vec <- c_r_1

#### First, the two sequential models

x_vals <- c(0,1.0)
numsims <- length(x_vals)
numsims1 <- numsims + 1
CS_pre <- as.vector(rep(0,numsims1))
CS_post <- as.vector(rep(0,numsims1))
p_R_pre <- matrix(data = 0, nrow = J, ncol = numsims1)
p_R_post <- matrix(data = 0, nrow = J, ncol = numsims1)
shares_pre <- matrix(data = 0, nrow = J, ncol = numsims1)
shares_post <- matrix(data = 0, nrow = J, ncol = numsims1)
p_W_pre <- matrix(data = 0, nrow = J, ncol = numsims1)
p_W_post <- matrix(data = 0, nrow = J, ncol = numsims1)


for (num in 1:numsims) {
  
  sigma_val <- x_vals[num]
  if (num == 1) {c_W_vec <- c_w_2}
  if (num == 2) {c_W_vec <- c_w_3}
  
  ## Pre-merger
  error <- rep(1,J)
  tol <- 0.01

  p_W0 <- price_w + .1
  p_R0 <- price_r
  iter <- 1
  
  while (max(error) > tol) {
    for (x in 1:J) {
      
      w_start <- p_W0[x]
      
      outtest <- optimize(f = Barg_NP_seq3_vert,
                          product_max = x, p_W = p_W0,
                          own_down = own_down_pre, own_up = own_up_pre, 
                          alpha= alpha, delta = delta, 
                          c_W = c_W_vec, c_R = c_R_vec, lambda = lambda, 
                          p_R0 = price_r, sigma = sigma_val, showAll = FALSE,
                          lower = 0, upper = 5.0)  
      
      pj_test <- outtest$minimum
      
      
      error[x] <- abs(pj_test - p_W0[x])
      print(p_W0)
      
      p_W0[x] <- pj_test
      
      # recover p_R at these p_W and update r_R0
      outtest_r <- BBoptim(f = Bert_foc_vert, par = p_R0, 
                           own_down = own_down_pre, own_up = own_up_pre,
                           alpha= alpha, 
                           delta = delta, c_R = c_R_vec,
                           p_W = p_W0, c_W = c_W_vec, sumFOC = TRUE,
                           control = list(trace=FALSE),
                           quiet = TRUE)
      
      p_R0 <- outtest_r$par
      
      iter <- iter + 1
      
    }
  }
  
  p_W1 <- p_W0
  p_R1 <- p_R0
  
  CS_pre[num] <- (1/alpha) * log(1 + sum(exp(delta - alpha * p_R1)) )
  p_R_pre[,num] <- p_R1
  p_W_pre[,num] <- p_W1
  shares_pre[,num] <- (exp(delta - alpha*p_R1))/(1+sum(exp(delta - alpha*p_R1)))

  ## Post-merger
  error <- rep(1,J)
  tol <- 0.01
  
  p_W0 <- p_W1
  p_R0 <- p_R1
  iter <- 1
  
  while (max(error) > tol) {
    for (x in 1:J) {
      
      w_start <- p_W0[x]
      
      outtest <- optimize(f = Barg_NP_seq3_vert,
                          product_max = x, p_W = p_W0,
                          own_down = own_down_post, own_up = own_up_post, 
                          alpha= alpha, delta = delta, 
                          c_W = c_W_vec, c_R = c_R_vec, lambda = lambda, 
                          p_R0 = p_R0, sigma = sigma_val, showAll = FALSE,
                          lower = 0, upper = 5.0)
      
      pj_test <- outtest$minimum
      
      error[x] <- abs(pj_test - p_W0[x])
      print(error)
      
      p_W0[x] <- pj_test
      
      # recover p_R at these p_W and update r_R0
      outtest_r <- BBoptim(f = Bert_foc_vert, par = p_R0, 
                           own_down = own_down_post, own_up = own_up_post,
                           alpha= alpha, 
                           delta = delta, c_R = c_R_vec,
                           p_W = p_W0, c_W = c_W_vec, sumFOC = TRUE,
                           control = list(trace=FALSE),
                           quiet = TRUE)
      
      p_R0 <- outtest_r$par
      
      iter <- iter + 1
      
    }
  }
  
  ## save values from the merger sim
  p_W1_post <- p_W0
  p_R1_post <- p_R0
  
  CS_post[num] <- (1/alpha) * log(1 + sum(exp(delta - alpha * p_R1_post)) )
  p_R_post[,num] <- p_R1_post
  p_W_post[,num] <- p_W1_post
  shares_post[,num] <- (exp(delta - alpha*p_R1_post))/(1+sum(exp(delta - alpha*p_R1_post)))
  
}


#### Next, the simultaneous model
num <- 3

## Pre-merger was done above
CS_pre[num] <- (1/alpha) * log(1 + sum(exp(delta - alpha * price_r1)) )
p_R_pre[,num] <- price_r1
p_W_pre[,num] <- price_w1
shares_pre[,num] <- (exp(delta - alpha*price_r1))/(1+sum(exp(delta - alpha*price_r1)))

## post merger

tol <- .0001
error <- 1

p_W0 <- price_w1
p_R0 <- price_r1

while (error > tol) {
  
  out1 <- BBoptim(f = Bert_foc_vert, par = p_R0, 
                  own_down = own_down_post, own_up = own_up_post,
                  alpha= alpha1, delta = delta1, c_R = c_r_1,
                  p_W = p_W0, c_W = c_w1, sumFOC = TRUE)
  
  p_R1 <- out1$par
  
  out2 <- BBoptim(par = as.numeric(p_W0), fn = Barg_foc_vert, 
                  own_down = own_down_post, own_up = own_up_post, 
                  alpha= alpha1, delta = delta1, 
                  c_W = c_w1, c_R = c_r_1, lambda = lambda, 
                  p_R = p_R1)
  
  p_W1 <- out2$par
  
  error <- max(abs(c(p_W1-p_W0,p_R1-p_R0)))
  print(error)
  
  p_W0 <- p_W1
  p_R0 <- p_R1
}

CS_post[num] <- (1/alpha) * log(1 + sum(exp(delta - alpha * p_R1)) )
p_R_post[,num] <- p_R1
p_W_post[,num] <- p_W1
shares_post[,num] <- (exp(delta - alpha*p_R1))/(1+sum(exp(delta - alpha*p_R1)))



#### Output table

df_tab1 <- data.frame( good = rep((1:J), times = numsims1) )

df_tab1$p_R_pre <- as.vector(  round(p_R_pre, 1) )
df_tab1$p_R_post <- as.vector( round(p_R_post, 1) )
df_tab1$p_W_pre <- as.vector(  round(p_W_pre, 1) )
df_tab1$p_W_post <- as.vector( round(p_W_post, 1) )
## VI good in output table:
df_tab1$p_W_post[c(1,5,9)] <- NA


df_tab1$shares_pre <- as.vector( round(shares_pre*100,1) )
df_tab1$shares_post <- as.vector( round(shares_post*100,1) )

df_tab1$pct_change_CS <- rep(round((CS_post - CS_pre)/ CS_pre *100, 1), each = J)
df_tab1$pct_change_p_R <- as.vector( round((p_R_post-p_R_pre)/p_R_pre*100, 1) )

