##################################################
## Filename: simulations.R
## Project: Vertical Sequential Model
## Author: mtp
## Date: 8.21.2023
## Task: Model Simulations
##################################################

## This script runs in about 17 minutes on a laptop with an AMD Ryzen 5 7640U
## processor and 16 GB of RAM.

rm(list=ls())  # clear workspace

library(BB)
library(rootSolve)
library(Matrix)
library(ggplot2)
library(beepr)
library(numDeriv)  # for jacobian() function
library(kableExtra)
library(dplyr)   # this is needed in conjunction with kable

options(knitr.kable.NA = '-')

##################################################
#### Load necessary functions
##################################################

source("my-functions.R")


##################################################
#### Define market primitives
##################################################

# structural parameters
alpha  <- 1.2
R <- 2
W <- 2
delta <- c(3.5, 3.5, 3.5, 3.5)
c_R_vec <- matrix(.20, nrow = (R*W), ncol = 1)
c_W_vec <- matrix(.20, nrow = (R*W), ncol = 1)

lambda <- 0.5
J <- length(c_R_vec)

# ownership
own_down_pre <- paste0("R",rep(c(1,2),each=2))
own_up_pre <- paste0("W",rep(c(1,2),2))
own_down_pre
own_up_pre

own_up_post <- own_up_pre
own_down_post <- own_down_pre
own_down_post[own_down_post == "R1"] <- "W1"


##################################################
#### 1. Baseline Model Figure
##################################################

x_vals <- c(0,0.2,0.4,0.6,0.8,1.0)

numsims <- length(x_vals)
CS_pre <- as.vector(rep(0,numsims))
CS_post <- as.vector(rep(0,numsims))
p_R_pre <- matrix(data = 0, nrow = J, ncol = numsims)
p_R_post <- matrix(data = 0, nrow = J, ncol = numsims)
shares_pre <- matrix(data = 0, nrow = J, ncol = numsims)
shares_post <- matrix(data = 0, nrow = J, ncol = numsims)
convflag_pre <- as.vector(rep(0,numsims))
convflag_post <- as.vector(rep(0,numsims))
profits_pre <- as.vector(rep(0,numsims))
profits_post <- as.vector(rep(0,numsims))
p_W_pre <- matrix(data = 0, nrow = J, ncol = numsims)
p_W_post <- matrix(data = 0, nrow = J, ncol = numsims)

useBBoptim <- FALSE

for (num in 1:numsims) {
  print("Scenario")
  print(num)
  sigma_val <- x_vals[num]
  
  ## Pre-merger
  error <- rep(1,J)
  tol <- 0.01
  
  p_W0 <- c_W_vec*2
  p_R0 <- as.numeric((c_R_vec+c_W_vec)*3)
  iter <- 1
  
  while (max(error) > tol) {
    for (x in 1:J) {
      
      w_start <- p_W0[x]
      
      if (useBBoptim == FALSE) {
      outtest <- optimize(f = Barg_NP_seq3_vert,
                         product_max = x, p_W = p_W0,
                         own_down = own_down_pre, own_up = own_up_pre, 
                         alpha= alpha, delta = delta, 
                         c_W = c_W_vec, c_R = c_R_vec, lambda = lambda, 
                         p_R0 = p_R0, sigma = sigma_val, showAll = FALSE,
                         maxJointProfits = FALSE,
                         lower = 0, upper = 2.0)
      
      pj_test <- outtest$minimum
      }
      if (useBBoptim == TRUE) {
        outtest <- BBoptim(par = w_start, fn = Barg_NP_seq3_vert,
                            product_max = x, p_W = p_W0,
                            own_down = own_down_pre, own_up = own_up_pre, 
                            alpha= alpha, delta = delta, 
                            c_W = c_W_vec, c_R = c_R_vec, lambda = lambda, 
                            p_R0 = p_R0, sigma = sigma_val, showAll = FALSE)
        
        pj_test <- outtest$par
      }
      
      error[x] <- abs(pj_test - p_W0[x])
      print(error)
      
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
  
  convflag_pre[num] <- outtest_r$convergence
  profits_pre[num] <- (p_R1[1] - p_W1[1] - c_R_vec[1,]) * shares_pre[,num][1] +
    (p_R1[2] - p_W1[2] - c_R_vec[2,]) * shares_pre[,num][1] + 
    (p_W1[1] - c_W_vec[1]) * shares_pre[,num][1] + 
    (p_W1[3] - c_W_vec[3]) * shares_pre[,num][3]
  
  ## Post-merger
  error <- rep(1,J)
  tol <- 0.01
  
  p_W0 <- p_W1
  p_R0 <- p_R1
  iter <- 1
  
  while (max(error) > tol) {
    for (x in 1:J) {
      
      w_start <- p_W0[x]
      
      if (useBBoptim == FALSE) {
      outtest <- optimize(f = Barg_NP_seq3_vert,
                         product_max = x, p_W = p_W0,
                         own_down = own_down_post, own_up = own_up_post, 
                         alpha= alpha, delta = delta, 
                         c_W = c_W_vec, c_R = c_R_vec, lambda = lambda, 
                         p_R0 = p_R0, sigma = sigma_val, showAll = FALSE,
                         maxJointProfits = FALSE,
                         lower = 0, upper = 2.0)
      
      pj_test <- outtest$minimum
      }
      if (useBBoptim == TRUE) {
        outtest <- BBoptim(par = w_start, fn = Barg_NP_seq3_vert,
                            product_max = x, p_W = p_W0,
                            own_down = own_down_post, own_up = own_up_post, 
                            alpha= alpha, delta = delta, 
                            c_W = c_W_vec, c_R = c_R_vec, lambda = lambda, 
                            p_R0 = p_R0, sigma = sigma_val, showAll = FALSE)
        
        pj_test <- outtest$par
      }
      
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
  
  convflag_post[num] <- outtest_r$convergence
  
  profits_post[num] <- (p_R1_post[1] - p_W1_post[1] - c_R_vec[1,]) * shares_post[,num][1] +
    (p_R1_post[2] - p_W1_post[2] - c_R_vec[2,]) * shares_post[,num][1] + 
    (p_W1_post[1] - c_W_vec[1]) * shares_post[,num][1] + 
    (p_W1_post[3] - c_W_vec[3]) * shares_post[,num][3]
  
}


## Create the figure

df_fig1 <- data.frame(sigma = x_vals )

df_fig1$cs_pre <- CS_pre
df_fig1$cs_post <- CS_post
df_fig1$pct_change_CS <- (df_fig1$cs_post-df_fig1$cs_pre)/df_fig1$cs_pre*100
shares_pre; shares_post; df_fig1

plot1 <- ggplot(data=df_fig1, aes(x=sigma, y=pct_change_CS)) +
  geom_point() +
  xlab("Sigma") + ylab("Percent CS change")
plot1



# check convergence flags
convflag_pre
convflag_post

# check if mergers are profitable
profits_pre
profits_post
is_profitable <- (profits_post > profits_pre)
is_profitable


#### Output table

df_tab1 <- data.frame(sigma = rep(x_vals, each = J),
                      good = rep((1:J), times = numsims) )

df_tab1$p_R_pre <- as.vector(  round(p_R_pre, 1) )
df_tab1$p_R_post <- as.vector( round(p_R_post, 1) )
df_tab1$p_W_pre <- as.vector(  round(p_W_pre, 1) )
df_tab1$p_W_post <- as.vector( round(p_W_post, 1) )
## VI good in output table:
df_tab1$p_W_post[c(1,5,9,13,17,21)] <- NA

df_tab1$shares_pre <- as.vector( round(shares_pre*100,1) )
df_tab1$shares_post <- as.vector( round(shares_post*100,1) )

df_tab1$pct_change_CS <- rep(round(df_fig1$pct_change_CS, 1), each = J)
df_tab1$pct_change_p_R <- as.vector( round((p_R_post-p_R_pre)/p_R_pre*100, 1) )




##################################################
#### 2. Proportion integrated figure
##################################################

x_vals <- c(c(1.5,2.5,3.5,4.5,5.5),   
            c(1.5,2.5,3.5,4.5,5.5))
numsims <- length(x_vals)
sigma_values <- rep( c(0,1), each = length(x_vals)/2)

CS_pre <- as.vector(rep(0,numsims))
CS_post <- as.vector(rep(0,numsims))
p_R_pre <- matrix(data = 0, nrow = J, ncol = numsims)
p_R_post <- matrix(data = 0, nrow = J, ncol = numsims)
p_W_pre <- matrix(data = 0, nrow = J, ncol = numsims)
p_W_post <- matrix(data = 0, nrow = J, ncol = numsims)
shares_pre <- matrix(data = 0, nrow = J, ncol = numsims)
shares_post <- matrix(data = 0, nrow = J, ncol = numsims)
profits_pre <- as.vector(rep(0,numsims))
profits_post <- as.vector(rep(0,numsims))

for (num in 1:numsims) {
  print("Scenario")
  print(num)
  
  delta_val <- delta
  delta1 <- x_vals[num]
  delta_val[1] <- delta1
  sigma_val <- sigma_values[num]
  
  ## Pre-merger
  error <- rep(1,J)
  tol <- 0.01
  
  p_W0 <- c_W_vec*4
  p_R0 <- as.numeric((c_R_vec+c_W_vec)*6)
  iter <- 1
  
  while (max(error) > tol) {
    for (x in 1:J) {
      
      w_start <- p_W0[x]
      
      outtest <- optimize(f = Barg_NP_seq3_vert,
                         product_max = x, p_W = p_W0,
                         own_down = own_down_pre, own_up = own_up_pre, 
                         alpha= alpha, delta = delta_val, 
                         c_W = c_W_vec, c_R = c_R_vec, lambda = lambda, 
                         p_R0 = p_R0, sigma = sigma_val, showAll = FALSE,
                         lower = 0, upper = 2.0)
      
      pj_test <- outtest$minimum
      
      error[x] <- abs(pj_test - p_W0[x])
      print(error)
      
      p_W0[x] <- pj_test
      
      # recover p_R at these p_W and update r_R0
      outtest_r <- BBoptim(f = Bert_foc_vert, par = p_R0, 
                           own_down = own_down_pre, own_up = own_up_pre,
                           alpha= alpha, 
                           delta = delta_val, c_R = c_R_vec,
                           p_W = p_W0, c_W = c_W_vec, sumFOC = TRUE,
                           control = list(trace=FALSE),
                           quiet = TRUE)
      
      p_R0 <- outtest_r$par
      
      iter <- iter + 1
      
    }
  }
  
  p_W1 <- p_W0
  p_R1 <- p_R0
  
  CS_pre[num] <- (1/alpha) * log(1 + sum(exp(delta_val - alpha * p_R1)) )
  p_R_pre[,num] <- p_R1
  p_W_pre[,num] <- p_W1
  shares_pre[,num] <- (exp(delta_val - alpha*p_R1))/(1+sum(exp(delta_val - alpha*p_R1)))
  
  profits_pre[num] <- (p_R1[1] - p_W1[1] - c_R_vec[1,]) * shares_pre[,num][1] +
    (p_R1[2] - p_W1[2] - c_R_vec[2,]) * shares_pre[,num][1] + 
    (p_W1[1] - c_W_vec[1]) * shares_pre[,num][1] + 
    (p_W1[3] - c_W_vec[3]) * shares_pre[,num][3]
  
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
                         alpha= alpha, delta = delta_val, 
                         c_W = c_W_vec, c_R = c_R_vec, lambda = lambda, 
                         p_R0 = p_R0, sigma = sigma_val, showAll = FALSE,
                         lower = 0, upper = 2.0)
      
      pj_test <- outtest$minimum
      
      error[x] <- abs(pj_test - p_W0[x])
      print(error)
      
      p_W0[x] <- pj_test
      
      # recover p_R at these p_W and update r_R0
      outtest_r <- BBoptim(f = Bert_foc_vert, par = p_R0, 
                           own_down = own_down_post, own_up = own_up_post,
                           alpha= alpha, 
                           delta = delta_val, c_R = c_R_vec,
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
  
  CS_post[num] <- (1/alpha) * log(1 + sum(exp(delta_val - alpha * p_R1_post)) )
  p_R_post[,num] <- p_R1_post
  p_W_post[,num] <- p_W1_post
  shares_post[,num] <- (exp(delta_val - alpha*p_R1_post))/(1+sum(exp(delta_val - alpha*p_R1_post)))
  
  profits_post[num] <- (p_R1_post[1] - p_W1_post[1] - c_R_vec[1,]) * shares_post[,num][1] +
    (p_R1_post[2] - p_W1_post[2] - c_R_vec[2,]) * shares_post[,num][1] + 
    (p_W1_post[1] - c_W_vec[1]) * shares_post[,num][1] + 
    (p_W1_post[3] - c_W_vec[3]) * shares_post[,num][3]
  
}


#### Create the figure
df_fig2 <- data.frame(int_pre_share = (shares_pre[1,]*100) )
df_fig2$cs_pre <- CS_pre
df_fig2$cs_post <- CS_post
df_fig2$sigma <- rep(c("Linear pricing","Two-part tariff"), each = numsims/2)

df_fig2$pct_change_CS <- (df_fig2$cs_post-df_fig2$cs_pre)/df_fig2$cs_pre*100


plot2 <- ggplot(data=df_fig2, 
                aes(x=int_pre_share, y=pct_change_CS, group=sigma, color=sigma)) +
  geom_point() +
  xlab("Pre-Merger Integrated Good Share") + ylab("Percent CS change") +
  scale_color_brewer(palette="Dark2") +
  theme(legend.position="bottom")
plot2



# check if mergers are profitable
profits_pre
profits_post
is_profitable <- (profits_post > profits_pre)
is_profitable

#### Output table

df_tab2 <- data.frame(sigma = rep(sigma_values, each = J),
                      good = rep((1:J), times = numsims) )

df_tab2$p_R_pre <- as.vector(  round(p_R_pre, 1) )
df_tab2$p_R_post <- as.vector( round(p_R_post, 1) )
df_tab2$p_W_pre <- as.vector(  round(p_W_pre, 1) )
df_tab2$p_W_post <- as.vector( round(p_W_post, 1) )
## VI good in output table:
df_tab2$p_W_post[((1:10) * 4 - 3)] <- NA

df_tab2$shares_pre <- as.vector( round(shares_pre*100,1) )
df_tab2$shares_post <- as.vector( round(shares_post*100,1) )

df_tab2$pct_change_CS <- rep(round(df_fig2$pct_change_CS, 1), each = J)
df_tab2$pct_change_p_R <- as.vector( round((p_R_post-p_R_pre)/p_R_pre*100, 1) )







##################################################
#### 4. Size of outside option
##################################################

# x_vals is all delta's, to scale down oo
x_vals <- c(c(2.5,3.5,4.5,5.5,6.0)*0.9,   
            c(1.5,2.5,3.5,4.5,5.5) )
numsims <- length(x_vals)
sigma_values <- rep( c(0,1), each = length(x_vals)/2)

#### Simulations

CS_pre <- as.vector(rep(0,numsims))
CS_post <- as.vector(rep(0,numsims))
p_R_pre <- matrix(data = 0, nrow = J, ncol = numsims)
p_R_post <- matrix(data = 0, nrow = J, ncol = numsims)
p_W_pre <- matrix(data = 0, nrow = J, ncol = numsims)
p_W_post <- matrix(data = 0, nrow = J, ncol = numsims)
shares_pre <- matrix(data = 0, nrow = J, ncol = numsims)
shares_post <- matrix(data = 0, nrow = J, ncol = numsims)
profits_pre <- as.vector(rep(0,numsims))
profits_post <- as.vector(rep(0,numsims))

for (num in 1:numsims) {
  print("Scenario")
  print(num)
  
  delta_val <- rep(x_vals[num], J)
  sigma_val <- sigma_values[num]
  
  ## Pre-merger
  error <- rep(1,J)
  tol <- 0.01
  
  p_W0 <- c_W_vec*2
  p_R0 <- as.numeric((c_R_vec+c_W_vec)*3)
  iter <- 1
  
  while (max(error) > tol) {
    for (x in 1:J) {
      
      w_start <- p_W0[x]
      
      outtest <- optimize(f = Barg_NP_seq3_vert,
                         product_max = x, p_W = p_W0,
                         own_down = own_down_pre, own_up = own_up_pre, 
                         alpha= alpha, delta = delta_val, 
                         c_W = c_W_vec, c_R = c_R_vec, lambda = lambda, 
                         p_R0 = p_R0, sigma = sigma_val, showAll = FALSE,
                         lower = 0, upper = 2.0)
      
      pj_test <- outtest$minimum
      
      error[x] <- abs(pj_test - p_W0[x])
      print(error)
      
      p_W0[x] <- pj_test
      
      # recover p_R at these p_W and update r_R0
      outtest_r <- BBoptim(f = Bert_foc_vert, par = p_R0, 
                           own_down = own_down_pre, own_up = own_up_pre,
                           alpha= alpha, 
                           delta = delta_val, c_R = c_R_vec,
                           p_W = p_W0, c_W = c_W_vec, sumFOC = TRUE,
                           control = list(trace=FALSE),
                           quiet = TRUE)
      
      p_R0 <- outtest_r$par
      
      iter <- iter + 1
      
    }
  }
  
  p_W1 <- p_W0
  p_R1 <- p_R0
  
  CS_pre[num] <- (1/alpha) * log(1 + sum(exp(delta_val - alpha * p_R1)) )
  p_R_pre[,num] <- p_R1
  p_W_pre[,num] <- p_W1
  shares_pre[,num] <- (exp(delta_val - alpha*p_R1))/(1+sum(exp(delta_val - alpha*p_R1)))
  
  profits_pre[num] <- (p_R1[1] - p_W1[1] - c_R_vec[1,]) * shares_pre[,num][1] +
    (p_R1[2] - p_W1[2] - c_R_vec[2,]) * shares_pre[,num][1] + 
    (p_W1[1] - c_W_vec[1]) * shares_pre[,num][1] + 
    (p_W1[3] - c_W_vec[3]) * shares_pre[,num][3]
  
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
                         alpha= alpha, delta = delta_val, 
                         c_W = c_W_vec, c_R = c_R_vec, lambda = lambda, 
                         p_R0 = p_R0, sigma = sigma_val, showAll = FALSE,
                         lower = 0, upper = 2.0)
      
      pj_test <- outtest$minimum
      
      error[x] <- abs(pj_test - p_W0[x])
      print(error)
      
      p_W0[x] <- pj_test
      
      # recover p_R at these p_W and update r_R0
      outtest_r <- BBoptim(f = Bert_foc_vert, par = p_R0, 
                           own_down = own_down_post, own_up = own_up_post,
                           alpha= alpha, 
                           delta = delta_val, c_R = c_R_vec,
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
  
  CS_post[num] <- (1/alpha) * log(1 + sum(exp(delta_val - alpha * p_R1_post)) )
  p_R_post[,num] <- p_R1_post
  p_W_post[,num] <- p_W1_post
  shares_post[,num] <- (exp(delta_val - alpha*p_R1_post))/(1+sum(exp(delta_val - alpha*p_R1_post)))
  
  profits_post[num] <- (p_R1_post[1] - p_W1_post[1] - c_R_vec[1,]) * shares_post[,num][1] +
    (p_R1_post[2] - p_W1_post[2] - c_R_vec[2,]) * shares_post[,num][1] + 
    (p_W1_post[1] - c_W_vec[1]) * shares_post[,num][1] + 
    (p_W1_post[3] - c_W_vec[3]) * shares_post[,num][3]
}

#### Create the figure
df_fig4 <- data.frame(oo_share = ((1-colSums(shares_pre))*100) )
df_fig4$cs_pre <- CS_pre
df_fig4$cs_post <- CS_post
df_fig4$sigma <- rep(c("Linear pricing","Two-part tariff"), each = numsims/2)

df_fig4$pct_change_CS <- (df_fig4$cs_post-df_fig4$cs_pre)/df_fig4$cs_pre*100

plot4 <- ggplot(data=df_fig4, 
                aes(x=oo_share, y=pct_change_CS, group=sigma, color=sigma)) +
  geom_point() +
  xlab("Pre-Merger Outside Option Share") + ylab("Percent CS change") +
  scale_color_brewer(palette="Dark2") +
  theme(legend.position="bottom")
plot4


# check if mergers are profitable
profits_pre
profits_post
is_profitable <- (profits_post > profits_pre)
is_profitable


#### Output table

df_tab4 <- data.frame(sigma = rep(sigma_values, each = J),
                      good = rep((1:J), times = numsims) )

df_tab4$p_R_pre <- as.vector(  round(p_R_pre, 1) )
df_tab4$p_R_post <- as.vector( round(p_R_post, 1) )
df_tab4$p_W_pre <- as.vector(  round(p_W_pre, 1) )
df_tab4$p_W_post <- as.vector( round(p_W_post, 1) )
## VI good in output table:
df_tab4$p_W_post[((1:10) * 4 - 3)] <- NA

df_tab4$shares_pre <- as.vector( round(shares_pre*100,1) )
df_tab4$shares_post <- as.vector( round(shares_post*100,1) )

df_tab4$pct_change_CS <- rep(round(df_fig4$pct_change_CS, 1), each = J)
df_tab4$pct_change_p_R <- as.vector( round((p_R_post-p_R_pre)/p_R_pre*100, 1) )





##################################################
#### 5. Comparison to simultaneous equilibrium (S&T)
##################################################

c_R_vec <- as.numeric(c_R_vec)
c_W_vec <- as.numeric(c_W_vec)

## 5.1 and 5.2.
## Here, loop through x_vals (lambda barg weights) twice times Once for sigma=0, then
## for sigma=1. Then separately for simultaneous code.

x_vals1 <- rep(c(0.30,.40,0.50,0.60,.70), 2)
x_vals2 <- c(0.30,.40,0.50,0.60,.70)
numsims1 <- length(x_vals1)
numsims2 <- length(x_vals2)
sigma_values <- rep( c(0,1), each = length(x_vals1)/2)


CS_pre <- as.vector(rep(0,numsims1))
CS_post <- as.vector(rep(0,numsims1))
p_R_pre <- matrix(data = 0, nrow = J, ncol = numsims1)
p_R_post <- matrix(data = 0, nrow = J, ncol = numsims1)
p_W_pre <- matrix(data = 0, nrow = J, ncol = numsims1)
p_W_post <- matrix(data = 0, nrow = J, ncol = numsims1)
shares_pre <- matrix(data = 0, nrow = J, ncol = numsims1)
shares_post <- matrix(data = 0, nrow = J, ncol = numsims1)
convflag_up_pre <- matrix(data = 0, nrow = J, ncol = numsims1)
convflag_up_post <- matrix(data = 0, nrow = J, ncol = numsims1)
convflag_down_pre <- matrix(data = 0, nrow = J, ncol = numsims1)
convflag_down_post <- matrix(data = 0, nrow = J, ncol = numsims1)
profits_pre <- as.vector(rep(0,numsims1))
profits_post <- as.vector(rep(0,numsims1))

for (num in 1:numsims1) {
  print("Scenario")
  print(num)
  print("Pre")
  lambda <- x_vals1[num]
  sigma_val <- sigma_values[num]
  
  ## Pre-merger
  error <- rep(1,J)
  tol <- 0.003
  eps <- rep(0, J)
  
  p_W0 <- c_W_vec*5.0
  p_R0 <- as.numeric((c_R_vec+c_W_vec)*5.5)
  iter <- 1
  
  while (max(error) > tol) {
    for (x in 1:J) {
      
      w_start <- p_W0[x] - .001   # add epsilon
      
      checktest <- optimize(f = Barg_NP_seq3_vert,
                            product_max = x, p_W = p_W0,
                            own_down = own_down_pre, own_up = own_up_pre, 
                            alpha= alpha, delta = delta, 
                            c_W = c_W_vec, c_R = c_R_vec, lambda = lambda, 
                            p_R0 = p_R0, sigma = sigma_val, showAll = FALSE,
                            lower = 0, upper = 1.05)
      
      pj_test <- checktest$minimum
      
      error[x] <- abs(pj_test - p_W0[x])
      print(max(error))
      
      p_W0[x] <- pj_test

      if (convflag_up_pre[x] > 0){
        eps[x] <- .0001
      }
      
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
  
  profits_pre[num] <- (p_R1[1] - p_W1[1] - c_R_vec[1]) * shares_pre[,num][1] +
    (p_R1[2] - p_W1[2] - c_R_vec[2]) * shares_pre[,num][1] + 
    (p_W1[1] - c_W_vec[1]) * shares_pre[,num][1] + 
    (p_W1[3] - c_W_vec[3]) * shares_pre[,num][3]
  
  ## Post-merger
  print("Post")
  error <- rep(1,J)
  tol <- 0.005
  eps <- rep(0, J)
  
  p_W0 <- c_W_vec*5
  p_R0 <- as.numeric((c_R_vec+c_W_vec)*6.5)
  iter <- 1
  
  while (max(error) > tol) {
    for (x in 1:J) {
      
      w_start <- p_W0[x] - .01
      
      checktest <- optimize(f = Barg_NP_seq3_vert,
                            product_max = x, p_W = p_W0,
                            own_down = own_down_post, own_up = own_up_post, 
                            alpha= alpha, delta = delta, 
                            c_W = c_W_vec, c_R = c_R_vec, lambda = lambda, 
                            p_R0 = p_R0, sigma = sigma_val, showAll = FALSE,
                            lower = 0, upper = 1.05)
      
      pj_test <- checktest$minimum
      
      
      error[x] <- abs(pj_test - p_W0[x])
      print(error)
      
      p_W0[x] <- pj_test
      
      if (convflag_up_post[x] > 0){
        eps[x] <- -.01
      }
      
      
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
  
  profits_post[num] <- (p_R1_post[1] - p_W1_post[1] - c_R_vec[1]) * shares_post[,num][1] +
    (p_R1_post[2] - p_W1_post[2] - c_R_vec[2]) * shares_post[,num][1] + 
    (p_W1_post[1] - c_W_vec[1]) * shares_post[,num][1] + 
    (p_W1_post[3] - c_W_vec[3]) * shares_post[,num][3]
  
}



#### 5.3 Simultaneous Model (code copied from simulations-vertical-mergers.Rmd)

CS_pre_C <- as.vector(rep(0,numsims2))
CS_post_C <- as.vector(rep(0,numsims2))
profits_pre_C <- as.vector(rep(0,numsims2))
profits_post_C <- as.vector(rep(0,numsims2))

p_R_pre_C <- matrix(data = 0, nrow = J, ncol = numsims2)
p_R_post_C <- matrix(data = 0, nrow = J, ncol = numsims2)
p_W_pre_C <- matrix(data = 0, nrow = J, ncol = numsims2)
p_W_post_C <- matrix(data = 0, nrow = J, ncol = numsims2)
shares_pre_C <- matrix(data = 0, nrow = J, ncol = numsims2)
shares_post_C <- matrix(data = 0, nrow = J, ncol = numsims2)

for (num in 1:numsims2) {
  
  lambda <- x_vals2[num]
  
  tol <- .0001
  error <- 1
  
  p_W0 <- matrix(.25, nrow = (R*W), ncol = 1)  # dummy wholesale prices
  p_r_start <- as.numeric((c_R_vec+c_W_vec)*6)
  p_R0 <- p_r_start
  
  while (error > tol) {
    
    out1 <- BBoptim(f = Bert_foc ,par = p_R0, 
                    own_down = own_down_pre, alpha= alpha, 
                    delta = delta, cost = c_R_vec,
                    p_W = p_W0, sumFOC = TRUE)
    
    p_R1 <- out1$par
    
    out2 <- BBoptim(par = as.numeric(p_W0), fn = Barg_foc, 
                    own_down = own_down_pre, own_up = own_up_pre, 
                    alpha= alpha, delta = delta, 
                    c_W = c_W_vec, c_R = c_R_vec, lambda = lambda, 
                    p_R = p_R1)
    
    p_W1 <- out2$par
    
    error <- max(abs(c(p_W1-p_W0,p_R1-p_R0)))
    print(error)
    
    p_W0 <- p_W1
    p_R0 <- p_R1
  }
  
  print(p_R1)
  print(p_W1)
  
  shares_pre1 <- (exp(delta - alpha*p_R1))/(1+sum(exp(delta - alpha*p_R1)))
  as.numeric(shares_pre1)
  sum(shares_pre1)
  
  
  tol <- .0001
  error <- 1
  
  p_W0_post <- matrix(.25, nrow = (R*W), ncol = 1)  # dummy wholesale prices
  p_R0_post <- p_r_start
  
  while (error > tol) {
    
    out1 <- BBoptim(f = Bert_foc_vert, par = p_R0_post, 
                    own_down = own_down_post, own_up = own_up_post,
                    alpha = alpha, delta = delta, c_R = c_R_vec, 
                    p_W = p_W0_post, c_W = c_W_vec, sumFOC = TRUE)
    
    p_R1_post <- out1$par
    
    out2 <- BBoptim(par = as.numeric(p_W0_post), fn = Barg_foc_vert, 
                    own_down = own_down_post, own_up = own_up_post, 
                    alpha= alpha, delta = delta, 
                    c_W = c_W_vec, c_R = c_R_vec, lambda = lambda, 
                    p_R = p_R1_post)
    
    p_W1_post <- out2$par
    
    error <- max(abs(c(p_W1_post-p_W0_post,p_R1_post-p_R0_post)))
    print(error)
    
    p_W0_post <- p_W1_post
    p_R0_post <- p_R1_post
  }
  
  p_R1
  p_W1
  p_R1_post
  p_W1_post
  
  # retail price change with vertical integration
  (p_R1_post-p_R1)/p_R1
  
  shares_post1 <- (exp(delta - alpha*p_R1_post))/(1+sum(exp(delta - alpha*p_R1_post)))
  
  CS_pre_C[num] <- (1/alpha) * log(1 + sum(exp(delta - alpha * p_R1)) )
  CS_post_C[num] <- (1/alpha) * log(1 + sum(exp(delta - alpha * p_R1_post)) )
  shares_pre_C[,num] <- shares_pre1
  shares_post_C[,num] <- shares_post1
  
  
  profits_pre_C[num] <- (p_R1[1] - p_W1[1] - c_R_vec[1]) * shares_pre1[1] +
    (p_R1[2] - p_W1[2] - c_R_vec[2]) * shares_pre1[1] + 
    (p_W1[1] - c_W_vec[1]) * shares_pre1[1] + 
    (p_W1[3] - c_W_vec[3]) * shares_pre1[3]
  
  profits_post_C[num] <- (p_R1_post[1] - p_W1_post[1] - c_R_vec[1]) * shares_post1[1] +
    (p_R1_post[2] - p_W1_post[2] - c_R_vec[2]) * shares_post1[1] + 
    (p_W1_post[1] - c_W_vec[1]) * shares_post1[1] + 
    (p_W1_post[3] - c_W_vec[3]) * shares_post1[3]
  
  p_R_pre_C[,num] <- p_R1
  p_R_post_C[,num] <- p_R1_post
  p_W_pre_C[,num] <- p_W1
  p_W_post_C[,num] <- p_W1_post
  
}


# check if mergers are profitable
is_profitable <- (profits_post > profits_pre)
is_profitable
is_profitable_C <- (profits_post_C > profits_pre_C)
is_profitable_C



#### Create the figure
df_fig5 <- data.frame(lambda = x_vals1)
df_fig5$cs_pre <- CS_pre
df_fig5$cs_post <- CS_post
df_fig5$model <- rep(c("sequential linear pricing","sequential two part tariff"), each = numsims1/2)

df_fig5_C <- data.frame(lambda = x_vals2)
df_fig5_C$cs_pre <- CS_pre_C
df_fig5_C$cs_post <- CS_post_C
df_fig5_C$model <- "simultaneous"


df_fig5 <- rbind.data.frame(df_fig5,df_fig5_C)

df_fig5$pct_change_CS <- (df_fig5$cs_post-df_fig5$cs_pre)/df_fig5$cs_pre*100

plot5 <- ggplot(data=df_fig5, 
                aes(x=lambda, y=pct_change_CS, group=model, color=model)) +
  geom_point() +
  xlab("Retailer bargaining power") + ylab("Percent CS change") +
  scale_color_brewer(palette="Dark2") +
  theme(legend.position="bottom")
plot5


#### Table showing underlying pre and post CS

df_fig5$cs_pre <- round(df_fig5$cs_pre, 2)
df_fig5$cs_post <- round(df_fig5$cs_post, 2)
df_fig5$pct_change_CS <- round(df_fig5$pct_change_CS, 1)




#### Table showing product-level details of bargaining varying
#### (Don't show two part tariff, since that is the same for every scenario.)

df_tab5_detail <- data.frame(model = rep(c("sequential linear pricing","simultaneous"), each = J*(numsims1/2+numsims2)/2),
                      good = rep((1:J), times = (numsims1/2+numsims2)) )

df_tab5_detail$p_R_pre <- as.vector(  round(cbind(p_R_pre[,1:5],p_R_pre_C), 1) )
df_tab5_detail$p_R_post <- as.vector( round(cbind(p_R_post[,1:5],p_R_post_C), 1) )
df_tab5_detail$p_W_pre <- as.vector(  round(cbind(p_W_pre[,1:5],p_W_pre_C), 1) )
df_tab5_detail$p_W_post <- as.vector( round(cbind(p_W_post[,1:5],p_W_post_C), 1) )
## VI good in output table:
df_tab5_detail$p_W_post[((1:10) * 4 - 3)] <- NA

df_tab5_detail$shares_pre <- as.vector( round(cbind(shares_pre[,1:5],shares_pre_C)*100,1) )
df_tab5_detail$shares_post <- as.vector( round(cbind(shares_post[,1:5],shares_post_C)*100,1) )

df_tab5_detail$pct_change_p_R <- as.vector( 
  round((df_tab5_detail$p_R_post-df_tab5_detail$p_R_pre)/df_tab5_detail$p_R_pre*100, 1) )

df_tab5_detail$lambda <- rep(c(x_vals1[1:5],x_vals2), each = J)



##################################################
#### 3. store vs. brand (GNL)
##################################################

### Use same baseline as before.
alpha  <- 1.2
delta <- c(3.5, 3.5, 3.5, 3.5)
c_R_vec <- c(0.2,0.2,0.2,0.2)
c_W_vec <- c(0.2,0.2,0.2,0.2)


## GNL parameters
K <- R+W
B1 <- 1 * matrix( c(own_down_pre == "R1",
                    own_down_pre == "R2",
                    own_up_pre == "W1",
                    own_up_pre == "W2"),
                  ncol = K, nrow = J)
a1 <- B1 * 0.5   # rows of a should sum to 1
mu1 <- rep(0.80, K) 


#### x_vals will determine nest membership
x_vals <- c(c(-1.0,-0.5,0,0.5,1.0),   
            c(-1.0,-0.5,0,0.5,1.0))
numsims <- length(x_vals)
sigma_values <- rep( c(0,1), each = length(x_vals)/2)

CS_pre <- as.vector(rep(0,numsims))
CS_post <- as.vector(rep(0,numsims))
p_R_pre <- matrix(data = 0, nrow = J, ncol = numsims)
p_R_post <- matrix(data = 0, nrow = J, ncol = numsims)
p_W_pre <- matrix(data = 0, nrow = J, ncol = numsims)
p_W_post <- matrix(data = 0, nrow = J, ncol = numsims)
shares_pre <- matrix(data = 0, nrow = J, ncol = numsims)
shares_post <- matrix(data = 0, nrow = J, ncol = numsims)
convflag_pre <- as.vector(rep(0,numsims))
convflag_post <- as.vector(rep(0,numsims))
profits_pre <- as.vector(rep(0,numsims))
profits_post <- as.vector(rep(0,numsims))

for (num in 1:numsims) {
  print("Scenario")
  print(num)
  print("Pre Merger")

  a1 <- cbind((-x_vals[num]/2+0.5)*B1[,1],(-x_vals[num]/2+0.5)*B1[,2],
              (x_vals[num]/2+0.5)*B1[,3], (x_vals[num]/2+0.5)*B1[,4])
  sigma_val <- sigma_values[num]
  
  ## Pre-merger
  error <- rep(1,J)
  tol <- 0.01
  
  p_W0 <- c_W_vec*4
  p_R0 <- as.numeric((c_R_vec+c_W_vec)*6)
  iter <- 1
  
  while (max(error) > tol) {
    for (x in 1:J) {
      
      outtest <- optimize(f = Barg_NP_seq3_vert_gnl,
                          product_max = x, p_W = p_W0,
                          own_down = own_down_pre, own_up = own_up_pre, 
                          alpha= alpha, delta = delta, 
                          c_W = c_W_vec, c_R = c_R_vec, lambda = lambda, 
                          p_R0 = p_R0, sigma = sigma_val, 
                          a_jk = a1, B=B1, mu=mu1, showAll = FALSE,
                          lower = 0, upper = 2.0)
      
      pj_test <- outtest$minimum
      
      error[x] <- abs(pj_test - p_W0[x])
      
      p_W0[x] <- pj_test
      
      # recover p_R at these p_W and update r_R0
      outtest_r <- BBoptim(f = Bert_foc_vert_gnl, par = p_R0, 
                           own_down = own_down_pre, own_up = own_up_pre,
                           alpha= alpha, 
                           delta = delta, c_R = c_R_vec,
                           p_W = p_W0, c_W = c_W_vec, 
                           a_jk = a1, B=B1, mu=mu1, sumFOC = TRUE,
                           control = list(trace=FALSE),
                           quiet = TRUE)
      
      p_R0 <- outtest_r$par
      
      iter <- iter + 1
      
    }
  }
  
  p_W1 <- p_W0
  p_R1 <- p_R0
  
  CS_pre[num] <- (1/alpha) * share_calc_gnl(p=p_R1,delta=delta,alpha=alpha,
                                            a_jk=a1,B=B1,mu=mu1, 
                                            returnLogsum=TRUE)
  p_R_pre[,num] <- p_R1
  p_W_pre[,num] <- p_W1
  shares_pre[,num] <- (exp(delta - alpha*p_R1))/(1+sum(exp(delta - alpha*p_R1)))
  
  convflag_pre[num] <- outtest_r$convergence
  
  profits_pre[num] <- (p_R1[1] - p_W1[1] - c_R_vec[1]) * shares_pre[,num][1] +
    (p_R1[2] - p_W1[2] - c_R_vec[2]) * shares_pre[,num][1] + 
    (p_W1[1] - c_W_vec[1]) * shares_pre[,num][1] + 
    (p_W1[3] - c_W_vec[3]) * shares_pre[,num][3]
  
  ## Post-merger
  print("Scenario")
  print(num)
  print("Post Merger")
  
  error <- rep(1,J)
  tol <- 0.01
  
  p_W0 <- p_W1   # use pre-merger as starting values
  p_R0 <- p_R1
  
  iter <- 1
  
  ## Not required to loop over integrated good.
  error[1] <- 0
  
  while (max(error) > tol) {
    for (x in 2:J) {
      
      outtest <- optimize(f = Barg_NP_seq3_vert_gnl,
                          product_max = x, p_W = p_W0,
                          own_down = own_down_post, own_up = own_up_post, 
                          alpha= alpha, delta = delta, 
                          c_W = c_W_vec, c_R = c_R_vec, lambda = lambda, 
                          p_R0 = p_R0, sigma = sigma_val, 
                          a_jk = a1, B=B1, mu=mu1, showAll = FALSE,
                          lower = 0, upper = 2.0)
      
      pj_test <- outtest$minimum
      
      error[x] <- abs(pj_test - p_W0[x])
      
      p_W0[x] <- pj_test
      
      # recover p_R at these p_W and update r_R0
      outtest_r <- BBoptim(f = Bert_foc_vert_gnl, par = p_R0, 
                           own_down = own_down_post, own_up = own_up_post,
                           alpha= alpha, 
                           delta = delta, c_R = c_R_vec,
                           p_W = p_W0, c_W = c_W_vec, a_jk = a1, B=B1, mu=mu1,
                           sumFOC = TRUE,
                           control = list(trace=FALSE),
                           quiet = TRUE)
      
      p_R0 <- outtest_r$par
      
      iter <- iter + 1
      
    }
  }
  
  ## save values from the merger sim
  p_W1_post <- p_W0
  p_R1_post <- p_R0
  
  CS_post[num] <- (1/alpha) * share_calc_gnl(p=p_R1_post,delta=delta,alpha=alpha,
                                             a_jk=a1,B=B1,mu=mu1, 
                                             returnLogsum=TRUE)
  p_R_post[,num] <- p_R1_post
  p_W_post[,num] <- p_W1_post
  shares_post[,num] <- (exp(delta - alpha*p_R1_post))/(1+sum(exp(delta - alpha*p_R1_post)))
  
  convflag_post[num] <- outtest_r$convergence
  
  profits_post[num] <- (p_R1_post[1] - p_W1_post[1] - c_R_vec[1]) * shares_post[,num][1] +
    (p_R1_post[2] - p_W1_post[2] - c_R_vec[2]) * shares_post[,num][1] + 
    (p_W1_post[1] - c_W_vec[1]) * shares_post[,num][1] + 
    (p_W1_post[3] - c_W_vec[3]) * shares_post[,num][3]
  
}

## check if mergers are profitable
is_profitable <- (profits_post > profits_pre)
is_profitable


## Create the figure

df_fig3 <- data.frame(nest = x_vals )

df_fig3$cs_pre <- CS_pre
df_fig3$cs_post <- CS_post
df_fig3$pct_change_CS <- (df_fig3$cs_post-df_fig3$cs_pre)/df_fig3$cs_pre*100

df_fig3$sigma <- rep(c("Linear pricing","Two-part tariff"), each = numsims/2)


plot3 <- ggplot(data=df_fig3, 
                aes(x=nest, y=pct_change_CS, group=sigma, color=sigma)) +
  geom_point() +
  xlab("Nest Structure") + ylab("Percent CS change") +
  scale_color_brewer(palette="Dark2") +
  theme(legend.position="bottom")
plot3

