
### Autors: Pablo S. Ascandoni & Lukas Malte Kemeter

### Acknowledge:  This code was written for an assignment of the course 
#                 "Advanced Empirical Finance - Topics and Data Science"
#                 using the tools provided by the professor Stefan Voigt.

# This project replicates the main findings of Jobson, J. D., & Korkie, B. (1980) and performs additional analysis

# Data: Kenneth French website --> monthly returns for 10 different industries from July 1926 up to November 2020

### References:
# Jobson, J. D., & Korkie, B. (1980). Estimation for Markowitz efficient portfolios. Journal of the American Statistical Association, 75(371), 544-554.




 ################################  ASSIGNMENT 1  ##############################

 
# BEFORE YOU RUN THIS CODE:
# 1. Restart your current R session. The package MASS seems to interfere with the select command 
# 2. In line 380 (Exercise 8) we simulate 250 times with T = 1.000.000 draws so that will take approx. 15 mins to run. 
#    If you do not like to wait that long, go to exercise 8 and copy out this specific line 
#    As a result, the 4th histogram that uses this simulation data would then obviously not be generated anymore 


rm(list = ls()) # clear workspace environment



library(tidyverse)
library(lubridate)
library(tidyquant)
library(reshape)
library(quadprog)
library(knitr)
library(kableExtra)


data <- read.csv("10_Industry_Portfolios.csv", skip=11)                         # Kenneth French data

### Organise the dataset:

# Create date:
data <- data %>% mutate(date = ymd(parse_date_time(X, "%Y%m"))) %>%             # Create date variable Year-Month-Day
  select(-X) %>%                                                                # Eliminate previous  
  relocate(date)                                                                # Put date column as first

# Create panel shape with industry variable:
data <- data %>% pivot_longer(!date, names_to = "industry", values_to = "weighted_returns")


################################################################################

### 1) Industry Sharpe Ratios:

# Create TxN matrix:
data <- data %>% pivot_wider(names_from = "industry", values_from="weighted_returns") %>%
  select(-date)

#Obtain mu and sigma:
mu <- data %>% colMeans()
sigma <- data %>% cov()

# Create LaTeX Tables:
options(knitr.table.format ="latex")
knitr::kable(t(mu), digits = 2)
knitr::kable(sigma, digits = 2)

# Industry Sharpe-ratios : 
(mu / sqrt(diag(sigma)))
options(knitr.table.format ="latex")
knitr::kable(t(mu / sqrt(diag(sigma))), digits = 2)
# Which industry delivered the highest Sharpe-ratio? : 
(mu / sqrt(diag(sigma))) %>% which.max()


################################################################################

### 2) Create function that computes efficient frontier:

compute_efficient_frontier <- function(sigma_hat, mu_hat, sigma_true, mu_true){
  # Inputs: NxN sigma matrix and 1xN mu vector of returns
  # Output: DataFrame with columns c,mu & sd
  
  ## Minimum Variance Portfolio (MVP) weights, expected return and volatility :
  sigma_inv <- solve(sigma_hat)                                                 # solve() inverts the matrix
  N <- ncol(sigma_hat)
  iota <- rep(1,N)                                                              # vector Nx1 of 1s
  w_mvp = (sigma_inv %*% iota) / sum(sigma_inv %*% iota)                        # MVP weights
  
  ## Efficient Portfolio weights of 2*mu_mvp :
  mu_bar <- 2 * t(w_mvp) %*% mu_hat                                             # Desired expected return
  C <- as.numeric(t(iota) %*% sigma_inv %*% iota)
  D <- as.numeric(t(iota) %*% sigma_inv %*% mu_hat)
  E <- as.numeric(t(mu_hat) %*% sigma_inv %*% mu_hat)
  lambda_tilde <- as.numeric(2 * (mu_bar - D / C) / (E - D ^ 2 / C))
  w_eff <- w_mvp + (lambda_tilde / 2) * (sigma_inv %*% mu_hat - D / C * sigma_inv %*% iota)   # Efficient weights:
  
  ## Mutual Fund Theoremevaluated at true parameters:
  c <- seq(from = -0.1, to = 1.2, by = 0.01)                                    # sequence of numbers
  res <- tibble(c = c, mu = NA, sd = NA)                                        # tibble() creates dataframe with columns c,mu & sd
  for(i in seq_along(c)){                                                       # For loop that fills in tibble with mu's and sd's. seq_along()=range()
    w <- (1-c[i])*w_mvp + c[i]*w_eff                                            # Weights
    res$mu[i] <- t(w) %*% mu_true                                               # Portfolio expected return evaluated at true mu
    res$sd[i] <- sqrt(t(w) %*% sigma_true %*% w)                                # Portfolio volatility evaluated at true sigma
  }
  
  return(res)
}

################################################################################

#### 3) Plot efficient frontier:

res <- compute_efficient_frontier(sigma, mu, sigma, mu)

# Plot the Efficient Frontier:
ggplot(res, aes(x = sd, y = mu)) + 
  geom_point(size=0.7) +                                                        # Plot all sd/mu portfolio combinations 
  geom_point(data = res %>% filter(c %in% c(0,1)), color = "red", size = 2) +   # locate  mvp (c==0) and eff portfolio (c==1)
  geom_point(data = tibble(mu = mu, sd = sqrt(diag(sigma))), 
             aes(y = mu, x  = sd), color = "blue", size = 1) +                  # locate the individual assets 
  labs(x="Standard deviation", y="Expected return") +
  geom_text(aes(x = 11.5, y = 1.72, label = "EFF"), size=3, color = "red") +
  geom_text(aes(x = 4.3, y = 0.9, label = "MVP"), size=3, color = "red") +
  theme(legend.position="none") +
  ylim(0.8, 1.8) +
  xlim(3.5, 12) +
  theme_minimal() +
  theme(
    axis.line        = element_line(colour = "black", size = .5),
    axis.ticks       = element_line(colour = "black", size = .5),
    axis.text        = element_text(size = 10, colour = "black"),
    axis.title       = element_text(size = 10),
  )

ggsave(filename=paste("plot_ex3.png",sep=""))


################################################################################

#### 4) Efficient tangent portfolio:

# Tangent portfolio weights for risk-free rate=0:
rf <- 0
N <- ncol(sigma)
iota <- rep(1,N)  
w_tgc <- (solve(sigma) %*% (mu-rf)) / as.numeric(t(iota) %*% solve(sigma) %*% (mu-rf))
knitr::kable(t(w_tgc), digits = 2)

# Efficient tangent portfolio:
tgc <- t(c(t(w_tgc) %*% mu, sqrt(t(w_tgc) %*% sigma %*% w_tgc))) %>% as_tibble()
colnames(tgc) <- c("mu","sd")

# Plot the Efficient Frontier with tangent line:
ggplot(res, aes(x = sd, y = mu)) + 
  geom_point(size=0.7) +                                                        # Plot all sd/mu portfolio combinations 
  geom_point(data = res %>% filter(c %in% c(0,1)), color = "red", size = 2) +   # locate  mvp (c==0) and eff portfolio (c==1)
  geom_point(data = tibble(mu = mu, sd = sqrt(diag(sigma))), 
             aes(y = mu, x  = sd), color = "blue", size = 1) +                  # locate the individual assets 
  geom_point(data = tgc, color="purple", size=2) +
  geom_abline(slope=(tgc$mu-rf)/tgc$sd , intercept = rf, color="purple") + 
  labs(x="Standard deviation", y="Expected return") +
  geom_text(aes(x = 11.5, y = 1.72, label = "EFF"), size=3, color = "red") +
  geom_text(aes(x = 3.8, y = 1.05, label = "TGC"), size=3, color = "purple") +
  geom_text(aes(x = 4.3, y = 0.9, label = "MVP"), size=3, color = "red") +
  theme(legend.position="none") +
  ylim(0.8, 1.8) +
  xlim(3.5, 12) +
  theme_minimal() +
  theme(
    axis.line        = element_line(colour = "black", size = .5),
    axis.ticks       = element_line(colour = "black", size = .5),
    axis.text        = element_text(size = 10, colour = "black"),
    axis.title       = element_text(size = 10),
  )
ggsave(filename=paste("plot_ex4.png",sep=""))


# Sharpe ratio = slope:
sharp_ratio <- (tgc$mu-rf)/tgc$sd
sprintf("Sharpe ratio tangent portfolio: %s", round(sharp_ratio, 2))

################################################################################

### 5) Function that simulates returns:

library(MASS)     # From now onward select() disabled

simulate_returns <- function(T, mu, sigma, seed){
  # Simulate sample of returns from the real data
  set.seed(seed)
  as.data.frame(mvrnorm(n = T, mu, sigma))
}

################################################################################

#### 6) Simulate sample of returns and obtain frontier:

data_s <- simulate_returns(100, mu, sigma, 2020)

mu_s <- data_s %>% colMeans()
sigma_s <- data_s %>% cov()
res_s <- compute_efficient_frontier(sigma_s, mu_s, sigma, mu)

# Plot both Efficient Frontiers together:
ggplot(data = res, aes(x = sd, y = mu)) + 
  geom_point(data = res, color="blue", size=0.7) +                                                    # Plot all sd/mu portfolio combinations 
  geom_point(data = res %>% filter(c %in% c(0,1)), color = "red", size = 2) +   # locate  mvp (c==0) and eff portfolio (c==1)
  geom_point(data = res_s, size=0.7) +                                                    # Plot all sd/mu portfolio combinations 
  geom_point(data = res_s %>% filter(c %in% c(0,1)), color = "red", size = 2) + # locate  mvp (c==0) and eff portfolio (c==1)
  geom_abline(slope=(tgc$mu-rf)/tgc$sd , intercept = rf, color="purple") +
  geom_point(data = tgc, color="purple", size=2) +
  labs(x="Standard deviation", y="Expected return") +
  theme_minimal() +
  theme(
    axis.line        = element_line(colour = "black", size = .5),
    axis.ticks       = element_line(colour = "black", size = .5),
    axis.text        = element_text(size = 10, colour = "black"),
    axis.title       = element_text(size = 10),
  )
ggsave(filename=paste("plot_ex6.png",sep=""))


#### 7) Plot N efficient frontiers:                     # (the code takes time to plot the graph)

# Create function that plots the efficient frontiers of T samples and obtains SR:

compute_efficient_frontier_T <- function(res, mu, sigma, sample_size, iterations){
  
  rf <- 0
  plot_i <- ggplot(res, aes(x = sd, y = mu))
  sharp_ratios <- rep(NA, iterations)

  for (i in 1:iterations) {
    data_i <- simulate_returns(sample_size, mu, sigma, i)
    mu_i <- data_i %>% colMeans()
    sigma_i <- data_i %>% cov()
    res_i <- compute_efficient_frontier(sigma_i, mu_i, sigma, mu)
  
    # Tangent sharp ratio:
    N <- ncol(sigma_i)
    iota <- rep(1,N)  
    w_tgc_i <- (solve(sigma_i) %*% (mu_i-rf)) / as.numeric(t(iota) %*% solve(sigma_i) %*% (mu_i-rf))
    sharp_ratio_i <- (t(w_tgc_i) %*% mu)/ sqrt(t(w_tgc_i) %*% sigma %*% w_tgc_i)
    sharp_ratios[i] <- sharp_ratio_i
  
    plot_i <- plot_i + geom_point(data=res_i, size=0.5, alpha=0.05, color="black")                                                                
  }
  
  #Make sharp_ratios a data frame:
  sharp_ratios <- sharp_ratios %>% as.data.frame()
  
  # Tangent sharp ratio:
  N <- ncol(sigma)
  iota <- rep(1,N)  
  w_tgc <- (solve(sigma) %*% (mu-rf)) / as.numeric(t(iota) %*% solve(sigma) %*% (mu-rf))
  sharp_ratio <- (t(w_tgc) %*% mu)/ sqrt(t(w_tgc) %*% sigma %*% w_tgc)
  
  plot_i <- plot_i + geom_point(data=res, size=0.7, color="blue") +                                                                
    geom_point(data = res %>% filter(c %in% c(0,1)), color = "red", size = 2) +
    geom_point(data = tgc, color="purple", size=2) +
    geom_abline(slope=sharp_ratio , intercept = rf, color="purple") +
    geom_point(data = tgc, color="purple", size=2) +
    labs(x="Standard deviation", y="Expected return") +
    ylim(0.5, 1.8) +
    xlim(3.5, 12) +
    theme_minimal() +
    theme(
      axis.line        = element_line(colour = "black", size = .5),
      axis.ticks       = element_line(colour = "black", size = .5),
      axis.text        = element_text(size = 10, colour = "black"),
      axis.title       = element_text(size = 10),
    )
  
  print(plot_i)
  ggsave(plot_i, filename=paste("plot_ex7_", sample_size, ".png", sep=""))
  return(sharp_ratios)
}

sharp_ratio_100 <- compute_efficient_frontier_T(res, mu, sigma, 100, 250)
ggsave(filename=paste("plot_ex7.png",sep=""))

# Run it for different T:
sharp_ratio_250 <- compute_efficient_frontier_T(res, mu, sigma, 250, 250)
sharp_ratio_500 <- compute_efficient_frontier_T(res, mu, sigma, 500, 250)
sharp_ratio_1000 <- compute_efficient_frontier_T(res, mu, sigma, 1000, 250)
# sharp_ratio_10000 <- compute_efficient_frontier_T(res, mu, sigma, 10000, 250)
# sharp_ratio_1000000 <- compute_efficient_frontier_T(res, mu, sigma, 1000000, 250)

# Create data frame with all sharpe ratios for each T:
sharp_ratio_all <- data.frame("T_100"=sharp_ratio_100,"T_250"=sharp_ratio_250,"T_500"=sharp_ratio_500,"T_1000"=sharp_ratio_1000)
colnames(sharp_ratio_all) <- c("T_100","T_250","T_500","T_1000")

###############################################################################

### 8) Histogram of Sharpe Ratios:

T=100
hist1 <- ggplot(data=sharp_ratio_100, aes(sharp_ratio_100$.)) + 
  labs(x="Sharpe Ratios", y="# observations") +
  xlim(0, .25) +
  ylim(0, 110) +
  geom_histogram(bins=30, colour="black", fill="black", alpha=0.25)+
  geom_vline(xintercept=sharp_ratio, colour="red") +
  theme_minimal() +
  theme(
    axis.line        = element_line(colour = "black", size = .5),
    axis.ticks       = element_line(colour = "black", size = .5),
    axis.text        = element_text(size = 10, colour = "black"),
    axis.title       = element_text(size = 10),
  )
ggsave(hist1, filename=paste("plot_ex9_",T,".png",sep=""))


T=250
hist2 <- ggplot(data=sharp_ratio_250, aes(sharp_ratio_250$.)) + 
  labs(x="Sharpe Ratios", y="# observations") +
  xlim(0, .25) +
  ylim(0, 110) +
  geom_histogram(bins=30, colour="black", fill="black", alpha=0.25)+
  geom_vline(xintercept=sharp_ratio, colour="red") +
  theme_minimal() +
  theme(
    axis.line        = element_line(colour = "black", size = .5),
    axis.ticks       = element_line(colour = "black", size = .5),
    axis.text        = element_text(size = 10, colour = "black"),
    axis.title       = element_text(size = 10),
  )
ggsave(hist2, filename=paste("plot_ex9_",T,".png",sep=""))

T=500
hist3 <- ggplot(data=sharp_ratio_500, aes(sharp_ratio_500$.)) + 
  labs(x="Sharpe Ratios", y="# observations") +  
  xlim(0, .25) +
  ylim(0, 110) +
  geom_histogram(bins=30, colour="black", fill="black", alpha=0.25)+
  geom_vline(xintercept=sharp_ratio, colour="red") +
  theme_minimal() +
  theme(
    axis.line        = element_line(colour = "black", size = .5),
    axis.ticks       = element_line(colour = "black", size = .5),
    axis.text        = element_text(size = 10, colour = "black"),
    axis.title       = element_text(size = 10),
  )
ggsave(hist3, filename=paste("plot_ex9_",T,".png",sep=""))

T=1000
hist4 <- ggplot(data=sharp_ratio_1000, aes(sharp_ratio_1000$.)) + 
  labs(x="Sharpe Ratios", y="# observations") +
  xlim(0, .25) +
  ylim(0, 110) +
  geom_histogram(bins=30, colour="black", fill="black", alpha=0.25)+
  geom_vline(xintercept=sharp_ratio, colour="red") +
  theme_minimal() +
  theme(
    axis.line        = element_line(colour = "black", size = .5),
    axis.ticks       = element_line(colour = "black", size = .5),
    axis.text        = element_text(size = 10, colour = "black"),
    axis.title       = element_text(size = 10),
  )
ggsave(hist4, filename=paste("plot_ex9_",T,".png",sep=""))


# Compare three density plots

sharp_ratio_all <- melt(sharp_ratio_all)

ggplot(sharp_ratio_all, aes(x=value, fill=variable)) + geom_density(alpha=0.25) +
  geom_vline(xintercept=sharp_ratio, colour="red") + 
  labs(x="Sharpe Ratios", y="Density %") +
  theme_minimal() +
  theme(
    axis.line        = element_line(colour = "black", size = .5),
    axis.ticks       = element_line(colour = "black", size = .5),
    axis.text        = element_text(size = 10, colour = "black"),
    axis.title       = element_text(size = 10),
  ) 
ggsave(filename=paste("plot_ex9_densities.png",sep=""))

################################################################################

### 9) Implementation of alternative strategies (inspired by: Niels Eriksen and Rasmus Juel)

# Implementing no short sale T = 100
c <- seq(1, 250, 1)
no_short_sell_portfolios <- tibble(c,  mu = NA, sd = NA, sharpe = NA)
T = 100

for (i in c){    
  return_sample <- simulate_returns(T, mu, sigma, i) 
  sigma_sample <- cov(return_sample)
  A <- t(cbind(1, diag(N)))
  
  no_short_sell_eff <- solve.QP(Dmat = sigma_sample,
                                dvec = rep(0, N), 
                                Amat = t(A), 
                                bvec = c(1, rep(0, N)), 
                                meq = 1)
  
  w <- no_short_sell_eff$solution
  
  # Calculating return, sd and sharpe
  no_short_sell_portfolios$mu[i] <- t(w) %*% mu
  no_short_sell_portfolios$sd[i] <- sqrt(t(w) %*% sigma %*% w)
  no_short_sell_portfolios$sharpe[i] <- no_short_sell_portfolios$mu[i]/no_short_sell_portfolios$sd[i]
}

no_short_sell_portfolios %>% ggplot(aes(x= sharpe)) +
  geom_histogram(bins=30, colour="black", fill="black", alpha=0.25)+
  geom_vline(xintercept=sharp_ratio, colour="red") +
  xlim(0, .25) +
  ylim(0, 200) +
  labs(x="Sharpe Ratios No Short Selling", y="# observations") +
  theme_minimal() +
  theme(
    axis.line        = element_line(colour = "black", size = .5),
    axis.ticks       = element_line(colour = "black", size = .5),
    axis.text        = element_text(size = 10, colour = "black"),
    axis.title       = element_text(size = 10),
  )
ggsave(filename=paste("ex9_no_short_100.png",sep=""))

# Implement naive portfolio
w_naive <- rep(1/N, N)

# Calculating return, sd and sharpe
mu_naive <- t(w_naive) %*% mu
sigma_naive <- sqrt(t(w_naive) %*% sigma %*% w_naive)
naive_sharpe <- mu_naive/sigma_naive

naive_vec <- as.data.frame(c(mu_naive,sigma_naive,naive_sharpe))
options(knitr.table.format ="latex")
knitr::kable(t(naive_vec), digits = 2)
