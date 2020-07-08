library(dplyr)
library(tidyverse)
library(ggplot2)
library(MCMCpack)
library(wordcloud2)
wine_data <- read.csv('C:/Users/Rohit/Desktop/winemag-data-130k-v2.csv')
wine_data_Africa <- wine_data %>% filter(
  country == 'South Africa' &
  price  == 15 &
  variety == 'Sauvignon Blanc'
)%>%select(variety, points)
wine_data_Chile <- wine_data %>% filter(
  country == 'Chile' &
    price  == 15 &
    variety == 'Chardonnay'
)%>%select(variety, points)

wine_data_compare <- rbind(wine_data_Africa, wine_data_Chile)

wine_data_compare$variety <- factor(wine_data_compare$variety)

ggplot(wine_data_compare) + geom_boxplot(aes(variety, points, fill = variety)) + geom_jitter(aes(variety, points, shape = variety))+scale_fill_manual(values=c("blue", "red"))

tapply(wine_data_compare$points, wine_data_compare$variety, mean)


t.test(points ~ variety, data=wine_data_compare, var.equal = TRUE)

compare_2_gibbs <- function(y, ind, mu0 = 90, tau0 = 1/100, del0 = 0, gamma0 = 1/100,
                            a0 = 1, b0 = 50, maxiter = 10000)
{
  y1 <- y[ind == 'Sauvignon Blanc']
  y2 <- y[ind == 'Chardonnay']
  n1 <- length(y1)
  n2 <- length(y2)
  ##### starting values
  mu <- (mean(y1) + mean(y2)) / 2
  del <- (mean(y1) - mean(y2)) / 2
  mat_store <- matrix(0, nrow = maxiter, ncol = 3)
  #####
  ##### Gibbs sampler
  an <- a0 + (n1 + n2)/2
  for(s in 1 : maxiter)
  {
    ##update tau
    bn <- b0 + 0.5 * (sum((y1 - mu - del) ^ 2) + sum((y2 - mu + del) ^ 2))
    tau <- rgamma(1, an, bn)
    ##
    ##update mu
    taun <- tau0 + tau * (n1 + n2)
    mun <- (tau0 * mu0 + tau * (sum(y1 - del) + sum(y2 + del))) / taun
    mu <- rnorm(1, mun, sqrt(1/taun))
    ##
    ##update del
    gamman <- tau0 + tau*(n1 + n2)
    deln <- ( del0 * tau0 + tau * (sum(y1 - mu) - sum(y2 - mu))) / gamman
    del<-rnorm(1, deln, sqrt(1/gamman))
    ##
    ## store parameter values
    mat_store[s, ] <- c(mu, del, tau)
  }
  colnames(mat_store) <- c("mu", "del", "tau")
  return(mat_store)
}

fit <- compare_2_gibbs(wine_data_compare$points, wine_data_compare$variety)

raftery.diag(as.mcmc(fit))
apply(fit, 2, mean)
mean(1/sqrt(fit[, 3]))
plot(as.mcmc(fit))

y1_sim <- rnorm(10000, fit[, 1] + fit[, 2], sd = 1/sqrt(fit[, 3]))
y2_sim <- rnorm(10000, fit[, 1] - fit[, 2], sd = 1/sqrt(fit[, 3]))

mean(y1_sim > y2_sim)

ggplot(data.frame(y_sim_diff = y1_sim - y2_sim)) + stat_bin(aes(y_sim_diff))

ggplot(data.frame(y1_sim, y2_sim)) + geom_point(aes(y1_sim, y2_sim), alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0)


wine_data_italy <- wine_data %>% filter(
  country == 'Italy',
  price < 20
)

wine_data_italy <- wine_data_italy%>%select(region_1,points)

wine_data_italy <- wine_data_italy[!(wine_data_italy$region_1 == ""), ]

average_mean <- mean(wine_data_italy$points)

keepers <- names(table(wine_data_italy$region_1))[table(wine_data_italy$region_1) > 3]
wine_data_italy <- wine_data_italy[wine_data_italy$region_1 %in% keepers,]

wine_data_italy <- wine_data_italy%>%filter(
  region_1 != 'Sardinia'
)

wine_data_italy$region_1 <- factor(wine_data_italy$region_1)

wine_data_italy <- wine_data_italy%>%
  group_by(region_1)%>%
  summarise(average_point = mean(points))

above_average <- wine_data_italy%>%filter(
  average_point >= average_mean
)

ggplot(wine_data_italy) + geom_boxplot(aes(x = reorder(region_1, points, median), points,
                               fill = reorder(region_1, points, median)), show.legend=FALSE)

ggplot(wine_data_italy, aes(x = reorder(region_1, region_1, length))) + stat_count()

ggplot(data.frame(size = tapply(wine_data_italy$points, wine_data_italy$region_1, length),
                  mean_score = tapply(wine_data_italy$points, wine_data_italy$region_1, mean)), aes(size, mean_score)) +
  geom_point()

compare_m_gibbs <- function(y, ind, mu0 = 50, tau0 = 1/400,
                            a0 = 1, b0 = 50, alpha0 =1, beta0 = 50, maxiter = 500)
{
  ### weakly informative priors
  a0 <- 1/2 ; b0 <- 50 ## tau_w hyperparameters
  alpha0 <-1/2 ; beta0 <- 50 ## tau_b hyperparameters
  mu0<-50 ; tau0 <- 1/25
  ###
  ### starting values
  m <- nlevels(ind)
  ybar <- theta <- tapply(y, ind, mean)
  tau_w <- mean(1 / tapply(y, ind, var)) ##within group precision
  mu <- mean(theta)
  tau_b <-var(theta) ##between group precision
  n_m <- tapply(y, ind, length)
  alphan <- alpha0 + sum(n_m)/2
  ###
  ### setup MCMC
  theta_mat <- matrix(0, nrow=maxiter, ncol=m)
  mat_store <- matrix(0, nrow=maxiter, ncol=3)
  ###
  ### MCMC algorithm
  for(s in 1:maxiter)
  {
    # sample new values of the thetas
    for(j in 1:m)
    {
      taun <- n_m[j] * tau_w + tau_b
      thetan <- (ybar[j] * n_m[j] * tau_w + mu * tau_b) / taun
      theta[j]<-rnorm(1, thetan, 1/sqrt(taun))
    }
    #sample new value of tau_w
    ss <- 0
    for(j in 1:m){
      ss <- ss + sum((y[ind == j] - theta[j])^2)
    }
    betan <- beta0 + ss/2
    tau_w <- rgamma(1, alphan, betan)
    #sample a new value of mu
    taum <- m * tau_b + tau0
    mum <- (mean(theta) * m * tau_b + mu0 * tau0) / taum
    mu <- rnorm(1, mum, 1/ sqrt(taum))
    # sample a new value of tau_b
    am <- a0 + m/2
    bm <- b0 + sum((theta - mu)^2) / 2
    tau_b <- rgamma(1, am, bm)
    #store results
    theta_mat[s,] <- theta
    mat_store[s, ] <- c(mu, tau_w, tau_b)
  }
  colnames(mat_store) <- c("mu", "tau_w", "tau_b")
  return(list(params = mat_store, theta = theta_mat))
}
fit2 <- compare_m_gibbs(wine_data_italy$points, factor(as.numeric(wine_data_italy$region_1)))

apply(fit2$params, 2, mean)

theta_hat <- apply(fit2$theta, 2, mean)

sum(theta_hat > 86.5491850)

ggplot(data.frame(size = tapply(wine_data_italy$points, wine_data_italy$region_1, length), theta_hat = theta_hat),
       aes(size, theta_hat)) + geom_point()

theta_df <- data.frame(samples = as.numeric(fit2$theta),
                       region_1 = rep(1:ncol(fit2$theta), each = nrow(fit2$theta)))
ggplot(theta_df) + geom_boxplot(aes(x = reorder(region_1, samples, median), samples,
                                    fill = reorder(region_1, samples, median)), show.legend=FALSE)

