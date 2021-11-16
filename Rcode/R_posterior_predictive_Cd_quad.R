
falltime_quad = read.csv("falltime_quad.csv", header = TRUE)
falltime_linear = read.csv("falltime_linear.csv", header = TRUE)

# data
linear <-c(1.12180,1.10990,1.09560,1.08020,1.08140,1.08850,1.12060,1.11350,1.12530,1.07420,1.12420,1.12060,1.10040,1.10990,1.09920,1.13490,1.10510,1.09800,1.10870,1.11470,1.11230)
quad <- c(0.925500,0.906100,0.883100,0.858600,0.860500,0.871800,0.923500,0.911900,0.931300,0.849200,0.929400,0.923500,0.890800,0.906100,0.888900,0.947000,0.898400,0.886900,0.904200,0.913800,0.910000)

#d <- linear
d <- quad


library(rstan)

d = list(N= length(d), y = d)
fit1 <- stan(file = "normal.stan", data = d)

fit1

traceplot(fit1, inc_warmup = TRUE)
plot(fit1, pars = c("mu"))
plot(fit1, pars = c("sigma"))


res <- extract(fit1)
y_pred <- res$y_pred
hist(y_pred, breaks = 30)

library(bayesplot)
mcmc_hist(fit1, "y_pred")
mcmc_areas(fit1, "y_pred", prob = 0.95) +
  ggtitle("Posterior predictive distribution of quadratic model Cd", "with median and 95% interval")




