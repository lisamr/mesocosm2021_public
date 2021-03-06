#code needed for defining beta-binomial as a custom family
# copied from vignette, written by Paul Burkner: https://cran.r-project.org/web/packages/brms/vignettes/brms_customfamilies.html

library(brms)
beta_binomial2 <- custom_family(
  "beta_binomial2", dpars = c("mu", "phi"),
  links = c("logit", "log"), lb = c(NA, 0),
  type = "int", vars = "vint1[n]"
)
stan_funs <- "
  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
  }
"
stanvars <- stanvar(scode = stan_funs, block = "functions")

#You would fit your model, then expose stan functions.
# fit2 <- brm(
#   incidence | vint(size) ~ period + (1|herd), data = cbpp, 
#   family = beta_binomial2, stanvars = stanvars
# )
# expose_functions(fit2, vectorize = TRUE)

#post processing
log_lik_beta_binomial2 <- function(i, prep) {
  mu <- prep$dpars$mu[, i]
  phi <- prep$dpars$phi
  trials <- prep$data$vint1[i]
  y <- prep$data$Y[i]
  beta_binomial2_lpmf(y, mu, phi, trials)
}

#Loo and K-fold work now too.

posterior_predict_beta_binomial2 <- function(i, prep, ...) {
  mu <- prep$dpars$mu[, i]
  phi <- prep$dpars$phi
  trials <- prep$data$vint1[i]
  beta_binomial2_rng(mu, phi, trials)
}
# pp_check(fit2, nsamples = 50)


posterior_epred_beta_binomial2 <- function(prep) {
  mu <- prep$dpars$mu
  trials <- prep$data$vint1
  trials <- matrix(trials, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  mu * trials
}
#conditional_effects(fit2, conditions = data.frame(size = 1))

