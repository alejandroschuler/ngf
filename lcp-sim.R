library(tidyverse)
library(broom)
library(zeallot)
library(magrittr)
library(survival)

# DGP for y(t)
gen_y = function(t, x, t_treat, t_surv, t_cens) {
  if (t_cens <= t & t_cens < t_surv) {
    return(NA)
  } else if (t_surv <= t & t_surv <= t_cens) {
    return(0)
  } else {
    y_t = 10 #+ exp(rnorm(1,x[[1]])) # normal cost depends on this covariate
    if (abs(t-t_surv) < 12) { # cost spikes in the last 6 months
      y_t = 10 + y_t
      if (!is.na(t_treat) & (t_treat<=t)) {
        y_t = y_t - 10 # treatment negates the cost spike
      }
    }
    return(y_t)
  }
}

# DGP for T, C, Z
gen_data = function(id, p, t_max = 10*12) {
  # x = rnorm(p) %>%
  #   set_names(str_c("x", 1:p))
  t_treat = round(runif(1, 1, t_max))
  # t_treat = round(abs(rnorm(1,60,30)))
  # t_treat = if_else(rbinom(1, 1, 0.5)==1, 1, t_max+1)

  t_surv = round(runif(1, 1, t_max)) # + max(12*x[[1]], 0))
  # t_surv = t_max - 12
  # t_cens = round(runif(1, 1, t_max))
  t_cens = t_max # C = inf

  t_obs = min(t_surv, t_cens)
  t_treat = if_else(t_treat < t_obs, t_treat, NA_real_)
  d = t_surv <= t_cens

  y = 1:t_max %>% map_dbl(gen_y, x, t_treat, t_surv, t_cens)

  list(patient = enframe(x) %>%
         spread(name, value) %>%
         mutate(id=id, t_treat=t_treat, t_obs=t_obs, d),
       monthly = tibble(t=1:t_max, y=y, id=id))
}

# get the monthly data out of the full data structure
data_t = function(tt) {
  patient %>%
    mutate(w = !is.na(t_treat) & (t_treat <= tt),
           Tk = pmin(t_obs, tt), # t_obs = t_surv when d==1
           dk = (tt < t_obs) | d, # always observed if t < t_obs and if death was observed
           weight = dk/G_hat(Tk)) %>%
    filter(!weight==0) %>%
    left_join(filter(monthly, t==tt), by=c("id"))
}

# do the monthly regression
regression_t = function(tt) {
  data_t(tt) %>%
    # lm(y ~ w+x1+x2+x3+x4+x5, weights=weight, data=.) %>%
    lm(y ~ w, weights=weight, data=.) %>%
    tidy()
}
#
# patient %>%
#   left_join(monthly, by="id") %>%
#   filter(id<=20) %$%
#   qplot(t,y, group=factor(id), color=(!is.na(t_treat) & (t_treat < t_obs)), geom="line")
#
# 1:120 %>%
#   map_dfr(function(t) {
#     data_t(t) %>%
#       group_by(alive = t < t_obs, w) %>%
#       summarize(n=n(), mean_y=mean(y))  %>%
#       mutate(t=t)
#     }) %>%
#   ggplot(aes(t,n, color=alive)) + geom_point(aes(size=mean_y, shape=w), alpha=0.5)
#
# 1:120 %>%
#   map_dfr(function(t) {
#     data_t(t) %>%
#       group_by(alive = t < t_obs, w) %>%
#       summarize(N=n(), mean_y=mean(y)) %>%
#       mutate(t=t)
#   }) %>%
#   ggplot(aes(t, mean_y, color=alive)) + geom_point(aes(size=N, shape=w), alpha=0.5)
#
# 1:120 %>%
#   map_dfr(function(t) {
#     data_t(t) %>%
#       group_by(w) %>%
#       summarize(N=n(), mean_y=mean(y)) %>%
#       mutate(t=t)
#   }) %>%
#   ggplot(aes(t, mean_y)) + geom_line(aes(color=w))

# generate the data
1:10000 %>%
  map(gen_data, p=5) %>%
  transpose() %>%
  map(bind_rows) %->%
  c(patient, monthly)

# get the KM estimator for censoring
km_cens = survfit(Surv(t_obs, 1-d)~1, patient)
G_hat = stepfun(km_cens$time, c(1, km_cens$surv)) # P(not censored at time t)

data_t = function(tt) {
  patient %>%
    mutate(w = !is.na(t_treat) & (t_treat <= tt),
           Tk = pmin(t_obs, tt), # t_obs = t_surv when d==1
           dk = (tt < t_obs) | d, # always observed if t < t_obs and if death was observed
           G = G_hat(Tk),
           # weight = 1) %>%
           weight = dk/(G)) %>%
    filter(!weight==0) %>%
    left_join(filter(monthly, t==tt), by=c("id"))
}

# do the monthly regression
regression_t = function(tt) {
  data_t(tt) %>%
    # lm(y ~ w+x1+x2+x3+x4+x5, weights=weight, data=.) %>%
    lm(y ~ w, weights=weight, data=.) %>%
    tidy()
}

# do all the monthly regressions
result = 1:120 %>%
  map(function(t) regression_t(t) %>% mutate(t=t)) %>%
  bind_rows() %>%
  filter(term=="wTRUE")

result %$% qplot(t, estimate) # beta over time
result %$% sum(estimate) # sum of betas
