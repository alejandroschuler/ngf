library(tidyverse)
library(broom)
library(furrr)
library(zeallot)
library(tidymodels)

as_factor.logical = function(lgl, ...) base::factor(lgl, levels=c("TRUE", "FALSE"), ...)
as_factor.default = function(x, ...) as.factor(x, ...)
fct_to_lgl = function(x) as.logical(2-as.numeric(x)) # amazingly, 3x as fast as as.logical()

# To do:
# add multiple confounders to sim
# functionalize naming

tau = 100 # months
n = 1000
p = 1

natural = function(data) if_else(data[,"D"]==1, data[,"A"], pmax(data[,"A"], rbinom(nrow(data),1,0.05))) # 5% chance of being treated each interval, treatment carries over
treat_all = function(data) rep(1,nrow(data))
treat_none = function(data) rep(0,nrow(data))

rmultinom = function(...) { # represented as integers in the matrix data
  data = stats::rmultinom(...)
  which(data==1, arr.ind = T)[,1]
}

# true_init = function(treat_plan, n) {
#   data = matrix(0, nrow=n, ncol=5)
#   colnames(data) = c("L1", "L2", "A", "D", "Y")
#   data[,"L1"] = rep(0, n)
#   data[,"L2"] = rmultinom(n, 1, c(1, 1, 2))
#   data[,"A"] = rep(0, n) # everyone starts not treated
#   data[,"D"] =
# }

true_next = function(data, treat_plan) {
  data_next = data
  data_next[,"L1"] = data[,"L1"] + rnorm(nrow(data), 0,1) + rep(1,nrow(data)) # disease progression
  data_next[,"L2"] = data[,"L2"] # category of some kind
  data_next[,"A"] = treat_plan(data)
  data_next[,"D"] = pmax(data[,"D"], rbinom(nrow(data),1, gtools::inv.logit(data[,"L1"]-10))) # p(death) increases with disease progression
  data_next[,"Y"] = (1-data[,"D"])*(1.1^(data[,"L1"] - 2*data[,"A"] + abs(rnorm(nrow(data), 0,1)))) # age increases cost
  data_next
}

matrix_to_df = function(data_matrix, name) {
  as_tibble(data_matrix) %>%
    set_names(str_c("t",1:ncol(.),sep="_")) %>%
    mutate(i=1:nrow(.)) %>%
    gather(t, !!name, -i) %>%
    separate(t, into=c("trash","t"), sep="_") %>%
    mutate(t=as.integer(t)) %>%
    select(-trash)
}

sim_data = function(sim_next, treat_plan, n, tau, p) {
  data_array = array(0, c(n, tau, p+3)) # p Ls + A+ D+ Y
  dimnames(data_array)[[3]] = c(str_c("L",1:p), "A", "D", "Y")
  for (t in 1:(tau-1)){
    data_array[,t+1,] = sim_next(data_array[,t,], treat_plan)
  }
  dimnames(data_array)[[3]] %>%
    map(~matrix_to_df(data_array[,,.], .)) %>%
    reduce(inner_join, by=c("t","i"))
}

data = sim_data(true_next, natural, n, tau, p=1)

# references to "D" not good
modeling_data = function(data, x) {
  data %>%
    select(i,t,x=!!sym(x)) %>%
    inner_join(data %>%
                 filter(!D) %>%
                 mutate(t=t+1), # put x(t+1) on the same row as x(t), z1(t), z2(t) etc.
               by=c("t","i")) %>%
    select(-t, -i, -D)
}

model_specs = list(
  "L1" = linear_reg() %>% set_engine("lm"),
  "D" = logistic_reg() %>% set_engine("glm"),
  "Y" = rand_forest() %>% set_engine("ranger")
)

# model_specs = list(
#   "L1" = rand_forest() %>% set_engine("ranger") %>% set_mode("regression"),
#   "D" = logistic_reg() %>% set_engine("glm"),
#   "Y" = rand_forest() %>% set_engine("ranger") %>% set_mode("regression")
# )

honest_sd = function(model, data) { # get an honest estimate of the sd that isn't biased by overfitting
  if (model$spec$engine == "glm" | model$spec$engine == "lm") {
    return(sigma(model$fit))
  } else if (model$spec$engine == "ranger") {
    return(model$fit$prediction.error)
  } else {
    folds = vfold_cv(data, 5)
    folds$splits %>%
      map_dbl(function(fold) {
        training = analysis(fold)
        validation = assessment(fold)
        model$spec %>%
          fit(x~., training) %>%
          predict(new_data=validation) %>%
          mutate(truth=validation$x) %>%
          rmse(truth, .pred) %>%
          mutate(mse=.estimate^2) %>%
          pull(mse)
      }) %>%
      mean()
  }
}

add_element = function(obj, ...) {
  obj_class = class(obj)
  obj = c(obj, ...)
  class(obj) = obj_class
  return(obj)
}

make_model = function(data, model_specs) {
  model_specs %>% imap(function(model_spec, var_name) {
    model_data = modeling_data(data, var_name)
    if (model_spec$mode == "regression") {
      model_spec %>%
        fit(x~., model_data) %>%
        add_element(sd=honest_sd(., model_data))
    } else {
      model_spec %>%
        fit(x~., mutate(model_data, x=as_factor(x==1)))
    }
  })
}

sample_from = function(model, data) {
  if (model$spec$mode == "regression") {
    rnorm(nrow(data),
          mean = predict(model, as_tibble(data)) %>%
                  pull(.pred),
          sd = model$sd)
  } else {
    rbinom(nrow(data), 1,
          prob = predict(model, as_tibble(data), type="prob") %>%
                  pull(.pred_TRUE))
  }
}

build_model_next = function(model) {
  model = model
  function(data, treat_plan) {
    data_next = data
    for (x in names(model)) {
      if (x == "D") { # to keep people dead
        data_next[,x] = pmax(data[,x], sample_from(model[[x]], data))
      } else if (x == "Y") { # to make costs 0 for the dead
        data_next[,x] = (1-data[,"D"]) * sample_from(model[[x]], data)
      } else {
        data_next[,x] = sample_from(model[[x]], data)
      }
    }
    data_next[,"A"] = treat_plan(data)
    return(data_next)
  }
}

cost = function(data) {
  data %>%
    group_by(t) %>%
    summarize(Y_monthly_mean = mean(Y)) %>%
    pull(Y_monthly_mean) %>%
    sum()
}

causal_contrast = function(sim_next, treat, control, n=1000, tau=100) {
  cost(sim_data(sim_next, treat, n, tau, p=1)) - cost(sim_data(sim_next, control, n, tau, p=1))
}

model = data %>%
  make_model(model_spec)
model_next = model %>%
  build_model_next()
true_tau = causal_contrast(true_next, treat_all, treat_none)
est_tau = causal_contrast(model_next, treat_all, treat_none)

true_tau
est_tau

# boot_causal_contrast = function(data, treat, control, ...) {
#   boot_data = data %>%
#     count(i) %>%
#     sample_n(nrow(.), replace = T) %>%
#     left_join(data, by="i")
#
#   data %>%
#     make_model(model_spec) %>%
#     build_model_next(model_spec) %>%
#     causal_contrast(treat, control, n, ...)
# }
#
# plan(multiprocess)
# est_tau_boot = 1:100 %>%
#   future_map_dbl(~boot_causal_contrast(data, treat_all, treat_none))
# var(est_tau_boot)

gcomp_cost(data,
  var_names = list(
    t="t",
    i="i",
    Y="Y",
    D = "D",
    A = "A"),
  model_specs = list(
    L1 = linear_reg() %>% set_engine("lm"),
    D = logistic_reg() %>% set_engine("glm"),
    Y = rand_forest() %>% set_engine("ranger"),
    default_numeric = linear_reg() %>% set_engine("lm"),
    default_factor = logistic_reg() %>% set_engine("glm")
  )
)

prep_data = function(data, var_names) {
  data %>%
    rename
}

prep_model_specs = function(data, model_specs) {

}

gcomp_cost(data, var_names) {
  model = data %>%
    make_model(model_spec)
  model_next = model %>%
    build_model_next(model_spec)
  true_tau = causal_contrast(true_next, treat_all, treat_none)
  est_tau = causal_contrast(model_next, treat_all, treat_none)

  boot_causal_contrast = function(data, treat, control, ...) {
    boot_data = data %>%
      count(i) %>%
      sample_n(nrow(.), replace = T) %>%
      left_join(data, by="i")

    data %>%
      make_model(model_spec) %>%
      build_model_next(model_spec) %>%
      causal_contrast(treat, control, n, ...)
  }

  plan(multiprocess)
  est_tau_boot = 1:100 %>%
    future_map_dbl(~boot_causal_contrast(data, treat_all, treat_none))
}

gcomp_sim_test = function() {
  data = sim_data(true_next, natural, n, tau, p=1)

  return(true_tau, est_tau, var(est_tau_boot))
}

