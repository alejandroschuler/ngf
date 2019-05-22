library(tidyverse)
library(magrittr)
library(broom)
library(furrr)
library(zeallot)
library(tidymodels)

as_factor.logical = function(lgl, ...) base::factor(lgl, levels=c("TRUE", "FALSE"), ...)
as_factor.default = function(x, ...) as.factor(x, ...)
fct_to_lgl = function(x) as.logical(2-as.numeric(x)) # amazingly, 3x as fast as as.logical()

rcategorical = function(n, p) { # returns a factor vector
  mc2d::rmultinomial(n=n, size=1, prob=as.matrix(p)) %>%
    equals(1) %>%
    which(arr.ind = T) %>%
    as_tibble() %>%
    arrange(row) %>%
    pull(col) %>%
    as_factor() %>%
    recode(!!!set_names(names(p), 1:length(p)))
}

# To do:
# add multiple confounders to sim
# functionalize naming

##### DGP ####
n = 100 # months
true_init = list(
  disease= rep(0, n),
  grp = rcategorical(n, tibble(a=1, b=1, c=2)) %>% as.character(),
  treatment = rep(0, n), # everyone starts not treated
  death = rep(F, n),
  cost = rep(0,n)
)

natural_treatment = function(data) if_else(data[["death"]]==1,
                                           data[["treatment"]],
                                           pmax(data[["treatment"]],
                                                purrr::rbernoulli(length(data[["treatment"]]), p=0.05))) # 5% chance of being treated each interval, treatment carries over

true_next = function(data, treat_plan) {
  n = length(data[[1]])
  data_next = data
  data_next[["disease"]] = data[["disease"]] + rnorm(n, 0,1) + rep(1,n) # disease progression
  data_next[["grp"]] = data[["grp"]] # category of some kind
  data_next[["treatment"]] = treat_plan(data)
  data_next[["death"]] = pmax(data[["death"]], purrr::rbernoulli(n, p=gtools::inv.logit(data[["disease"]]-10))) # p(death) increases with disease progression
  data_next[["cost"]] = (1-data[["death"]])*(1.1^(data[["disease"]] - 2*(data[["treatment"]]) + abs(rnorm(n, 0,1)))) # disease increases cost
  data_next
}
##### DGP ####

matrix_to_df = function(data_matrix, name) {
  as_tibble(data_matrix) %>%
    set_names(str_c("t",1:ncol(.),sep="_")) %>%
    mutate(i=1:nrow(.)) %>%
    gather(t, !!name, -i) %>%
    separate(t, into=c("trash","t"), sep="_") %>%
    mutate(t=as.integer(t)) %>%
    select(-trash)
}

at_time = function(data, t) {
  data %>% map(~.[,t])
}

prep_sim_data = function(data_init, tau) {
  data = list()
  n = length(data_init[[1]])
  for (x in names(data_init)) {
    if (is.numeric(data_init[[x]])) {
      data[[x]] = matrix(double(), n, tau)
    } else if (is.logical(data_init[[x]])) {
      data[[x]] = matrix(logical(), n, tau)
    } else {
      data[[x]] = matrix(character(), n, tau)
    }
    data[[x]][,1] = data_init[[x]]
  }
  # data_array = model_specs %>%
  #   map_lgl(~.$mode == "regression") %>%
  #   if_else(list(matrix(double(), n, tau)), list(matrix(factor(), n, tau)))
  return(data)
}

sim_data = function(sim_next, treat_plan, data) {
  tau = ncol(data[[1]])
  for (t in 1:(tau-1)){
    data_next = sim_next(at_time(data, t), treat_plan)
    for (x in names(data)) {
      data[[x]][,t+1] = data_next[[x]]
    }
  }
  data %>%
    imap(matrix_to_df)%>%
    reduce(inner_join, by=c("t","i"))
  # dimnames(data_array)[[3]] %>%
  #   map(~matrix_to_df(data_array[,,.], .)) %>%
  #   reduce(inner_join, by=c("t","i"))
}

data = sim_data(true_next, natural_treatment, prep_sim_data(true_init, 100))
cols=list(
  time="t",
  obs="i",
  treatment="treatment",
  death="death",
  cost="cost")

modeling_data = function(data, x, cols=list(time="t", obs="i", death="D")) {
  data %>%
    select(!!sym(cols$obs),!!sym(cols$time),x=!!sym(x)) %>%
    inner_join(data %>%
                 filter(!(!!sym(cols$death))) %>%
                 mutate(t=t+1), # put x(t+1) on the same row as x(t), z1(t), z2(t) etc.
               by=c(cols$obs, cols$time)) %>%
    select(-!!sym(cols$obs), -!!sym(cols$time), -!!sym(cols$death))
}

# data %>% filter(i==1) %>%
#   modeling_data("cost", cols)

model_specs = list(
  disease = rand_forest() %>%
    set_engine("ranger") %>%
    set_mode("regression"),
  grp = rand_forest() %>%
    set_engine("ranger") %>%
    set_mode("classification"),
  death = logistic_reg() %>%
    set_engine("glm"),
  cost = rand_forest() %>%
    set_engine("ranger") %>%
    set_mode("regression")
)

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

# add something to make this efficient for constant variables
make_model = function(data, model_specs, cols) {
  model_specs %>% imap(function(model_spec, var_name) {
    model_data = modeling_data(data, var_name, cols)
    if (model_spec$mode == "regression") {
      model_spec %>%
        fit(x~., model_data) %>%
        add_element(sd=honest_sd(., model_data))
    } else {
      model_spec %>%
        fit(x~., mutate(model_data, x=as_factor(x)))
    }
  })
}

model = data %>%
  make_model(model_specs, cols)

sample_from = function(model, data) {
  if (model$spec$mode == "regression") {
    rnorm(nrow(data),
          mean = predict(model, as_tibble(data)) %>%
                  pull(.pred),
          sd = model$sd)
  } else {
    rcategorical(nrow(data),
                 p = predict(model, as_tibble(data), type="prob") %>%
                    rename_all(str_extract, "(?<=_).*$"))
  }
}

p = predict(model$grp, data, type="prob")

function(fct, )

build_model_next = function(model, cols, death_level) {
  model = model
  function(data, treat_plan) {
    data_next = data
    for (x in names(model)) {
      if (x == cols$death) { # to keep people dead
        data_next[[x]] = pmax(data[[x]], sample_from(model[[x]], data))
      } else if (x == cols$cost) { # to make costs 0 for the dead
        data_next[[x]] = !(data[[cols$death]]==death_level) * sample_from(model[[x]], data)
      } else {
        data_next[[x]] = sample_from(model[[x]], data)
      }
    }
    data_next[[cols$treatment]] = treat_plan(data)
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

treat_all = function(data) rep(1,nrow(data))
treat_none = function(data) rep(0,nrow(data))

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
#
# gcomp_cost(data,
#   var_names = list(
#     t="t",
#     i="i",
#     Y="Y",
#     D = "D",
#     A = "A"),
#   model_specs = list(
#     L1 = linear_reg() %>% set_engine("lm"),
#     D = logistic_reg() %>% set_engine("glm"),
#     Y = rand_forest() %>% set_engine("ranger"),
#     default_numeric = linear_reg() %>% set_engine("lm"),
#     default_factor = logistic_reg() %>% set_engine("glm")
#   )
# )
#
# prep_data = function(data, var_names) {
#   data %>%
#     rename
# }
#
# prep_model_specs = function(data, model_specs) {
#
# }
#
# gcomp_cost(data, var_names) {
#   model = data %>%
#     make_model(model_spec)
#   model_next = model %>%
#     build_model_next(model_spec)
#   true_tau = causal_contrast(true_next, treat_all, treat_none)
#   est_tau = causal_contrast(model_next, treat_all, treat_none)
#
#   boot_causal_contrast = function(data, treat, control, ...) {
#     boot_data = data %>%
#       count(i) %>%
#       sample_n(nrow(.), replace = T) %>%
#       left_join(data, by="i")
#
#     data %>%
#       make_model(model_spec) %>%
#       build_model_next(model_spec) %>%
#       causal_contrast(treat, control, n, ...)
#   }
#
#   plan(multiprocess)
#   est_tau_boot = 1:100 %>%
#     future_map_dbl(~boot_causal_contrast(data, treat_all, treat_none))
# }
#
# gcomp_sim_test = function() {
#   data = sim_data(true_next, natural, n, tau, p=1)
#
#   return(true_tau, est_tau, var(est_tau_boot))
# }
#
