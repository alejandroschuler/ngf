library(tidyverse)
library(magrittr)
library(broom)
library(furrr)
library(zeallot)
library(tidymodels)

as_factor.logical = function(lgl, ...) base::factor(lgl, levels=c("TRUE", "FALSE"), ...)
as_factor.default = function(x, ...) as.factor(x, ...)
fct_to_lgl = function(x) as.logical(2-as.numeric(x)) # amazingly, 3x as fast as as.logical()

static = function() {
  x = 0
  class(x) = c("static")
  return(x)
}

rcategorical = function(n, p) { # returns a character vector (more portable than factor...)
  mc2d::rmultinomial(n=n, size=1, prob=as.matrix(p)) %>%
    equals(1) %>%
    which(arr.ind = T) %>%
    as_tibble() %>%
    arrange(row) %>%
    pull(col) %>%
    as.character() %>%
    recode(!!!set_names(names(p), 1:length(p)))
}

matrix_to_df = function(data_matrix, name) {
  as_tibble(data_matrix) %>%
    set_names(str_c("t",1:ncol(.),sep="_")) %>%
    mutate(obs=1:nrow(.)) %>%
    gather(t, !!name, -obs) %>%
    separate(t, into=c("trash","t"), sep="_") %>%
    mutate(t=as.integer(t)) %>%
    select(-trash)
}

at_time = function(data, t) {
  data %>% map(~.[,t])
}

prep_sim_data = function(data_init, tau, n=nrow(data_init), sample=T) {
  if (sample) {
    data_init = sample_n(data_init, n, replace = T)
  }
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
  return(data)
}

sim_data = function(data, transition_sampler, treat_plan) {
  tau = ncol(data[[1]])
  for (t in 1:(tau-1)){
    data_next = transition_sampler(at_time(data, t), treat_plan, var_specs)
    for (x in names(data)) {
      data[[x]][,t+1] = data_next[[x]]
    }
  }
  data %>%
    imap(matrix_to_df)%>%
    reduce(inner_join, by=c("t","obs"))
}

modeling_data = function(data, x, var_specs) {
  data %>%
    select(!!sym(var_specs$obs),!!sym(var_specs$time),x=!!sym(x)) %>%
    inner_join(data %>%
                 filter(!(!!sym(var_specs$death))) %>%
                 mutate(t=t+1), # put x(t+1) on the same row as x(t), z1(t), z2(t) etc.
               by=c(var_specs$obs, var_specs$time)) %>%
    select(-!!sym(var_specs$obs), -!!sym(var_specs$time), -!!sym(var_specs$death))
}

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

make_model = function(data, model_specs, var_specs) {
  model_specs %>% imap(function(model_spec, var_name) {
    if (class(model_spec)[[1]] == "static") {
      return(model_spec)
    }
    model_data = modeling_data(data, var_name, var_specs)
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

sample_from = function(model, data) {
  n = length(data[[1]])
  if (model$spec$mode == "regression") {
    rnorm(n,
      mean = predict(model, data) %>%
        pull(.pred),
      sd = model$sd)
  } else {
    rcategorical(n,
      p = predict(model, data, type="prob") %>%
        rename_all(str_extract, "(?<=_).*$"))
  }
}

build_sampler = function(model, var_specs, death_levels) {
  model = model # weird but doesn't work without this.
  function(data, treat_plan, var_specs) {
    data_next = data # takes care of carry-over unless otherwise specified
    dead_index = (data[[var_specs$death]]==death_levels$dead) # only bother predicting for those still alive
    data_next[[var_specs$cost]][dead_index] = 0 # dead cost 0
    data_next[[var_specs$treatment]] = treat_plan(data) # assign treatment
    if (all(dead_index)) {
      return(data_next) # exit early if nobody is around to predict on
    }
    for (x in names(model)) {
      if (class(model[[x]])[[1]] != "static") { # don't bother predicting for static covariates
        data_next[[x]][!dead_index] = sample_from(model[[x]], as_tibble(data)[!dead_index,])
      }
    }
    return(data_next)
  }
}

cost = function(data, var_specs) {
  data %>%
    group_by(t) %>%
    summarize(Y_monthly_mean = mean(!!sym(var_specs$cost))) %>%
    pull(Y_monthly_mean) %>%
    sum()
}
