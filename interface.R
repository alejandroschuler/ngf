# user-facing functions

# To do:
# - Docs

make_var_spec = function(time=NULL, observation=NULL, covariates=NULL, treatment=NULL, death=NULL, cost=NULL) {
  if (list(time, observation, treatment, death, cost) %>%
      map_lgl(is.null) %>%
      any()) {
    rlang::abort("the names of the time, observation, treatment, death, and cost variables must be sepcified")
  }
  list(
    time=time,
    observation=observation,
    treatment=treatment,
    death=death,
    cost=cost,
    covariates=covariates) %>%
    purrr::discard(is.null)
}

validate_ngf_cost_spec = function(spec = list()) {
  spec$vars = exec(make_var_spec, !!!(spec$vars)) # checks that the right variables are here and named correctly
  needed_model_names = spec$vars %>%
    extract(!(names(.) %in% c("time","observation","treatment"))) %>%
    unlist()
  if (is.null(spec$models)) {
    warning("No model specification provided. Defaulting to GLMs")
    spec$models = list(
      logistic_reg() %>% set_engine("glm"),
      linear_reg() %>% set_engine("lm")) %>%
    set_names(spec$vars[c("death","cost")])
  }
  # check mode of death is classification and mode of cost is regression
  if (any(!(names(spec$models) %in% unlist(spec$vars)))) {
    rlang::abort("Variables in the model spec are not identified in the variable spec")
  }
  return(spec)
}

new_ngf_cost_spec = function(x = list()) {
  stopifnot(is.list(x))
  structure(x, class="ngf_cost_spec")
}

ngf_cost_spec = function(var_spec, model_spec=NULL, dead_level=NULL) {
  list(
    vars=var_spec,
    models=model_spec,
    dead=dead_level) %>%
  validate_ngf_cost_spec() %>% # checks
  new_ngf_cost_spec() # classes the list
}

check_data = function(spec, data) {
  data = data %>%
    mutate_if(negate(is.numeric), as.character)

  death_col = data %>%
    arrange(!!sym(spec$vars$time)) %>%
    pull(spec$vars$death)
  if (length(unique(death_col)) > 2) {
    rlang::abort("Death column has more than two unique values")
  }
  cost_col = data %>%
    pull(spec$vars$cost)
  if (!is.numeric(cost_col)) {
    rlang::abort("Cost column is not numeric")
  }

  if (any(!(unlist(spec$vars) %in% names(data)))) {
    rlang::abort("Some variables in the variable spec are not found as column names in the data")
  }

  if (is.null(spec$dead)) {
    rlang::warn("No death level provided, guessing based on data")
    life_level = first(death_col)
    spec$dead = unique(death_col) %>%
      purrr::discard(~equals(.,life_level))  %>%
      as.character()
  }

  # guess model specs for covariates  based on data type
  for (covariate in spec$vars$covariates) {
    if (is.null(spec$models[[covariate]])) {
      rlang::warn(str_c("No model spec provided for ", covariate))
      if (is.numeric(data[[covariate]]))  {
        rlang::warn("Defaulting to linear regression based on datatype")
        spec$models[[covariate]] = linear_reg() %>% set_engine("lm")
      } else {
        rlang::warn("Defaulting to logistic regression based on datatype")
        # Make sure this works for multiclass regression?
        spec$models[[covariate]] = logistic_reg() %>% set_engine("glm")
      }
    }
  }

  list(spec, data)
}

fit.ngf_cost_spec = function(spec, data) {
  c(spec, data) %<-% check_data(spec, data)

  transition_model = spec %>%
    make_model(data)
  transition_sampler = transition_model %>%
    build_sampler(spec$vars, spec$dead)
  t0_data = data %>%
    filter(!!sym(spec$vars$time)==min(!!sym(spec$vars$time))) %>%
    select(spec$vars %$% c(covariates, treatment, death, cost))
  structure(
    list(model=transition_model, sampler=transition_sampler,
         t0_data=t0_data, vars=spec$vars,
         t_max=length(unique(data[[spec$vars$time]]))),
    class = "ngf_cost_fit")
}

predict.ngf_cost_fit = function(fit_model, policy, t_max=fit_model$t_max, ...) {
  fit_model$t0_data %>%
    prep_sim_data(t_max, ...) %>%
    sim_data(fit_model$sampler, policy) %>%
    cost(fit_model$vars)
}

causal_contrast = function(model, treat, control, ...) {
  predict(model, treat, ...) - predict(model, control, ...)
}

boot_sample = function(data, obs_var) {
  data %>%
    count(!!sym(obs_var)) %>%
    select(-n) %>%
    sample_n(nrow(.), replace = T) %>%
    left_join(data, by=obs_var)
}

estimate = function(spec, data, treat, control, B=100, alpha=0.95, ...) {
  delta = spec %>%
    fit(data) %>%
    causal_contrast(treat, control, ...)

  boot_est = 1:B %>%
    future_map_dbl(function(b) {
        fit(spec, boot_sample(data, spec$vars$observation)) %>%
        causal_contrast(treat, control, ...)
    })
  list(
    estimate = delta,
    conf_int = 2*delta - quantile(boot_est, (1+c(-alpha, alpha))/2) %>%
      set_names(rev(names(.)))
  )
}

