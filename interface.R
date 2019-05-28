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

validate_ngf_cost_spec = function(x = list()) {
  x$var_spec = exec(make_var_spec, !!!(x$var_spec)) # checks that the right variables are here and named correctly
  needed_model_names = x$var_spec %>%
    extract(!(names(.) %in% c("time","observation","treatment"))) %>%
    unlist()
  if (is.null(x$model_spec)) {
    warning("No model specification provided. Defaulting to GLMs")
    x$model_spec = list(
      logistic_reg() %>% set_engine("glm"),
      linear_reg() %>% set_engine("lm")) %>%
    set_names(x$var_spec[c("death","cost")])
  }
  # check mode of death is classification and mode of cost is regression
  if (any(!(names(x$model_spec) %in% unlist(x$var_spec)))) {
    rlang::abort("Variables in the model spec are not identified in the variable spec")
  }
  x
}

new_ngf_cost_spec = function(x = list()) {
  stopifnot(is.list(x))
  structure(x, class="ngf_cost_spec")
}

ngf_cost_spec = function(var_spec, model_spec=NULL, dead_level=NULL) {
  obj = new_ngf_cost_spec(list(
    var_spec=var_spec,
    model_spec=model_spec,
    dead_level=dead_level))
  validate_ngf_cost_spec(obj)
}

check_data = function(spec_obj, data) {
  death_col = data %>%
    pull(spec_obj$var_spec$death)
  if (length(unique(death_col)) > 2) {
    rlang::abort("Death column has more than two unique values")
  }
  cost_col = data %>%
    pull(spec_obj$var_spec$cost)
  if (!is.numeric(cost_col)) {
    rlang::abort("Cost column is not numeric")
  }

  if (any(!(unlist(spec_obj$var_spec) %in% names(data)))) {
    rlang::abort("Some variables in the variable spec are not found as column names in the data")
  }

  if (is.null(spec_obj$dead_level)) {
    rlang::warn("No death level provided, guessing based on data")
    life_level = first(death_col)
    spec_obj$dead_level = unique(death_col) %>%
      purrr::discard(~equals(.,life_level))  %>%
      as.character()
  }

  # guess model specs for covariates  based on data type
  for (covariate in spec_obj$var_spec$covariates) {
    if (is.null(spec_obj$model_spec[[covariate]])) {
      rlang::warn(str_c("No model spec provided for ", covariate))
      if (is.numeric(data[[covariate]]))  {
        rlang::warn("Defaulting to linear regression based on datatype")
        spec_obj$model_spec[[covariate]] = linear_reg() %>% set_engine("lm")
      } else {
        rlang::warn("Defaulting to logistic regression based on datatype")
        # Make sure this works for multiclass regression?
        spec_obj$model_spec[[covariate]] = logistic_reg() %>% set_engine("glm")
      }
    }
  }

  return(spec_obj)
}

fit.ngf_cost_spec = function(obj, data) {
  obj = check_data(obj, data)

  transition_model = data %>%
    make_model(obj$model_spec, obj$var_spec)
  transition_sampler = transition_model %>%
    build_sampler(obj$var_spec, obj$dead_level)
  t0_data = data %>%
    filter(!!sym(obj$var_spec$time)==min(!!sym(obj$var_spec$time))) %>%
    select(obj$var_spec %$% c(covariates, treatment, death, cost))
  structure(
    list(model=transition_model, sampler=transition_sampler,
         t0_data=t0_data, var_spec=obj$var_spec,
         t_max=length(unique(data[[obj$var_spec$time]]))),
    class = "ngf_cost_fit")
}

predict.ngf_cost_fit = function(obj, policy, t_max=obj$t_max, ...) {
  obj$t0_data %>%
    prep_sim_data(t_max, ...) %>%
    sim_data(obj$sampler, policy) %>%
    cost(obj$var_spec)
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
        fit(spec, boot_sample(data, spec$var_spec$obs)) %>%
        causal_contrast(treat, control, ...)
    })
  list(
    estimate = delta,
    conf_int = 2*delta - quantile(boot_est, (1+c(-alpha, alpha))/2) %>%
      set_names(rev(names(.)))
  )
}

