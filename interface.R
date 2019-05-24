# user-facing functions

var_spec = function(time=NULL, observation=NULL, covariates=NULL, treatment=NULL, death=NULL, cost=NULL) {
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

# add an input checker/clearner
# add default model specs (glms) there
# check all entries in var_spec have an entry in model_spec, except trt, t, and obs
validate_ngf_cost_spec = function(x=list()) {
  x
}

new_ngf_cost_spec = function(x = list()) {
  stopifnot(is.list(x))
  structure(x, class="ngf_cost_spec")
}

ngf_cost_spec = function(var_spec, model_spec, death_levels) {
  obj = new_ngf_cost_spec(list(
    var_spec=var_spec,
    model_spec=model_spec,
    death_levels=death_levels))
  validate_ngf_cost_spec(obj)
}

# build default death_levels from unique vals of death col if not already defined
fit.ngf_cost_spec = function(obj, data) {
  transition_model = data %>%
    make_model(obj$model_spec, obj$var_spec)
  transition_sampler = transition_model %>%
    build_sampler(obj$var_spec, obj$death_levels)
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
    conf_int = 2*delta - quantile(boot_est, (1+c(alpha, -alpha))/2)
  )
}

