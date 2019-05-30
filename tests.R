source("g-comp.R")
source("interface.R")

##### DGP #####
make_true_init = function(n=100) {
  tibble(
    disease= rep(0, n),
    grp = rcategorical(n, tibble(a=1, b=1, c=2)),
    treatment = rep(0, n), # everyone starts not treated
    death = rep(F, n),
    cost = rep(0,n)
  )
}

natural_treatment = function(data) if_else(data[["death"]]==1,
                                           data[["treatment"]],
                                           pmax(data[["treatment"]],
                                                purrr::rbernoulli(length(data[["treatment"]]), p=0.05))) # 5% chance of being treated each interval, treatment carries over

true_transition_sampler = function(data, treat_plan, var_specs) {
  n = length(data[[1]])
  data_next = data
  data_next[["disease"]] = data[["disease"]] + rnorm(n, 0,1) + rep(1,n) # disease progression
  data_next[["grp"]] = data[["grp"]] # category of some kind
  data_next[["treatment"]] = treat_plan(data)
  data_next[["death"]] = pmax(data[["death"]], purrr::rbernoulli(n, p=gtools::inv.logit(data[["disease"]]-10))) %>% as.logical() # p(death) increases with disease progression
  data_next[["cost"]] = (1-data[["death"]])*(1.3^(data[["disease"]] + 1.5*(data[["treatment"]]) + abs(rnorm(n, 0,1)))) # disease increases cost
  data_next
}

##### Testing data #####

data = make_true_init(1000) %>% # 1000 patients
  prep_sim_data(24) %>% # for 24 months
  sim_data(true_transition_sampler, natural_treatment)

##### Inputs  #####

model_specs = list(
  disease = rand_forest() %>%
    set_engine("ranger") %>%
    set_mode("regression"),
  grp = static(),
  death = rand_forest() %>%
    set_engine("ranger") %>%
    set_mode("classification"),
  cost = rand_forest() %>%
    set_engine("ranger") %>%
    set_mode("regression")
)

var_specs=make_var_spec(
  time="t", # time and obs are only used by data making for model fitting
  observation = "obs",
  covariates = c("grp", "disease"),
  treatment="treatment",
  death="death", # this is needed for data making too though
  cost="cost")

# make the definition of these functions clearer
treat_all = function(data) rep(1,length(data[[1]]))
treat_none = function(data) rep(0,length(data[[1]]))

##### Tests #####

true_delta = list(
    model=NULL,
    sampler=true_transition_sampler,
    t0_data=make_true_init(100000),
    vars=var_specs,
    t_max=24) %>%
  structure(class = "ngf_cost_fit") %>%
  causal_contrast(treat_all, treat_none, n=10000)

plan(multiprocess)
tic()
delta = ngf_cost_spec(var_specs, model_specs) %>%
  estimate(data, treat_all, treat_none, B=8, n=100, t_max=12)
toc()

# spec = ngf_cost_spec(var_specs, model_specs)
# fit_model = spec %>% fit(data)
# pred = fit_model %>% predict(treat_none, n=100)
