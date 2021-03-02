library(tidymodels)
source("functions.R")

IRLS <- function(mode = "classification", lambda=NULL, beta0 = NULL, 
                 eps = 0.0001, iter_max = 100) {
  
  args <- list(
      lambda = rlang::enquo(lambda),
      beta0 = rlang::enquo(beta0),
      eps = rlang::enquo(eps),
      iter_max = rlang::enquo(iter_max)
    )
  
  new_model_spec("IRLS",
                 args = args,
                 mode = mode,
                 eng_args = NULL,
                 method = NULL,
                 engine = NULL)
}
set_new_model("IRLS")
set_model_mode(model = "IRLS", mode = "classification")
set_model_engine("IRLS",
                 mode = "classification",
                 eng = "fit_logistic_lasso"
)
set_dependency("IRLS", eng = "fit_logistic_lasso", pkg = "base")
set_encoding(
  model = "IRLS",
  eng = "fit_logistic_lasso",
  mode = "classification",
  options = list(
    predictor_indicators = "traditional",
    compute_intercept = TRUE,
    remove_intercept = TRUE,
    allow_sparse_x = FALSE
  )
)
show_model_info("IRLS")

set_model_arg(
  model        = "IRLS", 
  eng          = "fit_logistic_lasso", 
  parsnip      = "lambda", 
  original     = "lambda", 
  func         = list(pkg = "foo", fun = "bar"),
  has_submodel = FALSE
)
set_model_arg(
  model        = "IRLS", 
  eng          = "fit_logistic_lasso", 
  parsnip      = "beta0", 
  original     = "beta0", 
  func         = list(pkg = "foo", fun = "bar"),
  has_submodel = FALSE
)
set_model_arg(
  model        = "IRLS", 
  eng          = "fit_logistic_lasso", 
  parsnip      = "eps", 
  original     = "eps", 
  func         = list(pkg = "foo", fun = "bar"),
  has_submodel = FALSE
)
set_model_arg(
  model        = "IRLS", 
  eng          = "fit_logistic_lasso", 
  parsnip      = "iter_max", 
  original     = "iter_max", 
  func         = list(pkg = "foo", fun = "bar"),
  has_submodel = FALSE
)
set_fit(
  model = "IRLS",
  eng = "fit_logistic_lasso",
  mode = "classification",
  value = list(
    interface = "matrix",
    protect = c("x", "y"),
    func = c(fun = "fit_logistic_lasso"),
    defaults = list()
  )
)
set_pred(
  model = "IRLS",
  eng = "fit_logistic_lasso",
  mode = "classification",
  type = "class",
  value = list(
    pre = NULL,
    post = NULL,
    func = c(fun = "predict_logistic_lasso"),
    args = list(
      object = expr(object$fit),
      new_x = expr(as.matrix(new_data[, names(object$fit$beta)]))
    )
  )
)


