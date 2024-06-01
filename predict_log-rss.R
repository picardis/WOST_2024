pred_log_rss <- function(mod, x1, x2) {

  pred_X <- model.matrix(mod, data = x1)
  ref_X <- model.matrix(mod, data = x2)
  B <- coef(mod)

  # Difference
  diff_X <- pred_X - ref_X

  # Linear predictor for difference
  diff_dist <- unname((diff_X %*% B)[, 1])

  # Get variance of difference
  var_pred <- diag(diff_X %*% vcov(mod) %*% t(diff_X))

  # Standard error for difference
  SE <- unname(sqrt(var_pred))

  out <- x1 %>%
    mutate(
      log_rss = diff_dist,
      se = SE,
      lwr = log_rss - 1.96 * se,
      upr = log_rss + 1.96 * se
    )

  return(out)

  }
