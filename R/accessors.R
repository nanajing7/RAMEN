#' Validate bootstrap ALS object
#'
#' @param object An object expected to inherit from
#'   "bootstrap_als_factorize_joint_cov".
#'
#' @return Invisibly returns TRUE if validation passes.
.validate_bootstrap_als_object <- function(object) {
  if (!inherits(object, "bootstrap_als_factorize_joint_cov")) {
    stop(
      "object must be of class 'bootstrap_als_factorize_joint_cov'.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}


#' Extract coefficient summary table from bootstrap ALS object
#'
#' @param object A bootstrap object returned by bootstrap_als_factorize_joint_cov().
#'
#' @return A data.frame of coefficient summaries.
#' @export
coef_summary_bootstrap_als <- function(object) {
  .validate_bootstrap_als_object(object)
  object$coef_summary
}


#' Extract coefficient draws from bootstrap ALS object
#'
#' @param object A bootstrap object returned by bootstrap_als_factorize_joint_cov().
#'
#' @return A matrix of bootstrap coefficient draws.
#' @export
coef_draws_bootstrap_als <- function(object) {
  .validate_bootstrap_als_object(object)
  object$coef_draws
}


#' Extract original coefficient estimates from bootstrap ALS object
#'
#' @param object A bootstrap object returned by bootstrap_als_factorize_joint_cov().
#'
#' @return A named numeric vector of original coefficient estimates.
#' @export
coef_hat_bootstrap_als <- function(object) {
  .validate_bootstrap_als_object(object)
  object$coef_hat
}


#' Extract bootstrap mean latent interaction matrix
#'
#' @param object A bootstrap object returned by bootstrap_als_factorize_joint_cov().
#'
#' @return A matrix of bootstrap means for U V^T.
#' @export
latent_mean_bootstrap_als <- function(object) {
  .validate_bootstrap_als_object(object)
  object$latent_summary$mean
}


#' Extract bootstrap standard deviation matrix for latent interactions
#'
#' @param object A bootstrap object returned by bootstrap_als_factorize_joint_cov().
#'
#' @return A matrix of bootstrap standard deviations for U V^T.
#' @export
latent_sd_bootstrap_als <- function(object) {
  .validate_bootstrap_als_object(object)
  object$latent_summary$sd
}


#' Extract bootstrap sign-probability matrix for latent interactions
#'
#' @param object A bootstrap object returned by bootstrap_als_factorize_joint_cov().
#'
#' @return A matrix giving Pr(UV^T > 0) across bootstrap draws.
#' @export
latent_sign_prob_bootstrap_als <- function(object) {
  .validate_bootstrap_als_object(object)
  object$latent_summary$sign.prob
}


#' Extract bootstrap confidence interval matrices for latent interactions
#'
#' @param object A bootstrap object returned by bootstrap_als_factorize_joint_cov().
#'
#' @return A list with elements conf.low and conf.high.
#' @export
latent_ci_bootstrap_als <- function(object) {
  .validate_bootstrap_als_object(object)
  list(
    conf.low = object$latent_summary$conf.low,
    conf.high = object$latent_summary$conf.high
  )
}


#' Extract bootstrap latent interaction draws
#'
#' @param object A bootstrap object returned by bootstrap_als_factorize_joint_cov().
#'
#' @return A list of latent interaction matrices, one per successful bootstrap replication.
#' @export
latent_draws_bootstrap_als <- function(object) {
  .validate_bootstrap_als_object(object)
  object$latent_draws
}


#' Extract bootstrap success information
#'
#' @param object A bootstrap object returned by bootstrap_als_factorize_joint_cov().
#'
#' @return A list containing requested replications, successful replications,
#'   success rate, and the logical success vector.
#' @export
bootstrap_diagnostics_als <- function(object) {
  .validate_bootstrap_als_object(object)

  list(
    B_requested = object$B_requested,
    B_successful = object$B_successful,
    success_rate = object$B_successful / object$B_requested,
    success = object$success
  )
}
