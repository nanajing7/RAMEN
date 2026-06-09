#' Print method for bootstrap_als_factorize_joint_cov objects
#'
#' @param x An object of class "bootstrap_als_factorize_joint_cov".
#' @param ... Unused.
#'
#' @return Invisibly returns x.
#' @export
print.bootstrap_als_factorize_joint_cov <- function(x, ...) {
  if (!inherits(x, "bootstrap_als_factorize_joint_cov")) {
    stop("x must be a bootstrap_als_factorize_joint_cov object.", call. = FALSE)
  }

  coef_n <- if (!is.null(x$coef_hat)) length(x$coef_hat) else NA_integer_

  latent_dim <- c(NA_integer_, NA_integer_)
  if (!is.null(x$latent_summary) && !is.null(x$latent_summary$mean)) {
    latent_dim <- dim(x$latent_summary$mean)
  }

  cat("Bootstrap inference object for als_factorize_joint_cov\n")
  cat("----------------------------------------------------\n")
  cat("Bootstrap type: ", x$bootstrap_type, "\n", sep = "")
  cat("Requested replications: ", x$B_requested, "\n", sep = "")
  cat("Successful replications: ", x$B_successful, "\n", sep = "")
  cat("Success rate: ",
      sprintf("%.1f%%", 100 * x$B_successful / x$B_requested),
      "\n", sep = "")
  cat("Number of coefficient targets: ", coef_n, "\n", sep = "")
  cat("Latent interaction matrix dimension: ",
      latent_dim[1], " x ", latent_dim[2], "\n", sep = "")

  invisible(x)
}






#' Summary method for bootstrap_als_factorize_joint_cov objects
#'
#' @param object An object of class "bootstrap_als_factorize_joint_cov".
#' @param ... Unused.
#'
#' @return An object of class "summary.bootstrap_als_factorize_joint_cov".
#' @export
summary.bootstrap_als_factorize_joint_cov <- function(object, ...) {
  if (!inherits(object, "bootstrap_als_factorize_joint_cov")) {
    stop("object must be a bootstrap_als_factorize_joint_cov object.", call. = FALSE)
  }

  latent_mean <- object$latent_summary$mean
  latent_sd <- object$latent_summary$sd
  latent_sign <- object$latent_summary$sign.prob

  latent_overview <- list(
    dim = dim(latent_mean),
    mean_of_mean = mean(latent_mean, na.rm = TRUE),
    sd_of_mean = stats::sd(as.vector(latent_mean), na.rm = TRUE),
    mean_of_sd = mean(latent_sd, na.rm = TRUE),
    max_sd = max(latent_sd, na.rm = TRUE),
    mean_sign_prob = mean(latent_sign, na.rm = TRUE)
  )

  out <- list(
    bootstrap_type = object$bootstrap_type,
    B_requested = object$B_requested,
    B_successful = object$B_successful,
    success_rate = object$B_successful / object$B_requested,
    coef_summary = object$coef_summary,
    latent_overview = latent_overview,
    latent_summary = object$latent_summary
  )

  class(out) <- "summary.bootstrap_als_factorize_joint_cov"
  out
}






#' Print method for summary.bootstrap_als_factorize_joint_cov objects
#'
#' @param x An object of class "summary.bootstrap_als_factorize_joint_cov".
#' @param digits Number of digits to print for numeric summaries.
#' @param ... Unused.
#'
#' @return Invisibly returns x.
#' @export
print.summary.bootstrap_als_factorize_joint_cov <- function(
    x,
    digits = 3,
    ...
) {
  if (!inherits(x, "summary.bootstrap_als_factorize_joint_cov")) {
    stop("x must be a summary.bootstrap_als_factorize_joint_cov object.", call. = FALSE)
  }

  cat("Summary of bootstrap inference for als_factorize_joint_cov\n")
  cat("---------------------------------------------------------\n")
  cat("Bootstrap type: ", x$bootstrap_type, "\n", sep = "")
  cat("Requested replications: ", x$B_requested, "\n", sep = "")
  cat("Successful replications: ", x$B_successful, "\n", sep = "")
  cat("Success rate: ",
      sprintf("%.1f%%", 100 * x$success_rate),
      "\n\n", sep = "")

  cat("Latent interaction overview\n")
  cat("--------------------------\n")
  cat("Dimension: ",
      x$latent_overview$dim[1], " x ", x$latent_overview$dim[2], "\n", sep = "")
  cat("Mean of latent means: ",
      formatC(x$latent_overview$mean_of_mean, digits = digits, format = "f"), "\n", sep = "")
  cat("SD of latent means: ",
      formatC(x$latent_overview$sd_of_mean, digits = digits, format = "f"), "\n", sep = "")
  cat("Mean bootstrap SD: ",
      formatC(x$latent_overview$mean_of_sd, digits = digits, format = "f"), "\n", sep = "")
  cat("Max bootstrap SD: ",
      formatC(x$latent_overview$max_sd, digits = digits, format = "f"), "\n", sep = "")
  cat("Mean sign probability: ",
      formatC(x$latent_overview$mean_sign_prob, digits = digits, format = "f"), "\n\n", sep = "")

  cat("Coefficient summary\n")
  cat("-------------------\n")
  print(x$coef_summary, row.names = FALSE)

  invisible(x)
}




#' Print method for bootstrap_temporal_bipartite_als objects
#'
#' @param x An object of class "bootstrap_temporal_bipartite_als".
#' @param ... Unused.
#'
#' @return Invisibly returns x.
#' @export
print.bootstrap_temporal_bipartite_als <- function(x, ...) {
  if (!inherits(x, "bootstrap_temporal_bipartite_als")) {
    stop("x must be a bootstrap_temporal_bipartite_als object.", call. = FALSE)
  }

  periods <- names(x$period_results)
  TT <- length(periods)
  B_req <- if (!is.null(x$params$B)) x$params$B else NA_integer_
  btype <- if (!is.null(x$params$bootstrap_type)) x$params$bootstrap_type else NA_character_

  total_req <- sum(x$success_df$B_requested, na.rm = TRUE)
  total_suc <- sum(x$success_df$B_successful, na.rm = TRUE)

  n_terms <- if (!is.null(x$coef_df) && nrow(x$coef_df) > 0) {
    length(unique(x$coef_df$term))
  } else NA_integer_

  cat("Temporal bootstrap inference object\n")
  cat("-----------------------------------\n")
  cat("Number of periods: ", TT, "\n", sep = "")
  cat("Periods: ", paste(periods, collapse = ", "), "\n", sep = "")
  cat("Bootstrap type: ", btype, "\n", sep = "")
  cat("Replications per period (requested): ", B_req, "\n", sep = "")
  cat("Total successful replications: ", total_suc, " / ", total_req, "\n", sep = "")
  if (total_req > 0) {
    cat("Overall success rate: ",
        sprintf("%.1f%%", 100 * total_suc / total_req), "\n", sep = "")
  }
  cat("Number of coefficient terms: ", n_terms, "\n", sep = "")
  cat("\nUse summary() for stacked mean + CI tables.\n")

  invisible(x)
}




#' Summary method for bootstrap_temporal_bipartite_als objects
#'
#' @param object An object of class "bootstrap_temporal_bipartite_als".
#' @param ... Unused.
#'
#' @return An object of class "summary.bootstrap_temporal_bipartite_als".
#' @export
summary.bootstrap_temporal_bipartite_als <- function(object, ...) {
  if (!inherits(object, "bootstrap_temporal_bipartite_als")) {
    stop("object must be a bootstrap_temporal_bipartite_als object.", call. = FALSE)
  }

  latent_overview <- NULL
  if (!is.null(object$latent_df) && nrow(object$latent_df) > 0) {
    latent_overview <- data.frame(
      period           = unique(object$latent_df$period),
      n_cells          = as.integer(tapply(object$latent_df$mean,
                                           object$latent_df$period, length)),
      mean_of_mean     = as.numeric(tapply(object$latent_df$mean,
                                           object$latent_df$period,
                                           mean, na.rm = TRUE)),
      mean_of_sd       = as.numeric(tapply(object$latent_df$sd,
                                           object$latent_df$period,
                                           mean, na.rm = TRUE)),
      mean_sign_prob   = as.numeric(tapply(object$latent_df$sign.prob,
                                           object$latent_df$period,
                                           mean, na.rm = TRUE)),
      row.names        = NULL,
      stringsAsFactors = FALSE
    )
  }

  out <- list(
    params          = object$params,
    success_df      = object$success_df,
    coef_df         = object$coef_df,
    latent_df       = object$latent_df,
    latent_overview = latent_overview
  )

  class(out) <- "summary.bootstrap_temporal_bipartite_als"
  out
}




#' Print method for summary.bootstrap_temporal_bipartite_als objects
#'
#' @param x An object of class "summary.bootstrap_temporal_bipartite_als".
#' @param digits Number of digits to print for numeric summaries.
#' @param ... Unused.
#'
#' @return Invisibly returns x.
#' @export
print.summary.bootstrap_temporal_bipartite_als <- function(
    x,
    digits = 3,
    ...
) {
  if (!inherits(x, "summary.bootstrap_temporal_bipartite_als")) {
    stop("x must be a summary.bootstrap_temporal_bipartite_als object.",
         call. = FALSE)
  }

  btype <- if (!is.null(x$params$bootstrap_type)) x$params$bootstrap_type else NA_character_
  B_req <- if (!is.null(x$params$B)) x$params$B else NA_integer_
  cl    <- if (!is.null(x$params$conf_level)) x$params$conf_level else NA_real_

  cat("Summary of temporal bootstrap inference\n")
  cat("---------------------------------------\n")
  cat("Bootstrap type: ", btype, "\n", sep = "")
  cat("Replications per period (requested): ", B_req, "\n", sep = "")
  cat("Confidence level: ",
      if (!is.na(cl)) sprintf("%.1f%%", 100 * cl) else NA, "\n", sep = "")

  cat("\nSuccess rates by period\n")
  cat("-----------------------\n")
  print(x$success_df, row.names = FALSE)

  cat("\nCoefficient mean + CI (stacked across periods)\n")
  cat("----------------------------------------------\n")
  print(x$coef_df, row.names = FALSE, digits = digits)

  if (!is.null(x$latent_overview)) {
    cat("\nLatent interaction overview by period\n")
    cat("-------------------------------------\n")
    print(x$latent_overview, row.names = FALSE, digits = digits)
    cat("\n(Full cell-level mean + CI available in $latent_df.)\n")
  }

  invisible(x)
}

