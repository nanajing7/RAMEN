
#' SVD initialization for latent factors
#'
#' @param Y Numeric matrix.
#' @param K Latent dimension.
#'
#' @return A list with U and V.
svd_init <- function(Y, K = 2) {
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"

  if (length(K) != 1 || K < 1 || K != as.integer(K)) {
    stop("K must be a positive integer.", call. = FALSE)
  }

  if (K > min(nrow(Y), ncol(Y))) {
    stop(
      sprintf("K = %d exceeds min(nrow(Y), ncol(Y)) = %d.", K, min(nrow(Y), ncol(Y))),
      call. = FALSE
    )
  }

  grand <- mean(Y)
  R <- sweep(sweep(Y, 1, rowMeans(Y), "-"), 2, colMeans(Y), "-") + grand
  s <- svd(R, nu = K, nv = K)

  U <- s$u[, seq_len(K), drop = FALSE] %*% diag(sqrt(s$d[seq_len(K)]), K, K)
  V <- s$v[, seq_len(K), drop = FALSE] %*% diag(sqrt(s$d[seq_len(K)]), K, K)

  rownames(U) <- rownames(Y)
  rownames(V) <- colnames(Y)
  colnames(U) <- paste0("k", seq_len(K))
  colnames(V) <- paste0("k", seq_len(K))

  list(U = U, V = V)
}
