# Steps 7 and 8: pin down the two gauge freedoms the observation model leaves in
# the latent factors, so that the plug-in variance components are well defined.
#
# Step 7 removes the common orthogonal rotation: U_t V_t' is unchanged by
# (U_t R, V_t R) for any orthogonal R, and so are all the initial-state and
# temporal penalties, so the MAP objective cannot select an orientation. A
# single rotation aligned to the Step-1 reference picks a reproducible one.
#
# Step 8 removes the reciprocal scale: U_t V_t' is unchanged by (c U_t, V_t / c),
# but sigma_U^2 and tau_U^2 scale by c^2 while sigma_V^2 and tau_V^2 scale by
# c^-2. Left alone, one factor drifts smaller and the other larger across outer
# iterations while the fit stays identical, corrupting the variance updates.
#
# ── Ordering: scale first, then rotate ──────────────────────────────────────
# The paper presents these as Step 7 (rotate) then Step 8 (rescale). This
# implementation applies them the other way round. The two orders produce
# identical fitted values and identical variance components — the choice is
# purely about which representative of the equivalence class gets reported —
# but only scale-then-rotate is idempotent:
#
#   * Rotation preserves Frobenius norms, so applying it after the rescaling
#     leaves the Step-8 scale convention intact. The result satisfies both
#     gauges simultaneously.
#   * Rescaling, however, changes the relative weight of U and V inside the
#     Procrustes objective sum_t{||U_t R - Uref||^2 + ||V_t R - Vref||^2},
#     whose solution is driven by sum_t U_t'Uref + sum_t V_t'Vref. Rescaling
#     after rotating turns that into c * sum U_t'Uref + (1/c) * sum V_t'Vref,
#     so a second pass would return a *different* rotation.
#
# With rotate-then-rescale the reported latent coordinates therefore keep
# turning slightly from one outer iteration to the next. Omega and the fitted
# values are both rotation-invariant, so this never affected convergence or any
# estimate — but it made the latent coordinates themselves irreproducible.
# Scale-then-rotate removes that wobble.
#
# NOTE: the paper's Steps 7 and 8 should be swapped to match.
#
# Both transformations leave the fitted values exactly unchanged, so this step
# never undoes the progress made by the inner loop. The rotation additionally
# leaves Q unchanged; the rescaling does not, and the spec does not claim it
# does — it runs after the inner loop has converged and before the penalties are
# refreshed, so no monotonicity is broken.

#' Identify the latent trajectories (Steps 7 and 8)
#'
#' Fixes both gauge freedoms of the latent factors: first the reciprocal scale,
#' so that the average squared coordinate magnitudes of the two node sets agree,
#' then the common orthogonal transform that best aligns the trajectories to the
#' Step-1 reference.
#'
#' Applying the scale before the rotation makes the operation idempotent, so the
#' reported coordinates settle rather than turning slightly on every outer
#' iteration; see the file header for why the reverse order does not. Fitted
#' values and variance components are the same either way.
#'
#' @param params Parameter list from the inner BCD.
#' @param reference List with `U` and `V`, the temporally aligned unperturbed
#'   initialization retained from Step 1.
#' @param anchor Which magnitudes fix the scale factor: the initial states
#'   alone (`"first"`) or the average over all periods (`"pooled"`). Pooling is
#'   not neutral in a dynamic model — if one node set genuinely moves more from
#'   period to period it also spreads further, and the pooled factor absorbs
#'   part of that asymmetry, pulling `tau_U^2 / tau_V^2` toward one. See
#'   `.scale_normalize_UV()`.
#'
#' @return A list with the transformed `params`, the rotation `R`, and the
#'   scale factor `c` that was applied.
#' @keywords internal
#' @noRd
.ame_identify <- function(params, reference, anchor = c("pooled", "first", "innovation")) {
  anchor <- match.arg(anchor)
  Tn <- length(params$U)

  if (length(reference$U) != Tn || length(reference$V) != Tn) {
    stop("reference must cover the same number of periods as params.",
         call. = FALSE)
  }

  # ── Step 8: reciprocal scale normalization ──────────────────────────────
  sn <- .scale_normalize_UV(params$U, params$V, anchor = anchor)

  # ── Step 7: one common rotation, fitted over U and V and all periods ─────
  # Frobenius norms are invariant under R, so the scale convention just set
  # survives this.
  R <- .orthogonal_procrustes(
    Alist = c(sn$U, sn$V),
    Blist = c(reference$U, reference$V)
  )

  U_new <- vector("list", Tn)
  V_new <- vector("list", Tn)
  for (t in seq_len(Tn)) {
    U_new[[t]] <- sn$U[[t]] %*% R
    V_new[[t]] <- sn$V[[t]] %*% R
    dimnames(U_new[[t]]) <- dimnames(params$U[[t]])
    dimnames(V_new[[t]]) <- dimnames(params$V[[t]])
  }

  params$U <- U_new
  params$V <- V_new

  list(params = params, R = R, c = sn$c)
}
