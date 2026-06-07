### Lotka-Volterra model for 2 mutualistic species
### dN1/dt = N1 * (r1 - alpha11*N1 + gamma12*N2)
### dN2/dt = N2 * (r2 + gamma21*N1 - alpha22*N2)
###
### Two qualitative outcomes:
###   (1) Stable positive equilibrium   (a11*a22 > gamma12*gamma21)
###   (2) Unstable positive equilibrium (a11*a22 < gamma12*gamma21)

library(deSolve)

# ─────────────────────────────────────────────
# ODE function
# ─────────────────────────────────────────────
LV_mutualism <- function(time, state, params) {
  N1 <- state["N1"]
  N2 <- state["N2"]
  
  r1     <- params["r1"]
  r2     <- params["r2"]
  a11    <- params["a11"]
  a22    <- params["a22"]
  g12    <- params["g12"]   # gamma12: benefit of N2 on N1
  g21    <- params["g21"]   # gamma21: benefit of N1 on N2
  
  dN1 <- N1 * (r1 - a11 * N1 + g12 * N2)
  dN2 <- N2 * (r2 + g21 * N1 - a22 * N2)
  
  return(list(c(dN1, dN2)))
}

# ─────────────────────────────────────────────
# Helper: run ODE for one set of initial conditions
# ─────────────────────────────────────────────
run_mut <- function(N1_0, N2_0, params,
                    t = seq(0, 150, by = 0.2)) {
  pop <- c(N1 = N1_0, N2 = N2_0)
  out <- tryCatch(
    ode(pop, t, LV_mutualism, parms = params,
        method = "lsoda", atol = 1e-8, rtol = 1e-8),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  df <- as.data.frame(out)
  # Trim trajectory once either population explodes beyond plotting range
  # (relevant for unstable case)
  df
}

# ─────────────────────────────────────────────
# Helper: vector field (normalised arrows)
# ─────────────────────────────────────────────
vector_field <- function(params, N_max, n_grid = 18, scale = 0.06) {
  N_seq <- seq(0.2, N_max, length.out = n_grid)
  grid  <- expand.grid(N1 = N_seq, N2 = N_seq)
  
  dN1_vec <- numeric(nrow(grid))
  dN2_vec <- numeric(nrow(grid))
  
  for (k in seq_len(nrow(grid))) {
    state_k <- c(N1 = grid$N1[k], N2 = grid$N2[k])
    derivs  <- LV_mutualism(0, state_k, params)[[1]]
    dN1_vec[k] <- derivs[1]
    dN2_vec[k] <- derivs[2]
  }
  
  mag          <- sqrt(dN1_vec^2 + dN2_vec^2)
  mag[mag == 0] <- 1
  dN1_norm <- scale * dN1_vec / mag
  dN2_norm <- scale * dN2_vec / mag
  
  list(grid = grid, dN1 = dN1_norm, dN2 = dN2_norm)
}

# ─────────────────────────────────────────────
# Helper: plot one scenario
# ─────────────────────────────────────────────
plot_scenario <- function(params, ic_list, title,
                          N_max    = NULL,
                          n_grid   = 18,
                          vf_scale = 0.06,
                          t_seq    = seq(0, 150, by = 0.2)) {
  
  r1  <- params["r1"];  r2  <- params["r2"]
  a11 <- params["a11"]; a22 <- params["a22"]
  g12 <- params["g12"]; g21 <- params["g21"]
  
  # ── Isocline geometry ─────────────────────────────────────────────────
  # N1-isocline (dN1/dt = 0, N1 > 0):
  #   r1 - a11*N1 + g12*N2 = 0
  #   N2 = (a11*N1 - r1) / g12          [positive slope: a11/g12]
  #   hits N1-axis (N2=0) at: N1 = r1/a11
  #   hits N2-axis (N1=0) at: N2 = -r1/g12  (can be negative if r1 > 0)
  #
  # N2-isocline (dN2/dt = 0, N2 > 0):
  #   r2 + g21*N1 - a22*N2 = 0
  #   N2 = (r2 + g21*N1) / a22          [positive slope: g21/a22]
  #   hits N2-axis (N1=0) at: N2 = r2/a22
  #   hits N1-axis (N2=0) at: N1 = -r2/g21  (negative if r2 > 0, so no crossing)
  #
  # Interior equilibrium exists when isoclines cross in positive quadrant:
  #   Solve: (a11*N1 - r1)/g12 = (r2 + g21*N1)/a22
  #   a22*(a11*N1 - r1) = g12*(r2 + g21*N1)
  #   N1* = (a22*r1 + g12*r2) / (a11*a22 - g12*g21)   [needs a11*a22 ≠ g12*g21]
  #   N2* = (g21*r1 + a11*r2) / (a11*a22 - g12*g21)
  #
  # Stability: det(J) > 0 and tr(J) < 0 at interior equilibrium
  #   tr(J) = -a11*N1* - a22*N2* < 0  always (both terms positive)
  #   det(J) = N1*N2*(a11*a22 - g12*g21)
  #   → STABLE   when a11*a22 > g12*g21
  #   → UNSTABLE when a11*a22 < g12*g21
  
  denom <- a11 * a22 - g12 * g21
  N1_eq <- (a22 * r1 + g12 * r2) / denom
  N2_eq <- (g21 * r1 + a11 * r2) / denom
  
  if (is.null(N_max)) {
    N_max <- max(N1_eq, N2_eq, abs(r1 / a11), abs(r2 / a22)) * 2.2
  }
  
  # Vector field
  vf <- vector_field(params, N_max, n_grid, vf_scale)
  
  # Trajectories
  traj_list <- lapply(ic_list, function(ic)
    run_mut(ic[1], ic[2], params, t = t_seq))
  
  # ── Base plot ──────────────────────────────────────────────────────────
  plot(NA,
       xlim = c(0, N_max), ylim = c(0, N_max),
       xlab = expression(N[1]), ylab = expression(N[2]),
       main = title, las = 1)
  
  # Vector field arrows
  arrows(
    x0 = vf$grid$N1, y0 = vf$grid$N2,
    x1 = vf$grid$N1 + vf$dN1,
    y1 = vf$grid$N2 + vf$dN2,
    length = 0.04,
    col    = adjustcolor("steelblue", alpha.f = 0.45),
    lwd    = 0.8
  )
  
  # N1-isocline: N2 = (a11/g12)*N1 - r1/g12
  #   abline(a = intercept, b = slope) where intercept is at N1=0
  abline(a = -r1 / g12,
         b =  a11 / g12,
         col = "tomato", lwd = 2, lty = 2)
  
  # N2-isocline: N2 = (g21/a22)*N1 + r2/a22
  abline(a = r2 / a22,
         b = g21 / a22,
         col = "forestgreen", lwd = 2, lty = 2)
  
  # Mark interior equilibrium if it falls in the positive quadrant
  if (!is.na(N1_eq) && N1_eq > 0 && N2_eq > 0 &&
      N1_eq < N_max && N2_eq < N_max) {
    eq_stable <- (a11 * a22 > g12 * g21)
    points(N1_eq, N2_eq,
           pch = if (eq_stable) 16 else 21,
           cex = 1.8, lwd = 2,
           col = "black",
           bg  = if (eq_stable) "black" else "white")
  }
  
  # Trajectories
  traj_cols <- c("black", "purple", "darkorange", "navy")
  for (j in seq_along(traj_list)) {
    df <- traj_list[[j]]
    if (is.null(df)) next
    # Clip to plot region to avoid runaway lines
    keep <- df$N1 <= N_max * 1.05 & df$N2 <= N_max * 1.05
    df   <- df[keep, ]
    if (nrow(df) < 2) next
    col_j <- traj_cols[((j - 1) %% length(traj_cols)) + 1]
    lines(df$N2 ~ df$N1, col = col_j, lwd = 1.8)
    points(df$N1[1], df$N2[1], pch = 16, cex = 1.2, col = col_j)
  }
  
  # Legend
  eq_label <- if (a11 * a22 > g12 * g21) "Stable eq. (●)" else "Unstable eq. (○)"
  legend("topright",
         legend = c("Trajectory",
                    expression(N[1] ~ "isocline (dN1/dt=0)"),
                    expression(N[2] ~ "isocline (dN2/dt=0)"),
                    eq_label),
         col    = c("black", "tomato", "forestgreen", "black"),
         lty    = c(1, 2, 2, NA),
         pch    = c(NA, NA, NA, if (a11 * a22 > g12 * g21) 16 else 21),
         lwd    = 2,
         bg     = "white",
         cex    = 0.80)
}

# ══════════════════════════════════════════════════════════════════════════════
# Parameter sets for the two qualitative outcomes
# ══════════════════════════════════════════════════════════════════════════════

# ── Condition summary ────────────────────────────────────────────────────────
#
# Interior equilibrium (N1*, N2* > 0) exists when:
#   N1* = (a22*r1 + g12*r2) / (a11*a22 - g12*g21)  > 0
#   N2* = (g21*r1 + a11*r2) / (a11*a22 - g12*g21)  > 0
#
# Jacobian at (N1*, N2*):
#   J = [ -a11*N1*    g12*N1* ]
#       [  g21*N2*   -a22*N2* ]
#
#   tr(J)  = -a11*N1* - a22*N2*  < 0  always (all params positive)
#   det(J) = N1*N2* * (a11*a22 - g12*g21)
#
# (1) STABLE equilibrium:   a11*a22 > g12*g21   → det(J) > 0 → stable node/spiral
# (2) UNSTABLE equilibrium: a11*a22 < g12*g21   → det(J) < 0 → saddle point
#     trajectories near the equilibrium are repelled; populations either
#     collapse to 0 (if mutualism is insufficient) or explode to +∞
#     (runaway mutualism — unbounded growth)
#
# ─────────────────────────────────────────────────────────────────────────────

# ── (1) Stable positive equilibrium ─────────────────────────────────────────
# Choose: a11=0.4, a22=0.4, g12=0.1, g21=0.1
#   a11*a22 = 0.16  >  g12*g21 = 0.01  ✓  stable
#   r1=1, r2=1 (both species viable alone)
#   N1* = (0.4*1 + 0.1*1) / (0.16 - 0.01) = 0.5/0.15 ≈ 3.33
#   N2* = (0.1*1 + 0.4*1) / 0.15          = 0.5/0.15 ≈ 3.33
params_stable <- c(r1=1,   r2=1,
                   a11=0.40, a22=0.40,
                   g12=0.10, g21=0.10)
ic_stable <- list(c(0.5, 5.0), c(5.0, 0.5), c(0.5, 0.5), c(6.0, 6.0))

# ── (2) Unstable positive equilibrium ───────────────────────────────────────
# Choose: a11=0.2, a22=0.2, g12=0.4, g21=0.4
#   a11*a22 = 0.04  <  g12*g21 = 0.16  ✓  unstable (saddle)
#   r1=1, r2=1
#   denom = 0.04 - 0.16 = -0.12  (negative!)
#   N1* = (0.2*1 + 0.4*1) / (-0.12) = 0.6/(-0.12) = -5  ← NEGATIVE, not valid
#
# When denom < 0 and r1,r2 > 0 the equilibrium lands in the negative quadrant.
# To get a positive equilibrium with a11*a22 < g12*g21 we need
# numerators also negative: requires r1 < 0 or r2 < 0 (obligate mutualists —
# species cannot survive without the partner).
#
# Use r1 = -0.5, r2 = -0.5 (negative intrinsic growth — obligate mutualists):
#   denom = 0.04 - 0.16 = -0.12
#   N1* = (0.2*(-0.5) + 0.4*(-0.5)) / (-0.12) = -0.3 / -0.12 = 2.5  ✓
#   N2* = (0.4*(-0.5) + 0.2*(-0.5)) / (-0.12) = -0.3 / -0.12 = 2.5  ✓
#
# This is the biologically meaningful "obligate mutualism" scenario:
# each species declines alone (r < 0) but can persist WITH the partner.
# The interior equilibrium is a saddle: trajectories starting ABOVE the
# separatrix grow without bound (runaway mutualism); those starting BELOW
# collapse to extinction.
params_unstable <- c(r1=-0.5, r2=-0.5,
                     a11=0.20, a22=0.20,
                     g12=0.40, g21=0.40)
ic_unstable <- list(c(0.5, 4.0), c(4.0, 0.5),   # below separatrix → extinction
                    c(4.1, 1.5), c(2.0, 4.5))    # above separatrix → explosion

# ══════════════════════════════════════════════════════════════════════════════
# Produce the two phase-plane plots side by side
# ══════════════════════════════════════════════════════════════════════════════
par(mfrow = c(1, 2), mar = c(4, 4, 3.5, 1))

plot_scenario(
  params   = params_stable,
  ic_list  = ic_stable,
  title    = "(1) Stable positive equilibrium",
  N_max    = 8,
  vf_scale = 0.10,
  t_seq    = seq(0, 200, by = 0.2)
)

plot_scenario(
  params   = params_unstable,
  ic_list  = ic_unstable,
  title    = "(2) Unstable positive equilibrium\n(obligate mutualism)",
  N_max    = 8,
  vf_scale = 0.10,
  t_seq    = seq(0, 30, by = 0.05)
)

par(mfrow = c(1, 1))
