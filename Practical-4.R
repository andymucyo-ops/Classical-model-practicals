### Lotka-Volterra model for 2 competing species
### dN1/dt = N1 * (r1 - alpha11*N1 - alpha12*N2)
### dN2/dt = N2 * (r2 - alpha21*N1 - alpha22*N2)
###
### Four qualitative outcomes:
###   (1) Coexistence
###   (2) Species 1 wins
###   (3) Species 2 wins
###   (4) Priority effect (founder control)

library(deSolve)

# ─────────────────────────────────────────────
# ODE function
# ─────────────────────────────────────────────
LV_competition <- function(time, state, params) {
  N1 <- state["N1"]
  N2 <- state["N2"]
  
  r1     <- params["r1"]
  r2     <- params["r2"]
  a11    <- params["a11"]
  a12    <- params["a12"]
  a21    <- params["a21"]
  a22    <- params["a22"]
  
  dN1 <- N1 * (r1 - a11 * N1 - a12 * N2)
  dN2 <- N2 * (r2 - a21 * N1 - a22 * N2)
  
  return(list(c(dN1, dN2)))
}

# ─────────────────────────────────────────────
# Helper: run ODE for a set of initial conditions
# ─────────────────────────────────────────────
run_lv <- function(N1_0, N2_0, params,
                   t = seq(0, 100, by = 0.5)) {
  pop  <- c(N1 = N1_0, N2 = N2_0)
  out  <- ode(pop, t, LV_competition, parms = params)
  as.data.frame(out)
}

# ─────────────────────────────────────────────
# Helper: vector field (normalised arrows)
# ─────────────────────────────────────────────
vector_field <- function(params, N_max, n_grid = 18, scale = 0.06) {
  N_seq <- seq(0.5, N_max, length.out = n_grid)
  grid  <- expand.grid(N1 = N_seq, N2 = N_seq)
  
  dN1_vec <- numeric(nrow(grid))
  dN2_vec <- numeric(nrow(grid))
  
  for (k in seq_len(nrow(grid))) {
    state_k <- c(N1 = grid$N1[k], N2 = grid$N2[k])
    derivs  <- LV_competition(0, state_k, params)[[1]]
    dN1_vec[k] <- derivs[1]
    dN2_vec[k] <- derivs[2]
  }
  
  mag       <- sqrt(dN1_vec^2 + dN2_vec^2)
  mag[mag == 0] <- 1   # avoid division by zero
  dN1_norm  <- scale * dN1_vec / mag
  dN2_norm  <- scale * dN2_vec / mag
  
  list(grid = grid, dN1 = dN1_norm, dN2 = dN2_norm)
}

# ─────────────────────────────────────────────
# Helper: plot one scenario
# ─────────────────────────────────────────────
plot_scenario <- function(params, ic_list, title,
                          N_max   = NULL,
                          n_grid  = 18,
                          vf_scale = 0.06) {
  
  r1  <- params["r1"];  r2  <- params["r2"]
  a11 <- params["a11"]; a12 <- params["a12"]
  a21 <- params["a21"]; a22 <- params["a22"]
  
  # Zero-growth isoclines (non-trivial)
  # N1-isocline: N2 = (r1 - a11*N1) / a12   → line in (N1, N2) space
  # N2-isocline: N2 = (r2 - a21*N1) / a22
  
  K1  <- r1 / a11   # N1 carrying capacity (N2 = 0)
  K2  <- r2 / a22   # N2 carrying capacity (N1 = 0)
  K12 <- r1 / a12   # N2 value where N1-isocline hits N2-axis
  K21 <- r2 / a21   # N1 value where N2-isocline hits N1-axis
  
  if (is.null(N_max)) N_max <- max(K1, K2, K12, K21) * 1.3
  
  # Vector field
  vf <- vector_field(params, N_max, n_grid, vf_scale)
  
  # Trajectories
  traj_list <- lapply(ic_list, function(ic)
    run_lv(ic[1], ic[2], params))
  
  # ── Plot ──
  plot(NA,
       xlim = c(0, N_max), ylim = c(0, N_max),
       xlab = expression(N[1]), ylab = expression(N[2]),
       main = title, las = 1)
  
  # Vector field
  arrows(
    x0 = vf$grid$N1, y0 = vf$grid$N2,
    x1 = vf$grid$N1 + vf$dN1,
    y1 = vf$grid$N2 + vf$dN2,
    length = 0.04,
    col    = adjustcolor("steelblue", alpha.f = 0.45),
    lwd    = 0.8
  )
  
  # N1-isocline: N2 = (r1 - a11*N1) / a12
  # passes through (K1, 0) and (0, K12)
  abline(a = K12, b = -a11 / a12, col = "tomato", lwd = 2, lty = 2)
  
  # N2-isocline: N2 = (r2 - a21*N1) / a22
  # passes through (K21, 0) and (0, K2)
  abline(a = K2, b = -a21 / a22, col = "forestgreen", lwd = 2, lty = 2)
  
  # Trajectories
  cols <- c("black", "purple", "darkorange", "navy")
  for (j in seq_along(traj_list)) {
    df <- traj_list[[j]]
    lines(df$N2 ~ df$N1,
          col = cols[((j - 1) %% length(cols)) + 1],
          lwd = 1.8)
    # starting point
    points(df$N1[1], df$N2[1], pch = 16, cex = 1.2,
           col = cols[((j - 1) %% length(cols)) + 1])
  }
  
  # Legend
  legend("topright",
         legend = c("Trajectory", expression(N[1] ~ "isocline"),
                    expression(N[2] ~ "isocline")),
         col    = c("black", "tomato", "forestgreen"),
         lty    = c(1, 2, 2),
         lwd    = 2,
         bg     = "white",
         cex    = 0.85)
}

# ══════════════════════════════════════════════
# Parameter sets for the four qualitative outcomes
# ══════════════════════════════════════════════

# ── Condition summary ──────────────────────────────────────────────────────
# Let K1 = r1/a11, K2 = r2/a22, K12 = r1/a12, K21 = r2/a21
#
# (1) Coexistence:       K21 > K1  AND  K12 > K2
#     i.e. a11 > a21  AND  a22 > a12
#     (intraspecific competition stronger than interspecific for each species)
#
# (2) Species 1 wins:    K1 > K21  AND  K2 < K12
#     (N1 isocline above N2 isocline everywhere)
#
# (3) Species 2 wins:    K1 < K21  AND  K2 > K12
#     (N2 isocline above N1 isocline everywhere)
#
# (4) Priority effect:   K1 < K21  AND  K2 < K12
#     (isoclines cross; winner depends on initial conditions)
# ─────────────────────────────────────────────────────────────────────────

# ── (1) Coexistence ─────────────────────────────────────────────────────
# Intraspecific >> interspecific competition ensures a wide basin of attraction
# The isocline intercepts give:
#   N1-isocline hits N1-axis at K1=r1/a11, hits N2-axis at K12=r1/a12
#   N2-isocline hits N1-axis at K21=r2/a21, hits N2-axis at K2=r2/a22
# Coexistence: N1-isocline ABOVE N2-isocline at N2-axis (K12 > K2)
#              AND N2-isocline ABOVE N1-isocline at N1-axis (K21 > K1)
#   → K12 > K2  AND  K21 > K1
#   → r1/a12 > r2/a22  AND  r2/a21 > r1/a11
#   → with r1=r2=1: 1/a12 > 1/a22 → a12 < a22  (inter < intra for N2 effect on N1)
#                   1/a21 > 1/a11  → a21 < a11  (inter < intra for N1 effect on N2)
# So coexistence: a11 > a21  AND  a22 > a12  ✓
#
# Chosen: a11=0.40, a22=0.40 (strong intra), a12=0.10, a21=0.10 (weak inter)
#   K1=2.5, K21=10, K2=2.5, K12=10
#   K21 > K1,  K12 > K2  → coexistence, interior equilibrium ≈ (2.27, 2.27)
params_coex <- c(r1=1, r2=1,
                 a11=0.40, a12=0.10,
                 a21=0.10, a22=0.40)
ic_coex <- list(c(0.5, 3.5), c(3.5, 0.5), c(8.0, 4.0), c(4.0, 8.0))


# ── (2) Species 1 wins ──────────────────────────────────────────────────
# K1=4, K21=2  ->  K1 > K21
# K2=3, K12=5  ->  K2 < K12
params_sp1 <- c(r1=1,    r2=0.75,
                a11=0.25, a12=0.20,
                a21=0.50, a22=0.25)
ic_sp1 <- list(c(0.5, 3.0), c(5.0, 4.0), c(0.5, 0.5), c(2.5, 2.5))

# ── (3) Species 2 wins ──────────────────────────────────────────────────
# K1=4, K21=6  ->  K1 < K21
# K2=3, K12=2  ->  K2 > K12
params_sp2 <- c(r1=1,    r2=0.75,
                a11=0.25, a12=0.50,
                a21=0.125, a22=0.25)
ic_sp2 <- list(c(3.0, 6.0), c(3.5, 0.5), c(0.5, 0.5), c(2.5, 2.5))

# ── (4) Priority effect ─────────────────────────────────────────────────
# K1=4, K21=6  ->  K1 < K21
# K2=4, K12=6  ->  K2 < K12
params_prio <- c(r1=1,  r2=1,
                 a11=0.25, a12=0.167,
                 a21=0.167, a22=0.25)
ic_prio <- list(c(0.3, 5.0), c(2.3, 0.3), c(0.5, 0.5), c(6.0, 4))

# ══════════════════════════════════════════════
# Produce the four phase-plane plots
# ══════════════════════════════════════════════
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

plot_scenario(params_coex,
              ic_list  = ic_coex,
              title    = "(1) Coexistence",
              vf_scale = 0.08)

plot_scenario(params_sp1,
              ic_list  = ic_sp1,
              title    = "(2) Species 1 wins",
              vf_scale = 0.08)

plot_scenario(params_sp2,
              ic_list  = ic_sp2,
              title    = "(3) Species 2 wins",
              vf_scale = 0.08)

plot_scenario(params_prio,
              ic_list  = ic_prio,
              title    = "(4) Priority effect (founder control)",
              vf_scale = 0.08)

par(mfrow = c(1, 1))  # reset layout