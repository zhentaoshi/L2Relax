
# require(CVXR)
# require(Rmosek)
# require(SparseM)

# CVX L2-relaxation primal problem
rL2_primal <- function(Sigma, tau = 0) {
  
  N <- nrow(Sigma)
  w_gamma <- Variable(N + 1)
  w <- w_gamma[1:N]
  gamm <- w_gamma[N + 1]

  objective <- Minimize(0.5 * sum_squares(w))
  constraints <- list(
    sum(w) == 1,
    Sigma %*% w + gamm <= tau,
    -Sigma %*% w - gamm <= tau
  )
  
  problem <- Problem(objective, constraints)
  result <- solve(problem, solver = "ECOS_BB")
  w_hat <- result$getValue(w_gamma)[1:N]

  return(w_hat)
}


rL2_primal_mosek <- function(Sigma, tau, tol = 1e-7) {
  N <- nrow(Sigma)
  # variable order: w_1, w_2, ..., w_N, gamma, t, s, r

  prob <- list(sense = "min")
  prob$dparam <- list(INTPNT_CO_TOL_REL_GAP = tol)

  prob$c <- c(rep(0, N + 1), 1 / 2, rep(0, 2))
  
  A_1 <- rbind(
    c(rep(1, N), 0), # sum of weight == 1
    cbind(Sigma, rep(1, N)) # ||Sigma_hat w + gamma||_\infty \leq tau
  )
  A_2 <- rbind(c(1 / 2, -1, 0), c(1 / 2, 0, -1)) # transformation of the squared l2 norm
  A <- Matrix::bdiag(A_1, A_2)
  prob$A <- as(A, "CsparseMatrix")
  prob$bc <- rbind(
    blc = c(1, tau * rep(1, N), 1 / 2, -1 / 2),
    buc = c(1, -tau * rep(1, N), 1 / 2, -1 / 2)
  )
  prob$bx <- rbind(
    blx = c(rep(-Inf, N + 1), 0, rep(-Inf, 2)),
    bux = rep(Inf, N + 4)
  )
  # conic constraint
  prob$cones <- matrix(list("QUAD", c(N + 4, 1:N, N + 3)))
  rownames(prob$cones) <- c("type", "sub")

  mosek_out <- Rmosek::mosek(prob, opts = list(verbose = 0))

  xx <- mosek_out$sol$itr$xx
  w_hat <- xx[1:N]
  gamma_hat <- xx[N + 1]
  status <- mosek_out$sol$itr$solsta

  return(list(
    w = w_hat,
    gamma = gamma_hat,
    status = status
  ))
}

# CVX L2-relaxation dual problem

rL2_dual <- function(Sigma, tau = 0) {
  N <- nrow(Sigma)
  A_hat <- (diag(N) - matrix(1, N, N) / N) %*% Sigma

  alpha <- Variable(N)
  objective <- Minimize(0.5 * sum_squares(A_hat %*% alpha)
    + sum(Sigma %*% alpha) / N
    + tau * cvxr_norm(alpha, p = 1))
  constraints <- list(sum(alpha) == 0)
  problem <- Problem(objective, constraints)

  result <- solve(problem, solver = "ECOS_BB")

  alpha_hat <- result$getValue(alpha)
  w_hat <- A_hat %*% alpha_hat + 1 / N

  return(list(alpha = alpha_hat, w = w_hat))
}

rL2_dual_mosek <- function(Sigma, tau = 0, tol = 1e-6) {
  N <- nrow(Sigma)
  A_hat <- (diag(N) - matrix(1, N, N) / N) %*% Sigma
  # variable order: alpha, alpha_tilde, alpha_+, alpha_-, t, s, r

  prob <- list(sense = "min")
  prob$dparam <- list(INTPNT_CO_TOL_REL_GAP = tol)

  prob$c <- c(
    rep(0, 2 * N),
    c(rep(1 / N, N) %*% Sigma + tau),
    c(rep(-1 / N, N) %*% Sigma + tau),
    1 / 2, 0, 0
  )
  A_1 <- rbind(
    c(rep(1, N), rep(0, 3 * N)),
    cbind(A_hat, -diag(N), Matrix::spMatrix(N, 2 * N)),
    cbind(diag(N), Matrix::spMatrix(N, N), -diag(N), diag(N))
  )
  A_2 <- rbind(c(1 / 2, -1, 0), c(1 / 2, 0, -1))
  A <- Matrix::bdiag(A_1, A_2)
  prob$A <- as(A, "CsparseMatrix")

  prob$bc <- rbind(
    blc = c(rep(0, 1 + 2 * N), 1 / 2, -1 / 2),
    buc = c(rep(0, 1 + 2 * N), 1 / 2, -1 / 2)
  )
  prob$bx <- rbind(
    blx = c(rep(-Inf, 2 * N), rep(0, 2 * N + 1), rep(-Inf, 2)),
    bux = rep(Inf, 4 * N + 3)
  )
  # Conic constraint
  prob$cones <- matrix(list("QUAD", c(4 * N + 3, (N + 1):(2 * N), 4 * N + 2)))
  rownames(prob$cones) <- c("type", "sub")
  
  mosek_out <- Rmosek::mosek(prob, opts = list(verbose = 0))
  
  xx <- mosek_out$sol$itr$xx
  alpha_hat <- xx[1:N]
  status <- mosek_out$sol$itr$solsta
  w_hat <- A_hat %*% alpha_hat + 1 / N
  
  return(list(
    w = w_hat,
    alpha = alpha_hat,
    status = status
  ))
  
}

# constrained OLS regression
conOLS <- function(X, Y) {
  N <- nrow(X)

  b <- Variable(N)
  objective <- Minimize(sum_squares(Y - t(b) %*% X))
  constraints <- list(sum(b) == 1)
  problem <- Problem(objective, constraints)

  result <- solve(problem, solver = "ECOS_BB")

  b_hat <- result$getValue(b)

  return(b_hat)
}