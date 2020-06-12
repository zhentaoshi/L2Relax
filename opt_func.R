library(CVXR)




# CVX L2-relaxation primal problem

rL2_primal <- function(Sigma, tau = 0) {
  N <- nrow(Sigma)
  w_gamma <- Variable(N + 1)
  w <- w_gamma(1:N)
  gamm <- w_gamma(N + 1)

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
