rm( list = ls() )
source("opt_func.R")

# dimension

set.seed(808)

TT = 20
N = 10

# prepare the raw data

e = matrix( rnorm(N*TT), TT)
Sigma = crossprod(e)/TT



###### compare the optimizers
###### When Sigma is of full rank
######

w_hat_L2 = rL2_primal(Sigma, tau = 0)
w_hat_dual = rL2_dual(Sigma, tau = 0)$w

w_hat_L2_mosek <- rL2_primal_mosek(Sigma, tau = 0)$w
w_hat_dual_mosek = rL2_dual_mosek(Sigma, tau = 0)$w

#######################################
########### Group structure
#######################################
N_K1 = round(N*0.4); N_K2 = N - N_K1;

Sigma_co = matrix( c(1, 0.3, 0.3, 2), 2, 2 )
Sigma_star = rbind( cbind( matrix(Sigma_co[1,1], N_K1, N_K1), matrix(Sigma_co[1,2], N_K1, N_K2)  ) , 
                    cbind( matrix(Sigma_co[2,1], N_K2, N_K1), matrix(Sigma_co[2,2], N_K2, N_K2)  ) )


w_hat_L2 = rL2_primal(Sigma_star, tau = 0)
w_hat_dual = rL2_dual(Sigma_star, tau = 0)$w


# alpha_hat is not equal in group
print( rL2_dual(Sigma_star, tau = 0)$alpha)
print( rL2_dual(Sigma_star, tau = 0.01)$alpha)

########## # verify an equation in the paper
K = 2
ff = matrix( c(N_K1,N_K2)/N, ncol = 1 )
FF = diag(ff)
ones_K = matrix(1, K, 1)
A_co = diag( as.vector( sqrt(ff) ) ) %*%( diag(K) - ones_K %*% t(ff) ) %*% Sigma_co 


print( t(A_co) %*% sqrt(ff) ) # verify an equation in the paper








