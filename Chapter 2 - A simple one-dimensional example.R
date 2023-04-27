library(rdist)
library(SLHD)
library(magick)
library(lhs)
library(MASS)
import::from(base_functions.R, one_D_func, updated_moments, run_emulation)

## 2.1 - Setup
n <- 5
init_points <- seq(0, 3, length.out=n)
D <- one_D_func(init_points)
beta_0 <- 0
sigma_u <- 0.9
theta <- 0.4

## 2.2 - Emulation of f(x)
new_points <- seq(from=0, to=3, length.out=300) + 0.001 # small value to avoid computational errors
save_path <- "Plots/Chapter 2 - A simple one-dimensional example/2.2 - Emulation of f(x)/"
new_moments <- run_emulation(init_points, new_points, D, beta_0, sigma_u, theta, true_func=TRUE,
                             func=one_D_func, xlim=c(0, 3), ylim=c(-3, 3), add_run_points=TRUE,
                             save_path=save_path)

## 2.3 - Investigating the hyperparameters
### 2.3.1 - sigma_u
sigma_u_vals <- c(0.45, 0.9, 1.8)
sigma_u_cols <- c("black", "red", "green")
pdf("Plots/Chapter 2 - A simple one-dimensional example/2.3 - Investigating the hyperparameters/sigma_u.pdf")
plot(init_points, D, pch=19, col="purple", xlim=c(0, 3), ylim=c(-4, 4), xlab="x", ylab="f(x)", 
     cex.lab=1.3, cex.axis=1.3)
for (i in 1:length(sigma_u_vals)){
  new_moments <- run_emulation(init_points, new_points, D, beta_0, sigma_u_vals[i], theta, plots=FALSE)
  E_fx_new <- new_moments$Exp
  Var_fx_new <- new_moments$Var
  if (i == 1){
    lines(new_points, E_fx_new, col="blue")
  }
  lines(new_points, E_fx_new + 3 * sqrt(Var_fx_new), col=sigma_u_cols[i])
  lines(new_points, E_fx_new - 3 * sqrt(Var_fx_new), col=sigma_u_cols[i])
}
dev.off()

### 2.3.2 - theta
theta_vals <- c(0.2, 0.8)
save_path <- "Plots/Chapter 2 - A simple one-dimensional example/2.3 - Investigating the hyperparameters/theta/"
for (i in 1:length(theta_vals)){
  new_moments <- run_emulation(init_points, new_points, D, beta_0, sigma_u, theta_vals[i],
                               add_run_points=TRUE, xlim=c(0, 3), ylim=c(-3, 3),
                               save_path=save_path, i=i)
}

### 2.3.3 - beta_0
new_points <- seq(-0.1, 3.1, length.out=320) + 0.001
save_path <- "Plots/Chapter 2 - A simple one-dimensional example/2.3 - Investigating the hyperparameters/beta_0/"
beta_0_vals <- c(-1, 1)
for (i in 1:length(beta_0_vals)){
  new_moments <- run_emulation(init_points, new_points, D, beta_0_vals[i], sigma_u, theta, 
                               add_run_points=TRUE, xlim=c(-0.1, 3.1), ylim=c(-3, 3),
                               save_path=save_path, i=i)
}