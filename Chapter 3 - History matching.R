library(rdist)
library(SLHD)
library(magick)
library(lhs)
library(MASS)
import::from(base_functions.R, one_D_func, updated_moments, run_emulation, implausibility)

plot_z <- function(z, epsilon, err, new_points, E_fx_new, xlim=NULL, ylim=NULL, save_path="", i=""){
  if (save_path != "") if (!dir.exists(save_path)) dir.create(save_path)
  if (save_path != "") pdf(paste(save_path, "Measurement", i, ".pdf", sep=""))
  xlim <- if (is.null(xlim)) c(min(new_points), max(new_points)) else xlim
  ylim <- if (is.null(ylim)) c(min(E_fx_new), max(E_fx_new)) else ylim
  plot(new_points, E_fx_new, col="blue", type="l", xlab="x", ylab="f(x)", cex.lab=1.3, cex.axis=1.3, 
       xlim=xlim, ylim=ylim)
  abline(h=z)
  abline(h=c(z - 3*sqrt(err), z + 3*sqrt(err)), lty=2)
  lines(new_points, E_fx_new + 3 * sqrt(epsilon), lty=2, col="blue")
  lines(new_points, E_fx_new - 3 * sqrt(epsilon), lty=2, col="blue")
  if (save_path != "") dev.off()
}

plot_implausibility <- function(z, epsilon, err, new_points, new_moments, save_path="", i=""){
  if (save_path != "") if (!dir.exists(save_path)) dir.create(save_path)
  if (save_path != "") pdf(paste(save_path, "Implausibility", i, ".pdf", sep=""))
  E_fx_new <- new_moments$Exp
  Var_fx_new <- new_moments$Var
  I <- implausibility(z, E_fx_new, Var_fx_new, epsilon, err)
  plot(new_points, I, type="l", ylim=c(0, 5), xlab="x", ylab="I(x)", cex.lab=1.3, cex.axis=1.3)
  abline(h=3, col="red")
  if (save_path != "") dev.off()
}

## Figure 3.1
n <- 5
beta_0 <- 0
sigma_u <- 0.9
theta <- 0.4
init_points <- seq(0, 3, length.out=n)
D <- one_D_func(init_points)
new_points <- seq(from=0, to=3, length.out=300) + 0.001
z <- 1.5
epsilon <- 0.005
err <- 0.001
new_moments <- run_emulation(init_points, new_points, D, beta_0, sigma_u, theta, plots=FALSE)
save_path <- "Plots/Chapter 3 - History matching/"
plot_z(z, epsilon, err, new_points, new_moments$Exp, xlim=c(0, 3), ylim=c(-2.5, 2.5),
       save_path=save_path)

## 3.1 - Implausibility
save_path <- "Plots/Chapter 3 - History matching/3.1 - Implausibility/"
plot_implausibility(z, epsilon, err, new_points, new_moments, save_path=save_path)

## 3.2 - Running further waves
save_path <- "Plots/Chapter 3 - History matching/3.2 - Running further waves/"
init_points <- append(init_points, c(0.4, 1.4, 1.6, 2.75, 2.85)) 
D <- append(D, one_D_func(init_points[6:10]))
new_moments <- run_emulation(init_points, new_points, D, beta_0, sigma_u, theta, add_run_points=TRUE, 
                             xlim=c(0, 3), ylim=c(-3, 3), save_path=save_path)
plot_implausibility(z, epsilon, err, new_points, new_moments, save_path=save_path)