library(rdist)
library(SLHD)
library(magick)
library(lhs)
library(MASS)
import::from(base_functions.R, two_D_func, updated_moments, run_emulation, get_points,
             implausibility, PI, EI, EI2, EI3)
import::from("Chapter 4 - Extension to higher dimensions.R", run_history_match, plot_wave_points,
             variance_sampling, run_next_wave, get_non_implausible_ind)
import::from("Chapter 5 - Bayesian optimization.R", curr_optima, run_implausibility_optimization,
             run_myopic_optimization, update_new_points, scale)

borehole <- function(points, dim=8, scaling=list(lb=rep(0, dim), ub=rep(1, dim)),
                     set_values=rep(-1, 8)){
  points <- scale(points, dim=dim, scaling=scaling, up=TRUE)
  r_w <- if (set_values[1] == -1) points$r_w else set_values[1]
  r <- if (set_values[2] == -1) points$r else set_values[2]
  T_u <- if (set_values[3] == -1) points$T_u else set_values[3]
  H_u <- if (set_values[4] == -1) points$H_u else set_values[4]
  T_l <- if (set_values[5] == -1) points$T_l else set_values[5]
  H_l <- if (set_values[6] == -1) points$H_l else set_values[6]
  L <- if (set_values[7] == -1) points$L else set_values[7]
  K_w <- if (set_values[8] == -1) points$K_w else set_values[8]
  numerator <- 2 * pi * T_u * (H_u - H_l)
  denominator <- (log(r) - log(r_w)) * (1 + (2 * L * T_u) / ((log(r) - log(r_w)) * r_w^2 * K_w) + T_u/T_l)
  return (numerator / denominator)
}

nir_area <- function(new_points, lb, ub, dim=8){
  domain_area <- prod(ub - lb)
  nir_area <- 1
  for (i in 1:dim){
    nir_area <- nir_area * (max(new_points[,i]) - min(new_points[,i]))
  }
  return (100 * nir_area / domain_area)
}

contour_plots <- function(param_pairs, param_names, lb, ub, type="Initial", lb_new=rep(0, 8), 
                          ub_new=rep(0, 8), optimal_point=NULL){
  save_path <- "Plots/Chapter 6 - Optimizing the borehole function/"
  save_path <- if (type == "Initial") paste(save_path, "6.2 - Setup/", sep="") else paste(save_path, "6.3 - Results/", sep="")
  if (!dir.exists(save_path)) dir.create(save_path)
  for (i in 1:4){
    x <- if (type == "Initial") seq(0, 1, length.out=50) else seq(lb_new[param_pairs[i, 1]], ub_new[param_pairs[i, 1]], length.out=50) 
    y <- if (type == "Initial") seq(0, 1, length.out=50) else seq(lb_new[param_pairs[i, 2]], ub_new[param_pairs[i, 2]], length.out=50)
    points <- expand.grid(x, y)
    set_values <- rep(-1, 8)
    for (j in 1:8){
      set_values[j] <- if (type == "Initial") mean(c(lb[j], ub[j])) else optimal_point[j]
    }
    set_values[param_pairs[i,]] <- -1
    names(points) <- param_names[param_pairs[i,]]
    f <- 0
    if (type == "Initial"){
      f <- borehole(points, dim=2, scaling=list(lb=lb[param_pairs[i,]], ub=ub[param_pairs[i,]]),
                    set_values=set_values)
    }
    else{
      init_points <- as.data.frame(get_points(50, c(lb_new[param_pairs[i, 1]], lb_new[param_pairs[i, 2]]), c(ub_new[param_pairs[i, 1]], ub_new[param_pairs[i, 2]])))
      names(init_points) <- param_names[param_pairs[i,]]
      f <- borehole(init_points, dim=2, scaling=list(lb=lb[param_pairs[i,]], ub=ub[param_pairs[i,]]), set_values=set_values)
    }
    z <- if (type == "Initial") matrix(f,  nrow=50) else matrix(updated_moments(init_points, points, f, 350, 5, 0.3, dim=2)$Exp, nrow=50)
    img_path <- paste(save_path, param_names[param_pairs[i, 1]], " ", param_names[param_pairs[i, 2]], 
                      if(type != "Initial") " Exp" else "", ".png", sep="")
    png(img_path, width=504, height=504)
    filled.contour(x=x, y=y, z=z, xlab=param_names[param_pairs[i, 1]], 
                   ylab=param_names[param_pairs[i, 2]], color.palette=terrain.colors, cex.lab=1.3,
                   cex.axis=1.3)
    dev.off()
    image_write(image_crop(image_read(img_path), "504x460+0+44"), 
                path=img_path, format="png")
  }
}

wave_points_plots <- function(init_points, param_pairs, param_names, lb, ub){
  cols <- c(rep("yellow", 50), rep("orange", 50), rep("red", 50), rep("blue", 50), rep("black", 50))
  points <- scale(init_points[151:400,1:8], dim=8, scaling=list(lb=lb, ub=ub), up=FALSE)
  for (i in 1:4){
    pair_points <- points[, param_pairs[i,]]
    save_path <- paste("Plots/Chapter 6 - Optimizing the Borehole function/6.3 - Results/",
                       param_names[param_pairs[i, 1]], " ", param_names[param_pairs[i, 2]],
                       " Waves", ".pdf", sep="")
    pdf(save_path)
    plot(pair_points[,1], pair_points[,2], xlab=param_names[param_pairs[i, 1]], 
         ylab=param_names[param_pairs[i, 2]], cex.lab=1.3, cex.axis=1.3, col=cols, pch=19)
    dev.off()
  }
}

## 6.2 - Setup
## 6.2.1 - Initial runs and emulator points
lb <- c(0.05, 100, 63070, 990, 63.1, 700, 1120, 1500)
ub <- c(0.15, 50000, 115600, 1110, 116, 820, 1680, 15000)
param_names <- c("r_w", "r", "T_u", "H_u", "T_l", "H_l", "L", "K_w")
param_pairs <- matrix(c(1, 7,   2, 8,   3, 5,   4, 6), nrow=4, byrow=TRUE)
contour_plots(param_pairs, param_names, lb, ub)

init_points <- as.data.frame(get_points(150, rep(0, 8), rep(1, 8), dim=8))
names(init_points) <- param_names
D <- borehole(init_points, scaling=list(lb=lb, ub=ub))
pdf("Plots/Chapter 6 - Optimizing the Borehole function/6.2 - Setup/Pairs.pdf", width=8, height=6)
pairs(init_points, pch=19, col="blue", cex.axis=1.3, cex.lab=1.3)
dev.off()

## 6.2.2 - Hyperparameter values and emulator diagnostics
beta_0 <- 100
sigma_u <- 60
theta <- 0.8

diagnostic_runs <- as.data.frame(get_points(20, rep(0, 8), rep(1, 8), dim=8, maximin=TRUE))
names(diagnostic_runs) <- param_names
diagnostic_runs_f_x <- borehole(diagnostic_runs, scaling=list(lb=lb, ub=ub))
diagnostic_runs_moments <- updated_moments(scale(init_points, dim=8, 
                                                 scaling=list(lb=lb, ub=ub), up=FALSE),
                                           diagnostic_runs, D, beta_0, sigma_u, theta, dim=8)
prediction_errors <- abs(diagnostic_runs_f_x - diagnostic_runs_moments$Exp) / sqrt(diagnostic_runs_moments$Var)
prediction_errors
credible_intervals <- diagnostic_runs_moments$Exp + 3 * sqrt(diagnostic_runs_moments$Var) & diagnostic_runs_f_x > diagnostic_runs_moments$Exp - 3 * sqrt(diagnostic_runs_moments$Var)
credible_intervals

## 6.3 - Results
optimization <- run_implausibility_optimization(borehole, init_points, D, 10000, 
                                                c(beta_0, sigma_u, theta), lb, ub,
                                                rep(50, 5), 1e-1, optima=NULL, 
                                                scaled=TRUE, dim=8, param_names=param_names,
                                                new_points=new_points)
new_points <- optimization$new_points
new_points$Exp <- optimization$new_moments$Exp
new_points$Var <- optimization$new_moments$Var
new_points$I <- optimization$I
write.csv(new_points, "Borehole new points.csv")
init_points <- optimization$init_points
init_points$f_x <- optimization$D
write.csv(init_points, "Borehole init points.csv")
optimal_point <- optimization$optimal_point$Point
optimal_point$f_x <- optimization$optimal_point$Value
write.csv(optimal_point, "Borehole optimal point.csv")

lb_new <- rep(0, 8)
ub_new <- rep(1, 8)
for (i in 1:8){
  print(param_names[i])
  new_range <- range(new_points[,i])
  lb_new[i] <- (new_range[1] - lb[i]) / (ub[i] - lb[i])
  ub_new[i] <- (new_range[2] - lb[i]) / (ub[i] - lb[i])
  print("New range:")
  print(new_range)
  print("% of original range:")
  print(100 * (new_range[2] - new_range[1]) / (ub[i] - lb[i]))
  print("--------")
}

final_nir_area <- nir_area(new_points, lb, ub)
contour_plots(param_pairs, param_names, lb, ub, type='Final', lb_new=lb_new, ub_new=ub_new,
              optimal_point=as.matrix(optimal_point))

pdf("Plots/Chapter 6 - Optimizing the Borehole function/6.3 - Results/f(x) vs run.pdf")
plot(1:405, init_points$f_x, type="l", xlab="Run", ylab="f(x)", cex.lab=1.3, cex.axis=1.3)
abline(v=c(150, 200, 250, 300, 350, 400), lty=2)
dev.off()

wave_points_plots(init_points, param_pairs, param_names, lb, ub)

T_u <- data.frame(T_u=seq(0, 1, length.out=500))
set_values <- c(0.15, 100, -1, 1110, 116, 700, 1120, 15000)
borehole_T_u <- borehole(T_u, dim=1, scaling=list(lb=63070, ub=115600), set_values=set_values)
pdf("Plots/Chapter 6 - Optimizing the Borehole function/6.3 - Results/T_u.pdf")
plot(seq(63070, 115600, length.out=500), borehole_T_u, type="l", xlab="T_u", ylab="f(x)",
     cex.axis=1.3, cex.lab=1.3)
dev.off()

r <- data.frame(r=seq(0, 1, length.out=500))
set_values[2] <- -1
set_values[3] <- 115600
borehole_r <- borehole(r, dim=1, scaling=list(lb=100, ub=50000), set_values=set_values)
pdf("Plots/Chapter 6 - Optimizing the Borehole function/6.3 - Results/r.pdf")
plot(seq(100, 50000, length.out=500), borehole_r, type="l", xlab="r", ylab="f(x)",
     cex.axis=1.3, cex.lab=1.3)
dev.off()