library(rdist)
library(SLHD)
library(magick)
library(lhs)
library(MASS)
import::from(base_functions.R, two_D_func, updated_moments, run_emulation, get_points,
             implausibility, PI, EI, EI2, EI3)
import::from("Chapter 4 - Extension to higher dimensions.R", run_history_match, plot_wave_points,
             variance_sampling, run_next_wave, get_non_implausible_ind)

optim_func_one_D <- function(points){
  return (exp(points/5) * sin(6*points))
}

two_peak_func <- function(points, scaling=NULL){
  f <- if (is.null(nrow(points))) cos(points[1]) + sin(points[2]) else cos(points[,1]) + sin(points[,2])
  return (f)
}

curr_optima <- function(D, init_points, dim=1){
  optimal_point <- list(Point=if (dim == 1) init_points[which.max(D)] else init_points[which.max(D),1:dim],
                        Value=max(D))
  return (optimal_point)
}

run_myopic_optimization <- function(acq_func, model_func, init_points, D, n_new_points, hyperparams, 
                                    lb, ub, xlim=NULL, ylim=NULL, zlim=NULL, dim=1, m=10, 
                                    plot_acq=FALSE, scaled=FALSE, new_points=NULL, param_names=NULL,
                                    colors=NULL, save_path="", ...){
  fx_top <- max(D)
  n <- length(D)
  beta_0 <- hyperparams[1]
  sigma_u <- hyperparams[2]
  theta <- hyperparams[3]
  
  colors <- if (is.null(colors)) rep("blue", n + m) else append(colors, rep("blue", m))
  if (is.null(new_points)){
    new_points <- seq(lb[1], ub[1], length.out=n_new_points)
    if (dim == 2){
      x <- if (scaled) seq(0, 1, length.out=sqrt(n_new_points)) else seq(lb[1], ub[1], length.out=sqrt(n_new_points))
      y <- if (scaled) seq(0, 1, length.out=sqrt(n_new_points)) else seq(lb[2], ub[2], length.out=sqrt(n_new_points))
      new_points <- as.matrix(expand.grid(x, y))
    }
    else if (dim > 2){
      new_points <- if (scaled) as.data.frame(get_points(n_new_points, rep(0, dim), rep(1, dim), dim=dim, maximin=FALSE)) else as.data.frame(get_points(n_new_points, lb, ub, dim=dim, maximin=FALSE))
      if (!is.null(param_names)) names(new_points) <- param_names
    }
  }
  for (i in 1:m){
    new_moments <- run_emulation(init_points, new_points, D, beta_0, sigma_u, theta, xlim=xlim, ylim=ylim,
                                 zlim=zlim, save_path=save_path, i=i, add_run_points=TRUE, 
                                 run_cols=colors, true_func=TRUE, func=model_func, dim=dim)
    E_fx_new <- new_moments$Exp
    Var_fx_new <- new_moments$Var
    
    acq <- acq_func(E_fx_new, Var_fx_new, fx_top, ...)
    next_point <- if (dim == 1) new_points[which.max(acq)] else (if (dim == 2) t(as.matrix(new_points[which.max(acq),])) else new_points[which.max(acq), 1:dim]) 
    D <- if (dim == 1) append(D, model_func(next_point)) else append(D, model_func(next_point, scaling=list(lb=lb, ub=ub))) 
    init_points <- if (dim == 1) append(init_points, next_point) else rbind(init_points, next_point)
    fx_top <- max(D)
    
    colors[length(D)] <- "black"
    if (i > 1){
      colors[length(D) - 1] <- "orange"
    }
    
    if (plot_acq){
      if (dim == 1){
        if (save_path != "") pdf(paste(save_path, "Acquisition", i, ".pdf", sep=""))
        plot(new_points, acq, type="l", ylim=c(0, max(acq)))
        points(init_points, rep(0, length(D)), col=colors, pch=19)
        if (save_path != "") dev.off()
      }
      else if (dim == 2){
        img_path <- paste(save_path, "Acquisition Function", i, ".png", sep="")
        if (save_path != "") png(img_path, width=504, height=504)
        filled.contour(x=x, y=y, z=acq, xlab="x", ylab="y", color.palette=topo.colors,
                       plot.axes={axis(1);axis(2);points(init_points[,1], init_points[,2], pch=19, 
                                                         col=colors)},
                       cex.axis=1.3, cex.lab=1.3)
        if (save_path != "") dev.off()
        if (save_path != "") image_write(image_crop(image_read(img_path), "504x460+0+44"), 
                                         path=img_path, format="png")
      }
    }
  }
  new_moments <- run_emulation(init_points, new_points, D, beta_0, sigma_u, theta, xlim=xlim, ylim=ylim,
                               zlim=zlim, save_path=save_path, i=m+1, add_run_points=TRUE,
                               run_cols=colors, dim=dim)
  
  optimal_point <- curr_optima(D, init_points, dim=dim)
  optimization <- list(optimal_point, init_points, D, new_moments)
  names(optimization) <- c("optimal_point", "init_points", "D", "new_moments")
  return (optimization)
}

if (sys.nframe() == 0){
  ## 5.1 - Acquisition functions
  n <- 6
  beta_0 <- 0
  sigma_u <- 0.5
  theta <- 0.25
  lb <- 1
  ub <- 3
  init_points <- seq(1, 3, length.out=n)
  D <- optim_func_one_D(init_points)
  n_new_points <- 500
  
  ## 5.1.1 - Probability of Improvement
  save_path <- "Plots/Chapter 5 - Bayesian optimization/5.1 - Acquisition functions/PI/"
  optimal_point_pi <- run_myopic_optimization(PI, optim_func_one_D, init_points, D, n_new_points,
                                              c(beta_0, sigma_u, theta), lb, ub, xlim=c(lb, ub),
                                              ylim=c(-2, 2), m=10, save_path=save_path,
                                              plot_acq=TRUE)
  save_path <- "Plots/Chapter 5 - Bayesian optimization/5.1 - Acquisition functions/PI-gamma/"
  optimal_point_pi_gamma <- run_myopic_optimization(PI, optim_func_one_D, init_points, D, n_new_points, 
                                                    c(beta_0, sigma_u, theta), lb, ub, xlim=c(lb, ub),
                                                    ylim=c(-2, 2), m=10, save_path=save_path,
                                                    plot_acq=TRUE, gamma=0.01)
  
  ## 5.1.2 - Expected Improvement
  save_path <- "Plots/Chapter 5 - Bayesian optimization/5.1 - Acquisition functions/EI/"
  optimal_point_ei <- run_myopic_optimization(EI, optim_func_one_D, init_points, D, n_new_points,
                                              c(beta_0, sigma_u, theta), lb, ub, xlim=c(lb, ub),
                                              ylim=c(-2, 2), m=10, save_path=save_path, plot_acq=TRUE)
  
  ## 5.1.3 - Generalized Expected Improvement
  n <- 10
  beta_0 <- 0
  sigma_u <- 0.5
  theta <- 1.4
  lb <- c(-pi/4, -pi/4)
  ub <- c(9*pi/4, 9*pi/4)
  init_points <- get_points(n, lb, ub, hypercube=TRUE, maximin=TRUE)
  D <- two_peak_func(init_points)
  n_new_points <- 2500
  save_path <- "Plots/Chapter 5 - Bayesian optimization/5.1 - Acquisition functions/GEI/q=1/"
  optimal_point <- run_myopic_optimization(EI, two_peak_func, init_points, D, n_new_points, 
                                           c(beta_0, sigma_u, theta), lb, ub, zlim=c(-2, 2), dim=2, 
                                           m=10, plot_acq=TRUE, save_path=save_path)
  save_path <- "Plots/Chapter 5 - Bayesian optimization/5.1 - Acquisition functions/GEI/q=2/"
  optimal_point_ei2 <- run_myopic_optimization(EI2, two_peak_func, init_points, D, n_new_points, 
                                               c(beta_0, sigma_u, theta), lb, ub, zlim=c(-2, 2), dim=2, 
                                               m=10, plot_acq=TRUE, save_path=save_path)
  save_path <- "Plots/Chapter 5 - Bayesian optimization/5.1 - Acquisition functions/GEI/q=3/"
  optimal_point_ei3 <- run_myopic_optimization(EI3, two_peak_func, init_points, D, n_new_points, 
                                               c(beta_0, sigma_u, theta), lb, ub, zlim=c(-2, 2), dim=2, 
                                               m=10, plot_acq=TRUE, save_path=save_path)
}

scale <- function(points, dim=8, scaling=list(lb=rep(0, dim), ub=rep(1, dim)), up=TRUE){
  lb <- scaling$lb
  ub <- scaling$ub
  for (i in 1:dim){
    points[,i] <- if (up) lb[i] + points[,i] * (ub[i] - lb[i]) else (points[,i] - lb[i]) / (ub[i] - lb[i])
  }
  return (points)
}

update_new_points <- function(n_new_points, new_points, I, lb, ub, scaled=FALSE, 
                              param_names=NULL, dim=2){
  lb_new <- rep(0, dim)
  ub_new <- rep(0, dim)
  non_implausible <- get_non_implausible_ind(I)
  non_implausible_points <- if (dim == 2) new_points[sqrt(n_new_points) * (non_implausible[,2] - 1) + non_implausible[,1],] else new_points[non_implausible,1:dim]
  lb <- if (scaled) rep(0, dim) else lb
  ub <- if (scaled) rep(1, dim) else ub
  for (k in 1:dim){
    mn <- min(non_implausible_points[,k]) 
    mx <- max(non_implausible_points[,k])
    lb_new[k] <- if (mn - lb[k] / 20 < lb[k]) lb[k] else mn - lb[k] / 20
    ub_new[k] <- if (mx + ub[k] / 20 > ub[k]) ub[k] else mx + ub[k] / 20
  }
  if (dim == 2){
    return (get_points(n_new_points, lb_new, ub_new, hypercube=FALSE))
  }
  else{
    new_points <- as.data.frame(get_points(n_new_points, lb_new, ub_new, dim=dim, maximin=FALSE))
    if (!is.null(param_names)){
      names(new_points) <- param_names
    }
    return (new_points)
  }
}

run_implausibility_optimization <- function(model_func, init_points, D, n_new_points, hyperparams, lb,
                                            ub, waves_config, epsilon, optima, xlim=NULL, ylim=NULL, 
                                            zlim=NULL, dim=1, plot_acq=FALSE, scaled=FALSE,
                                            param_names=NULL, new_points=NULL, save_path="", ...){
  fx_top <- max(D)
  beta_0 <- hyperparams[1]
  sigma_u <- hyperparams[2]
  theta <- hyperparams[3]

  
  if (is.null(new_points)){
    if (dim == 2){
      x <- if (scaled) seq(0, 1, length.out=sqrt(n_new_points)) else seq(lb[1], ub[1], length.out=sqrt(n_new_points))
      y <- if (scaled) seq(0, 1, length.out=sqrt(n_new_points)) else seq(lb[2], ub[2], length.out=sqrt(n_new_points))
      new_points <- as.matrix(expand.grid(x, y))
    }
    else{
      new_points <- if (scaled) as.data.frame(get_points(n_new_points, rep(0, dim), rep(1, dim), dim=dim, maximin=FALSE)) else as.data.frame(get_points(n_new_points, lb, ub, dim=dim, maximin=FALSE))
      if (!is.null(param_names)) names(new_points) <- param_names
    }
  }
  colors <- rep("blue", nrow(init_points) + sum(waves_config))
  new_moments <- run_emulation(init_points, new_points, D, beta_0, sigma_u, theta, xlim=xlim, ylim=ylim,
                               zlim=zlim, save_path=save_path, 
                               run_cols=colors, true_func=TRUE, func=model_func, dim=dim)
  I <- run_history_match(fx_top, epsilon, 0, new_points, new_moments, optima=optima, optim=TRUE,
                         save_path=save_path, dim=dim)
  save_path_1 <- if (save_path != "") paste(save_path, "Waves/", sep="") else ""
  for (j in 1:length(waves_config)){
    new_points <- update_new_points(n_new_points, new_points, I, lb, ub, scaled=scaled, dim=dim,
                                    param_names=param_names)
    new_moments <- run_emulation(init_points, new_points, D, beta_0, sigma_u, theta, 
                                 xlim=xlim, ylim=ylim, i=j + 1,
                                 save_path=save_path_1, dim=dim)
    I <- run_history_match(fx_top, epsilon, 0, new_points, new_moments, optima=optima, optim=TRUE,
                           save_path=save_path_1, i=j+1, dim=dim)
    next_wave <- run_next_wave(waves_config[j], model_func, I, new_moments, init_points, new_points, 
                               D, beta_0, sigma_u, theta, fx_top, epsilon, 0, optima=optima,
                               optim=TRUE, save_path=save_path_1, i=j + 1, zlim=zlim, dim=dim, 
                               scaling=list(lb=lb, ub=ub))
    init_points <- next_wave$init_points
    D <- next_wave$D
    I <- next_wave$I
    fx_top <- max(D)
    print(paste("Wave ", j, " done", sep=""))
  }
  new_points <- update_new_points(n_new_points, new_points, I, lb, ub, scaled=scaled, dim=dim,
                                  param_names=param_names)
  save_path_2 <- if (save_path != "") paste(save_path, "EI/", sep="") else ""
  optimization <- run_myopic_optimization(EI, model_func, init_points, D, 1, hyperparams, lb, ub,
                                          xlim=xlim, ylim=ylim, dim=dim, m=4,
                                          plot_acq=plot_acq, scaled=scaled, new_points=new_points,
                                          colors=colors, param_names=param_names,
                                          save_path=save_path_2, ...)
  init_points <- optimization$init_points
  D <- optimization$D
  save_path_3 <- if (save_path != "") paste(save_path, "PI/", sep="") else ""
  optimization <- run_myopic_optimization(PI, model_func, init_points, D, 1, hyperparams, lb, ub,
                                          xlim=xlim, ylim=ylim, dim=dim, m=1,
                                          plot_acq=plot_acq, scaled=scaled, new_points=new_points,
                                          param_names=param_names, save_path=save_path_3, ...)
  init_points <- optimization$init_points
  D <- optimization$D
  fx_top <- optimization$optimal_point$Value
  I <- run_history_match(fx_top, epsilon, 0, new_points, optimization$new_moments,
                         optima=optima, optim=TRUE, plot=FALSE, dim=dim)
  new_points <- update_new_points(n_new_points, new_points, I, lb, ub, scaled=scaled,
                                  dim=dim, param_names=param_names)
  optimization$new_moments <- run_emulation(init_points, new_points, D, beta_0, sigma_u, theta, 
                                            add_run_points=FALSE, xlim=xlim, ylim=ylim, i="FINAL",
                                            save_path=save_path, dim=dim)
  optimization$I <- run_history_match(fx_top, epsilon, 0, new_points, optimization$new_moments,
                                      optima=optima, optim=TRUE, save_path=save_path, dim=dim,
                                      i="FINAL")
  optimization$optimal_point$Point <- if (scaled) scale(optimization$optimal_point$Point, dim=dim, scaling=list(lb=lb, ub=ub)) else optimization$optimal_point$Point 
  optimization$init_points <- if (scaled) scale(init_points, dim=dim, scaling=list(lb=lb, ub=ub)) else init_points
  optimization$new_points <- if (scaled) scale(new_points, dim=dim, scaling=list(lb=lb, ub=ub)) else new_points
  return (optimization)
}

run_wave_config_analysis <- function(waves_configs, funcs_hyperparams, funcs_bounds,
                                     funcs_maximums, funcs_maximums_loc, func_names){
  mean_vars <- matrix(0, nrow=nrow(waves_configs), ncol=length(func_names))
  final_nir_areas <- matrix(0, nrow=nrow(waves_configs), ncol=length(func_names))
  final_maximum_distances <- matrix(0, nrow=nrow(waves_configs), ncol=length(func_names))
  final_loc_distances <- matrix(0, nrow=nrow(waves_configs), ncol=length(func_names))
  n <- 25
  for (f in 1:length(func_names)){
    func <- if (f == 1) trid else (if (f == 2) sphere else sum_squares)
    hyperparameters <- c(funcs_hyperparams[f,])
    lb <- funcs_bounds[f,1:2]
    ub <- funcs_bounds[f,3:4]
    init_points <- get_points(n, lb, ub, hypercube=TRUE, maximin=TRUE)
    D <- func(init_points)
    maximum_loc <- c(funcs_maximums_loc[f,])
    maximum <- funcs_maximums[f]
    total_domain_area <- (ub[1] - lb[1]) * (ub[2] - lb[2]) 
    for (w in 1:nrow(waves_configs)){
      wave_config <- c(waves_configs[w,])
      save_path <- paste("Plots/Chapter 5 - Bayesian optimization/5.2 - Multi-points optimization/", 
                         func_names[f], "/", paste(wave_config[1], "_", wave_config[2], "_",
                                                   wave_config[3], "_", wave_config[4], sep=""),
                         "/", sep="")
      if (!dir.exists(save_path)) dir.create(save_path)
      optimization <- run_implausibility_optimization(func, init_points, D, 4900,
                                                      hyperparameters, lb, ub, wave_config,  
                                                      1e-4, maximum_loc, dim=2,
                                                      save_path=save_path)
      new_points <- get_points(10000, lb, ub, hypercube=FALSE)
      final_I <- run_history_match(optimization$optimal_point$Value, 1e-4, 0, new_points,
                                   updated_moments(optimization$init_points, new_points,
                                                   optimization$D, hyperparameters[1],
                                                   hyperparameters[2], hyperparameters[3], dim=2),
                                   plot=FALSE, optim=TRUE)
      mean_vars[w, f] <- mean(optimization$new_moments$Var)
      final_nir_areas[w, f] <- 100 * length(final_I[final_I <= 3]) / 10000
      final_maximum_distances[w, f] <- abs(maximum - optimization$optimal_point$Value)
      final_loc_distances[w, f] <- sqrt((optimization$optimal_point$Point[1] - maximum_loc[1])^2 + (optimization$optimal_point$Point[2] - maximum_loc[2])^2) 
    }
  }
  final_output <- matrix(0, nrow=nrow(waves_configs), ncol=4)
  for (w in 1:nrow(waves_configs)){
    final_output[w, 1] <- mean(mean_vars[w,])
    final_output[w, 2] <- mean(final_nir_areas[w,])
    final_output[w, 3] <- mean(final_maximum_distances[w,])
    final_output[w, 4] <- mean(final_loc_distances[w,])
  }
  colnames(final_output) <- c("Var", "NIR", "Max distance", "Max loc distance")
  return (list(var=mean_vars, nir=final_nir_areas, max_d=final_maximum_distances,
               loc=final_loc_distances, final=final_output))
}

trid <- function(points, scaling=NULL){
  f <- if (is.null(nrow(points))) sum((points - 1)^2) - prod(points) else (points[,1] - 1)^2 + (points[,2] - 1)^2 - points[,1] * points[,2]
  return (-1 * f)
}

sphere <- function(points, scaling=NULL){
  f <- if (is.null(nrow(points))) sum(points^2) else (points[,1]^2 + points[,2]^2)
  return (-1 * f)
}

sum_squares <- function(points, scaling=NULL){
  f <- if (is.null(nrow(points))) (points[1]^2 + 2*points[2]^2) else (points[,1]^2 + 2*points[,2]^2)
  return (-1 * f)
}

if (sys.nframe() == 0){
  ## 5.2 - Multi-points optimization
  ## 5.2.2 - Optimizing a simple two-dimensional function
  n <- 10
  sigma_u <- 0.3
  theta <- 1.4
  beta_0 <- 0
  lb <- rep(0, 2)
  ub <- rep(2*pi, 2)
  init_points <- get_points(n, lb, ub, hypercube=TRUE, maximin=TRUE)
  D <- two_D_func(init_points)
  n_new_points <- 2500
  optima <- c(pi/2, pi)
  waves_config <- c(5, 5, 5, 5, 5)
  save_path <- "Plots/Chapter 5 - Bayesian optimization/5.2 - Multi-points optimization/Initial/"
  optimization <- run_implausibility_optimization(two_D_func, init_points, D, n_new_points,
                                                  c(beta_0, sigma_u, theta), lb, ub, waves_config,  
                                                  1e-4, optima, zlim=c(-2, 2), dim=2,
                                                  save_path=save_path)
  
  # 5.2.3 - Investigating different wave configurations
  waves_configs <- matrix(c(8, 8, 8, 8,     12, 8, 7, 5,      5, 7, 8, 12), nrow=3, byrow=TRUE)
  func_names <- c("Trid", "Sphere", "Sum Squares")
  funcs_hyperparams <- matrix(c(-12, 10, 1.6,    -17, 10, 1.85,    -25, 10, 2), nrow=3, byrow=TRUE)
  funcs_bounds <- matrix(c(-4, -4, 4, 4,    -5.12, -5.12, 5.12, 5.12,     -5.12, -5.12, 5.12, 5.12),
                         nrow=3, byrow=TRUE)
  funcs_maximums <- c(2, 0, 0)
  funcs_maximums_loc <- matrix(c(2, 2,    0, 0,     0, 0), nrow=3, byrow=TRUE)
  
  analysis <- run_wave_config_analysis(waves_configs, funcs_hyperparams, funcs_bounds,
                                       funcs_maximums, funcs_maximums_loc, func_names)
}