library(rdist)
library(SLHD)
library(magick)
library(lhs)
library(MASS)
import::from(base_functions.R, two_D_func, updated_moments, run_emulation, get_points, implausibility)

design <- function(n, lb, ub, file_name, hypercube=TRUE, maximin=TRUE){
  points <- get_points(n, lb, ub, hypercube=hypercube, maximin=maximin)
  x <- points[,1]
  y <- points[,2]
  pdf(file_name)
  plot(x, y, pch=19, col="blue", xlim=c(0, 2*pi), ylim=c(0, 2*pi), xlab="x", ylab="y", cex.lab=1.3, 
       cex.axis=1.3)
  if (hypercube){
    abline(h=seq(0, 2*pi, length.out=n+1), v=seq(0, 2*pi, length.out=n+1), col="black")
  }
  dev.off()
  return (points)
}

if (sys.nframe() == 0){
  ## 4.1 - Basic design in two dimensions
  save_path <- "Plots/Chapter 4 - Extension to higher dimensions/4.1 - Basic design in two dimensions/"
  lb <- rep(0, 2)
  ub <- rep(2*pi, 2)
  init_points_grid <- design(16, lb, ub, paste(save_path, "Grid design.pdf", sep=""), hypercube=FALSE)
  new_points <- design(1600, lb, ub, paste(save_path, "Grid design emulator.pdf", sep=""),
                       hypercube=FALSE)
  n <- 16
  sigma_u <- 0.3
  theta <- 1.4
  beta_0 <- 0
  D <- two_D_func(init_points_grid)
  new_moments_grid <- run_emulation(init_points_grid, new_points, D, beta_0, sigma_u, theta, 
                                    save_path=save_path, dim=2, zlim=c(-2, 2), true_func=TRUE,
                                    func=two_D_func)
  
  ## 4.2 - Improving points selection
  save_path <- "Plots/Chapter 4 - Extension to higher dimensions/4.2 - Improving points selection/"
  init_points_rlhd <- design(16, lb, ub, paste(save_path, "Random LHD.pdf", sep=""), hypercube=TRUE,
                             maximin=FALSE)
  new_moments_rlhd <- run_emulation(init_points_rlhd, new_points, two_D_func(init_points_rlhd), beta_0, 
                                    sigma_u, theta, save_path=save_path, dim=2, zlim=c(-2, 2),
                                    func=two_D_func, i=2)
  init_points_mlhd <- design(16, lb, ub, paste(save_path, "Maximin LHD.pdf", sep=""), hypercube=TRUE,
                             maximin=TRUE)
  new_moments_mlhd <- run_emulation(init_points_mlhd, new_points, two_D_func(init_points_mlhd), beta_0, 
                                    sigma_u, theta, save_path=save_path, dim=2, zlim=c(-2, 2),
                                    func=two_D_func, i=3)
  mean_var_grid <- mean(new_moments_grid$Var)
  mean_var_rlhd <- mean(new_moments_rlhd$Var)
  mean_var_mlhd <- mean(new_moments_mlhd$Var)
}


## 4.3 - Optimizing wave design
run_history_match <- function(z, epsilon, err, new_points, new_moments, plot=TRUE, dim=2, 
                              save_path="", i="", optima=NULL, optim=FALSE){
  x <- sort(unique(new_points[,1]), decreasing=FALSE)
  y <- sort(unique(new_points[,2]), decreasing=FALSE)
  E_fx_new <- new_moments$Exp
  Var_fx_new <- new_moments$Var
  I <- implausibility(z, E_fx_new, Var_fx_new, epsilon, err, optim=optim)
  if (plot & dim == 2){
    img_path <- paste(save_path, "Implausibility", i, ".png", sep="")
    if (save_path != "") if (!dir.exists(save_path)) dir.create(save_path)
    if (save_path != "") png(img_path, width=504, height=504)
    if (is.null(optima)) filled.contour(x=x, y=y, z=I, zlim=c(0, 3), xlab="x", ylab="y", cex.lab=1.3, cex.axis=1.3) else filled.contour(x=x, y=y, z=I, zlim=c(0, 3), plot.axes={axis(1);axis(2);points(optima[1], optima[2], pch=19, col="blue")}, xlab="x", ylab="y",  cex.axis=1.3, cex.lab=1.3)
    if (save_path != "") dev.off()
    if (save_path != "") image_write(image_crop(image_read(img_path), "504x460+0+44"), path=img_path,
                                     format="png")
  }
  return (I)
}

plot_wave_points <- function(wave_points, n_wave_points, new_points, I, i="", save_path=""){
  x <- sort(unique(new_points[,1]), decreasing=FALSE)
  y <- sort(unique(new_points[,2]), decreasing=FALSE)
  img_path <- paste(save_path, "Wave points", i, ".png", sep="")
  if (save_path != "") if (!dir.exists(save_path)) dir.create(save_path)
  if (save_path != "") png(img_path, width=504, height=504)
  x_wave_points <- if (n_wave_points == 1) wave_points[1] else wave_points[,1]
  y_wave_points <- if (n_wave_points == 1) wave_points[2] else wave_points[,2]
  filled.contour(x=x, y=y, z=I, zlim=c(0, 3), 
                 plot.axes={axis(1);axis(2);points(x_wave_points, y_wave_points, pch=19, col="blue")},
                 xlab="x", ylab="y", cex.axis=1.3, cex.lab=1.3)
  if (save_path != "") dev.off()
  if (save_path != "") image_write(image_crop(image_read(img_path), "504x460+0+44"), path=img_path,
                                   format="png")
}

get_non_implausible_ind <- function(I){
  non_implausible_thresholds <- c(5, 7, 10, 12, 15, 17, 20)
  non_implausible <- which(I <= 3, arr.ind=TRUE)
  p <- 1
  while (length(non_implausible) == 0){
    non_implausible <- which(I <= non_implausible_thresholds[p], arr.ind=TRUE)
    p <- p + 1
  }
  return (non_implausible)
}

sdd <- function(n_wave_points, I, new_points, save_path="", i=""){
  diff <- function(v){
    p <- t(combn(v, 2))
    return (sum(abs(p[,1] - p[,2])))
  }
  
  x <- sort(unique(new_points[,1]), decreasing=FALSE)
  y <- sort(unique(new_points[,2]), decreasing=FALSE)
  non_implausible <- get_non_implausible_ind(I)
  non_implausible_points <- new_points[length(x) * (non_implausible[,2] - 1) + non_implausible[,1],]
  
  if (nrow(non_implausible_points) <= n_wave_points){
    return (non_implausible_points)
  }
  else{
    wave_points <- matrix(0, nrow=n_wave_points, ncol=2)
    pairwise_dists <- pdist(non_implausible_points)
    selected_points_idx <- rep(0, n_wave_points)
    for (j in 1:n_wave_points){
      idx <- 0 # index of new point
      if (j == 1){
        idx <- sample.int(n_wave_points, size=1)
      }
      else if (j == 2){
        idx <- which.max(pairwise_dists[selected_points_idx[1],])
      }
      else{
        distances <- t(pairwise_dists[selected_points_idx,])
        sum_diffs <- apply(distances, 1, diff)
        sum_diffs[selected_points_idx] <- Inf
        idx <- which.min(sum_diffs)
      }
      selected_points_idx[j] <- idx
      wave_points[j,] <- non_implausible_points[idx,]
    }
    
    plot_wave_points(wave_points, n_wave_points, new_points, I, i=i, save_path=save_path)
    
    return (wave_points)
  }
  
}

implausibility_sampling <- function(n_wave_points, I, new_points, Var_fx_new, save_path="", i=""){
  x <- sort(unique(new_points[,1]), decreasing=FALSE)
  y <- sort(unique(new_points[,2]), decreasing=FALSE)
  non_implausible <- get_non_implausible_ind(I)
  non_implausible_points <- new_points[length(x) * (non_implausible[,2] - 1) + non_implausible[,1],]
  weights <- I[non_implausible]
  wave_points_ind <- sample(1:nrow(non_implausible_points), n_wave_points, prob=weights)
  wave_points <- non_implausible_points[wave_points_ind,]
  plot_wave_points(wave_points, n_wave_points, new_points, I, i=i, save_path=save_path)
  return (wave_points)
}

variance_sampling <- function(n_wave_points, I, init_points, new_points, D, Var_fx_new, dim=2,  
                             save_path="", i=""){
  n <- nrow(init_points)
  x <- if (dim == 2) sort(unique(new_points[,1]), decreasing=FALSE) else NULL
  y <- if (dim == 2) sort(unique(new_points[,2]), decreasing=FALSE) else NULL
  non_implausible <- get_non_implausible_ind(I)
  non_implausible_points <- if (dim == 2) new_points[length(x) * (non_implausible[,2] - 1) + non_implausible[,1],] else new_points[non_implausible,1:dim]
  n_wave_points <- if (length(non_implausible) < n_wave_points) length(non_implausible) else n_wave_points
  for (j in 1:n_wave_points){
    variances <- Var_fx_new[non_implausible]
    next_point_ind <- which.max(variances)
    next_point <- non_implausible_points[next_point_ind,]
    if (dim > 2) names(next_point) <- names(init_points)
    init_points <- rbind(init_points, next_point)
    D <- append(D, 0)
    if (j != n_wave_points){
      Var_fx_new <- updated_moments(init_points, new_points, D, beta_0, sigma_u, theta,
                                    calc_var_only=TRUE, dim=dim)$Var
    }
  }
  if (dim == 2) plot_wave_points(init_points[n + 1:n_wave_points, 1:dim], n_wave_points, new_points, I, i=i, save_path=save_path)
  return (init_points[n + 1:n_wave_points, 1:dim])
}

run_next_wave <- function(n_wave_points, func, I, new_moments, init_points, new_points, D, beta_0, 
                          sigma_u, theta, z, epsilon, err, algorithm="variance", zlim=NULL,
                          save_path="", i="", optima=NULL, optim=FALSE, dim=2, 
                          scaling=list(lb=rep(0, dim), ub=rep(1, dim))){
  Var_fx_new <- new_moments$Var
  next_wave_points <- NULL
  if (algorithm == "variance"){
    next_wave_points <- variance_sampling(n_wave_points, I, init_points, new_points, D, 
                                         Var_fx_new, save_path=save_path, i=i, dim=dim)
  }
  else if (algorithm == "implausibility"){
    next_wave_points <- implausibility_sampling(n_wave_points, I, new_points, 
                                                Var_fx_new, save_path=save_path, i=i)
  }
  else if (algorithm == "sdd"){
    next_wave_points <- sdd(n_wave_points, I, new_points, save_path=save_path, i=i)
  }
  init_points <- rbind(init_points, next_wave_points)
  D <- append(D, func(next_wave_points, scaling=scaling))
  new_moments <- run_emulation(init_points, new_points, D, beta_0, sigma_u, theta, save_path=save_path,
                               i=i, add_run_points=TRUE, dim=dim, zlim=zlim)
  I <- run_history_match(z, epsilon, err, new_points, new_moments, save_path=save_path, i=i, 
                         optima=optima, optim=optim, dim=dim)
  wave <- list(init_points=init_points, D=D, new_moments=new_moments, I=I, NIR_size=length(I[I<=3]), NIR_prop=length(I[I<=3])/length(I))
  return (wave)
}

if (sys.nframe() == 0){
  save_path <- "Plots/Chapter 4 - Extension to higher dimensions/4.3 - Optimizing wave design/"
  ## Figure 4.7
  z <- 0.1
  epsilon <- 2e-3
  err <- 0
  I <- run_history_match(z, epsilon, err, new_points, new_moments_mlhd, save_path=save_path)
  length(I[I<=3])
  
  ### 4.3.1 - SDD
  save_path <- "Plots/Chapter 4 - Extension to higher dimensions/4.3 - Optimizing wave design/SDD/"
  wave2_s <- run_next_wave(5, two_D_func, I, new_moments_mlhd, init_points_mlhd, new_points, D, beta_0,
                           sigma_u, theta, z, epsilon, err, algorithm="sdd", save_path=save_path, i=2,
                           zlim=c(-2, 2))
  wave3_s <- run_next_wave(5, two_D_func, wave2_s$I, wave2_s$new_moments, wave2_s$init_points,
                           new_points, wave2_s$D, beta_0, sigma_u, theta, z, epsilon, err,
                           algorithm="sdd", save_path=save_path, i=3, zlim=c(-2, 2))
  wave4_s <- run_next_wave(5, two_D_func, wave3_s$I, wave3_s$new_moments, wave3_s$init_points,
                           new_points, wave3_s$D, beta_0, sigma_u, theta, z, epsilon, err, 
                           algorithm="sdd", save_path=save_path, i=4, zlim=c(-2, 2))
  
  ### 4.3.2 - Implausibility sampling
  save_path <- "Plots/Chapter 4 - Extension to higher dimensions/4.3 - Optimizing wave design/Implausibility sampling/"
  wave2_i <- run_next_wave(5, two_D_func, I, new_moments_mlhd, init_points_mlhd, new_points, D, beta_0,
                           sigma_u, theta, z, epsilon, err, algorithm="implausibility",
                           save_path=save_path, i=2, zlim=c(-2, 2))
  wave3_i <- run_next_wave(5, two_D_func, wave2_i$I, wave2_i$new_moments, wave2_i$init_points,
                           new_points, wave2_i$D, beta_0, sigma_u, theta, z, epsilon, err,
                           algorithm="implausibility", save_path=save_path, i=3, zlim=c(-2, 2))
  wave4_i <- run_next_wave(5, two_D_func, wave3_i$I, wave3_i$new_moments, wave3_i$init_points,
                           new_points, wave3_i$D, beta_0, sigma_u, theta, z, epsilon, err, 
                           algorithm="implausibility", save_path=save_path, i=4, zlim=c(-2, 2))
  
  ### 4.3.3 - Variance sampling
  save_path <- "Plots/Chapter 4 - Extension to higher dimensions/4.3 - Optimizing wave design/Variance sampling/"
  wave2_v <- run_next_wave(5, two_D_func, I, new_moments_mlhd, init_points_mlhd, new_points, D, beta_0,
                           sigma_u, theta, z, epsilon, err, algorithm="variance", save_path=save_path,
                           i=2, zlim=c(-2, 2))
  wave3_v <- run_next_wave(5, two_D_func, wave2_v$I, wave2_v$new_moments, wave2_v$init_points,
                           new_points, wave2_v$D, beta_0, sigma_u, theta, z, epsilon, err,
                           algorithm="variance", save_path=save_path, i=3, zlim=c(-2, 2))
  wave4_v <- run_next_wave(5, two_D_func, wave3_v$I, wave3_v$new_moments, wave3_v$init_points,
                           new_points, wave3_v$D, beta_0, sigma_u, theta, z, epsilon, err, 
                           algorithm="variance", save_path=save_path, i=4, zlim=c(-2, 2))
}