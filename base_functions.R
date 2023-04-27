one_D_func <- function(points){
  return (points*sin(5*points))
} 

two_D_func <- function(points, scaling=NULL){
  f <- if (is.null(nrow(points))) sin(points[1]) - cos(points[2]) else sin(points[,1]) - cos(points[,2])
  return (f)
}

updated_moments <- function(init_points, new_points, D, beta_0, sigma_u, theta, calc_var_only=FALSE, 
                            dim=1){
  n <- if (dim == 1) length(init_points) else nrow(init_points)
  m <- if (dim == 1) length(new_points) else nrow(new_points)

  E_fx <- rep(beta_0, times=m)
  
  Var_fx <- diag(sigma_u^2, nrow=m, ncol=m)
  
  E_D <- rep(beta_0, times=n)
  
  new_points_dist <- if (dim <= 2) cdist(new_points, init_points) else cdist(new_points[,1:dim], 
                                                                             init_points[,1:dim])
  Cov_fx_D <- matrix(0, nrow=m, ncol=n)
  for (i in 1:m){
    Cov_fx_D[i,] <- sigma_u^2 * exp(-1 * ((new_points_dist[i,])^2) / theta^2)
  }
  
  init_points_dist <- if (dim <=2) pdist(init_points) else pdist(init_points[,1:dim])
  delta <- 10e-6
  Var_D <- diag(delta, nrow=n, ncol=n)
  for (i in 1:n){
    Var_D[i,] <- Var_D[i,] + sigma_u^2 * exp(-1 * ((init_points_dist[i,])^2) / theta^2)
  }
  Var_D_inv <- solve(Var_D)
  
  if (dim <= 2){
    E_fx_new <- if (!calc_var_only) E_fx + Cov_fx_D%*%Var_D_inv%*%(D - E_D) else 0
    Var_fx_new <- Var_fx - Cov_fx_D%*%Var_D_inv%*%t(Cov_fx_D)
    if (dim == 2){
      E_fx_new <- if (!calc_var_only) matrix(E_fx_new, nrow=sqrt(m)) else 0
      Var_fx_new <- matrix(Var_fx_new, nrow=m)
    }
    Var_fx_new <- if (dim == 1) diag(Var_fx_new) else matrix(diag(Var_fx_new), nrow=sqrt(m))
    return (list(Exp=E_fx_new, Var=Var_fx_new))
  }
  else{
    Exp <- rep(0, m) 
    Var <- rep(0, m)
    for (i in 1:m){
      if(!calc_var_only) Exp[i] <- beta_0 + Cov_fx_D[i,]%*%Var_D_inv%*%(D - E_D)
      Var[i] <- sigma_u^2 - Cov_fx_D[i,]%*%Var_D_inv%*%t(Cov_fx_D)[,i]
    }
    return (list(Exp=Exp, Var=Var))
  }
}

run_emulation <- function(init_points, new_points, D, beta_0, sigma_u, theta, plots=TRUE, 
                          add_run_points=FALSE, run_cols=NULL, variance=TRUE, true_func=FALSE, 
                          func=NULL, xlim=NULL, ylim=NULL, zlim=NULL, save_path="", i="", dim=1){
  emulation <- updated_moments(init_points, new_points, D, beta_0, sigma_u, theta, dim=dim)
  E_fx_new <- if (dim <= 2) emulation$Exp
  Var_fx_new <- if (dim <= 2) emulation$Var
  
  if (plots == TRUE & dim <= 2){
    if (save_path != "") if (!dir.exists(save_path)) dir.create(save_path)
    if (dim == 2){
      x <- sort(unique(new_points[,1]), decreasing=FALSE)
      y <- sort(unique(new_points[,2]), decreasing=FALSE)
      if (true_func & !is.null(func)){
        img_path <- paste(save_path, "Function.png", sep="")
        if (save_path != "") png(img_path, width=504, height=504)
        z <- matrix(func(new_points), nrow=length(x))
        filled.contour(x=x, y=y, z=z, xlab="x", ylab="y", cex.lab=1.3,
                       cex.axis=1.3, color.palette=terrain.colors)
        if (save_path != "") dev.off()
        if (save_path != "") image_write(image_crop(image_read(img_path), "504x460+0+44"), 
                                         path=img_path, format="png")
        zlim <- if (is.null(zlim)) c(min(z), max(z)) else zlim
      }
      img_path <- paste(save_path, "Expectation", i, ".png", sep="")
      if (save_path != "") png(img_path, width=504, height=504)
      zlim <- if (is.null(zlim)) range(E_fx_new, finite=TRUE) else zlim
      if (add_run_points){
        filled.contour(x=x, y=y, z=E_fx_new, xlab="x", ylab="y", 
                       plot.axes={axis(1);axis(2);points(init_points[,1], init_points[,2], pch=19, col=if (is.null(run_cols)) "blue" else run_cols)},
                       cex.axis=1.3, cex.lab=1.3, color.palette=terrain.colors, zlim=zlim)
      }
      else{
        filled.contour(x=x, y=y, z=E_fx_new, xlab="x", ylab="y", zlim=zlim, cex.axis=1.3, 
                       cex.lab=1.3, color.palette=terrain.colors)
      }
      if (save_path != "") dev.off()
      if (save_path != "") image_write(image_crop(image_read(img_path), "504x460+0+44"), 
                                       path=img_path, format="png")

      # Variance contour plot
      img_path <- paste(save_path, "Variance", i, ".png", sep="")
      if (save_path != "") png(img_path, width=504, height=504)
      filled.contour(x=x, y=y, z=3*sqrt(Var_fx_new), xlab="x", ylab="y", 
                     cex.axis=1.3, cex.lab=1.3, color.palette=cm.colors)
      if (save_path != "") dev.off()
      if (save_path != "") image_write(image_crop(image_read(img_path), "504x460+0+44"), 
                                       path=img_path, format="png")
    }
    else if (dim == 1){
      if (save_path != "") pdf(paste(save_path, "Emulation", i, ".pdf", sep=""))
      xlim <- if (is.null(xlim)) c(min(points), max(points)) else xlim
      ylim <- if (is.null(ylim)) c(min(E_fx_new), max(E_fx_new)) else ylim
      plot(new_points, E_fx_new, col="blue", xlim=xlim, ylim=ylim, xlab="x", ylab="f(x)", cex.lab=1.3, 
           cex.axis=1.3, type="l")
      if (add_run_points){
        points(init_points, D, pch=19, col=if (is.null(run_cols)) "blue" else run_cols)
      }
      if (variance){
        lines(new_points, E_fx_new + 3 * sqrt(Var_fx_new), col="red")
        lines(new_points, E_fx_new - 3 * sqrt(Var_fx_new), col="red")
      }
      if (true_func == TRUE & !is.null(func)){
        lines(new_points, func(new_points), col="green")
      }
      
      if (save_path != "") dev.off()
    }
  }
  return (emulation)
}

implausibility <- function(z, E_fx_new, Var_fx_new, epsilon, err, optim=FALSE){
  I_sq <- (E_fx_new - z)^2 / (Var_fx_new + epsilon + err)
  I <- sqrt(I_sq)
  if (optim){
    idx <- which(E_fx_new > z, arr.ind=TRUE)
    I[idx] <- 0
  }
  return (I)
}

get_points <- function(n, lb, ub, dim=2, hypercube=TRUE, maximin=TRUE){
  if (hypercube){
    D <- if (maximin) maximinSLHD(t = 1, m = n, k = dim, power=5, nstarts=1)$StandDesign else randomLHS(n, dim)
    for (i in 1:dim){
      D[,i] <- lb[i] + D[,i] * (ub[i] - lb[i])
    }
    return (D)
  }
  else{
    x <- seq(lb[1], ub[1], length.out=sqrt(n))
    y <- seq(lb[2], ub[2], length.out=sqrt(n))
    points <- as.matrix(expand.grid(x, y))
    return (points)
  }
}

PI <- function(E_fx_new, Var_fx_new, fx_top, gamma=0){
  i <- (E_fx_new - fx_top - gamma) / sqrt(Var_fx_new)
  return (pnorm(i))
}

EI <- function(E_fx_new, Var_fx_new, fx_top, gamma=0){
  Z <- (E_fx_new - fx_top - gamma) / sqrt(Var_fx_new)
  E <- sqrt(Var_fx_new) * dnorm(Z) + (E_fx_new - fx_top - gamma) * pnorm(Z)
  return (E)
}

EI2 <- function(E_fx_new, Var_fx_new, fx_top, ...){
  Z <- (E_fx_new - fx_top) / sqrt(Var_fx_new)
  E2 <- Var_fx_new * ((Z^2 + 1) * pnorm(Z) + Z * dnorm(Z))
  return (E2)
}

EI3 <- function(E_fx_new, Var_fx_new, fx_top, ...){
  Z <- (E_fx_new - fx_top) / sqrt(Var_fx_new)
  E3 <- Var_fx_new^(3/2) * ((Z^3 + 3 * Z) * pnorm(Z) + (Z^2 + 2) * dnorm(Z))
  return (E3)
}