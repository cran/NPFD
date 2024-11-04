#' N-Power Fourier Deconvolution
#'
#'
#' @description Estimates the density \eqn{f_y}, given vectors \eqn{x} and \eqn{z}, where \eqn{f_z} results from the convolution of \eqn{f_x} and \eqn{f_y}.
#'
#' @param x Vector of observations for \eqn{x}.
#' @param z Vector of observations for \eqn{z}.
#' @param mode Deconvolution mode (\code{empirical} or \code{denspr}). If \code{empirical}, the Fourier transforms of \eqn{x} and \eqn{z} are estimated using the empirical form. If \code{denspr}, they are calculated based on the density estimations using \code{densprf} (see the package \strong{siggenes}).
#' @param dfx Degrees of freedom for the estimation of \eqn{f_x} if mode is set to \code{denspr}.
#' @param dfz Degrees of freedom for the estimation of \eqn{f_z} if mode is set to \code{denspr}.
#' @param Lx Number of points for \eqn{f_x}-grid if mode is set to \code{denspr}.
#' @param Lz Number of points for \eqn{f_z}-grid if mode is set to \code{denspr}.
#' @param Ly Number of points for \eqn{f_y}-grid.
#' @param N Possible power values.
#' @param FT.grid Vector of grid for Fourier transformation of \eqn{f_x} and \eqn{f_z}.
#' @param lambda Smoothing parameter.
#' @param eps Tolerance for convergence.
#' @param delta Small margin value.
#' @param error Error model (\code{unknown}, \code{normal}, \code{laplacian}). If \code{unknown}, the Fourier transform of \eqn{x} is calculated based on the mode. If \code{normal}, the exact form of the Fourier transform of a centered normal distribution with standard deviation sigma is used for \eqn{x}. If \code{laplacian}, the exact form of the Fourier transform of a centered Laplace distribution with standard deviation sigma is used for \eqn{x}.
#' @param sigma Standard deviation for normal or Laplacian error.
#' @param calc.error Logical indicating whether to calculate error (10 x ISE between \eqn{f_z} and \eqn{f_x * f_y}).
#' @param plot Logical indicating whether to plot \eqn{f_z} vs. \eqn{f_x * f_y} if \code{calc.error} is \code{TRUE}.
#' @param legend Logical indicating whether to include a legend in the plot if \code{calc.error} is \code{TRUE}.
#' @param positive Logical indicating whether to enforce non-negative density estimation.
#'
#' @return A list with the following components:
#' \item{\code{x}}{A vector of \eqn{x}-values of the resulting density estimation.}
#' \item{\code{y}}{A vector of \eqn{y}-values of the resulting density estimation.}
#' \item{\code{N}}{The power used in the deconvolution process.}
#' \item{\code{error}}{The calculated error if \code{calc.error = TRUE}.}
#'
#' @author Akin Anarat \email{akin.anarat@hhu.de}
#'
#' @references Anarat A., Krutmann, J., and Schwender, H. (2024). A nonparametric statistical method for deconvolving densities in the analysis of proteomic data. Submitted.
#'
#' @importFrom stats dnorm integrate sd
#' @importFrom VGAM dlaplace
#' @importFrom graphics lines par
#'
#' @export
#'
#' @examples
#' # Deconvolution when mixed data and data from an independent experiment are provided:
#' set.seed(123)
#' x <- rnorm(1000)
#' y <- rgamma(1000, 10, 2)
#' z <- x + y
#'
#' f <- function(x) dgamma(x, 10, 2)
#'
#' independent.x <- rnorm(100)
#'
#' fy.NPFD <- deconvolve(independent.x, z, calc.error = TRUE, plot = TRUE)
#' x.x <- fy.NPFD$x
#' fy <- f(x.x)
#'
#' # Check power and error values
#' fy.NPFD$N
#' fy.NPFD$error
#'
#' # Plot density functions
#' plot(NULL, xlim = range(y), ylim = c(0, max(fy, fy.NPFD$y)), xlab = "x", ylab = "Density")
#' lines(x.x, fy, col = "blue", lwd = 2)
#' lines(fy.NPFD, col = "orange", lwd = 2)
#' legend("topright", legend = c(expression(f[y]), expression(f[y]^{NPFD})),
#'        col = c("blue", "orange"), lwd = c(2, 2))
#'
#' # For replicated mixed data:
#' set.seed(123)
#' x1 <- VGAM::rlaplace(1000, 0, 1/sqrt(2))
#' x2 <- VGAM::rlaplace(1000, 0, 1/sqrt(2))
#' y <- rgamma(1000, 10, 2)
#' z1 <- z <- x1 + y
#' z2 <- x2 + y
#'
#' x <- createSample(z1, z2)
#'
#' fy.NPFD <- deconvolve(x, z, mode = "denspr", calc.error = TRUE, plot = TRUE)
#' x.x <- fy.NPFD$x
#' fy <- f(x.x)
#'
#' # Check power and error values
#' fy.NPFD$N
#' fy.NPFD$error
#'
#' # Plot density functions
#' plot(NULL, xlim = range(y), ylim = c(0, max(fy, fy.NPFD$y)), xlab = "x", ylab = "Density")
#' lines(x.x, fy, col = "blue", lwd = 2)
#' lines(fy.NPFD, col = "orange", lwd = 2)
#' legend("topright", legend = c(expression(f[y]), expression(f[y]^{NPFD})),
#'        col = c("blue", "orange"), lwd = c(2, 2))
#'
#' # When the distribution of x is asymmetric and the sample size is very small:
#' set.seed(123)
#' x <- rgamma(5, 4, 2)
#' y <- rgamma(1000, 10, 2)
#' z <- x + y
#'
#' fy.NPFD <- deconvolve(x, z, mode = "empirical", lambda = 2)
#' x.x <- fy.NPFD$x
#' fy <- f(x.x)
#'
#' # Check power value
#' fy.NPFD$N
#'
#' # Plot density functions
#' plot(NULL, xlim = range(y), ylim = c(0, max(fy, fy.NPFD$y)), xlab = "x", ylab = "Density")
#' lines(x.x, fy, col = "blue", lwd = 2)
#' lines(fy.NPFD, col = "orange", lwd = 2)
#' legend("topright", legend = c(expression(f[y]), expression(f[y]^{NPFD})),
#'        col = c("blue", "orange"), lwd = c(2, 2))

deconvolve <- function(x = NULL, z, mode = c("empirical", "denspr"), dfx = 5, dfz = 5,
                       Lx = 10^2, Lz = 10^2, Ly = 10^2, N = 1:100, FT.grid = seq(0, 100, 0.1),
                       lambda = 1, eps = 10^-3, delta = 10^-2, error = c("unknown", "normal", "laplacian"),
                       sigma = NULL, calc.error = FALSE, plot = FALSE, legend = TRUE, positive = FALSE) {

  mode <- match.arg(mode)
  error <- match.arg(error)

  if (!is.null(x) && error != "unknown") {
    warning("No assumptions about error distribution can be made if x is provided; error is set to 'unknown'.")
    error <- "unknown"
  }
  if (!is.null(x) && !is.null(sigma)) {
    warning("No assumptions about sigma can be made if x is provided; argument sigma is not used.")
    sigma <- NULL
  }

  K <- 0

  # Deconvolution process based on the selected mode
  if (mode == "denspr") {
    for (i in seq(N)) {  # Iterate over different values of k in N
      k <- N[i]

      # Scaling and translation of x and z based on the current k
      a <- 1 / sqrt(k)
      if (error == "unknown") {
        bx <- (1/k - 1/sqrt(k)) * mean(x)
        xr <- a * x + bx  # Transformed x
      } else {
        xr <- a * c(-sigma, sigma) / sqrt(2)
      }

      bz <- (1/k - 1/sqrt(k)) * mean(z)
      zr <- a * z + bz  # Transformed z

      # Define the x and z intervals based on the error model
      if (error == "unknown") {
        x1 <- min(xr)
        x2 <- max(xr)
        x.x <- seq(x1, x2, length.out = Lx)
      } else if (error == "normal") {
        if (is.null(sigma)) stop("Argument sigma is missing.")
        x1 <- -6 * a * sigma
        x2 <- -x1
        x.x <- seq(x1, x2, length.out = Lx)
      } else if (error == "laplacian") {
        if (is.null(sigma)) stop("Argument sigma is missing.")
        x1 <- -8 * a * sigma
        x2 <- -x1
        x.x <- seq(x1, x2, length.out = Lx)
      }

      z1 <- min(zr)
      z2 <- max(zr)
      z.x <- seq(z1, z2, length.out = Lz)

      # Define the density functions for x and z
      if (error == "unknown") {
        fx <- densprf(xr, df = dfx)
      } else if (error == "normal") {
        fx <- function(x) dnorm(x, 0, sd(xr))
      } else if (error == "laplacian") {
        fx <- function(x) dlaplace(x, 0, sd(xr) / sqrt(2))
      }

      fz <- densprf(zr, df = dfz)

      # Normalize the density functions over the defined intervals
      fx <- fx(x.x) / integrate(fx, x1, x2)$value
      fz <- fz(z.x) / integrate(fz, z1, z2)$value

      x0 <- (x2 - x1) / length(x.x) * sum(fx)
      z0 <- (z2 - z1) / length(z.x) * sum(fz)

      # Fourier transform based approach for deconvolution
      for (l in FT.grid) {
        ax <- (x2 - x1) / length(x.x) * sum(fx * cos(l * x.x))
        bx <- (x2 - x1) / length(x.x) * sum(fx * sin(l * x.x))
        az <- (z2 - z1) / length(z.x) * sum(fz * cos(l * z.x))
        bz <- (z2 - z1) / length(z.x) * sum(fz * sin(l * z.x))
        FT <- ((az + 1i * bz) / (ax + 1i * bx) * (x0 / z0))^k
        r <- FT.grid[which(FT.grid == l) - 1]
        if (abs(FT) > 1.001) break
        else if (abs(FT) < eps) {
          ax <- (x2 - x1) / length(x.x) * sum(fx * cos((l + delta) * x.x))
          bx <- (x2 - x1) / length(x.x) * sum(fx * sin((l + delta) * x.x))
          az <- (z2 - z1) / length(z.x) * sum(fz * cos((l + delta) * z.x))
          bz <- (z2 - z1) / length(z.x) * sum(fz * sin((l + delta) * z.x))
          FT2 <- ((az + 1i * bz) / (ax + 1i * bx) * (x0 / z0))^k
          if (abs(FT2) < eps) {
            K <- 1
            break
          }
        }
      }
      if (K == 1) break
    }

    # Rescale and recompute the density functions with the final k
    k <- k * lambda
    a <- 1 / sqrt(k)

    if (error == "unknown") {
      bx <- (1/k - 1/sqrt(k)) * mean(x)
      xr <- a * x + bx
    } else {
      xr <- a * c(-sigma, sigma) / sqrt(2)
    }

    bz <- (1/k - 1/sqrt(k)) * mean(z)
    zr <- a * z + bz

    # Recompute intervals for x and z based on the final k
    if (error == "unknown") {
      x1 <- min(xr)
      x2 <- max(xr)
      x.x <- seq(x1, x2, length.out = Lx)
    } else if (error == "normal") {
      x1 <- -6 * a * sigma
      x2 <- -x1
      x.x <- seq(x1, x2, length.out = Lx)
    } else if (error == "laplacian") {
      x1 <- -8 * a * sigma
      x2 <- -x1
      x.x <- seq(x1, x2, length.out = Lx)
    }

    z1 <- min(zr)
    z2 <- max(zr)
    z.x <- seq(z1, z2, length.out = Lz)

    # Recompute the density functions
    if (error == "unknown") {
      fx <- densprf(xr, df = dfx)
    } else if (error == "normal") {
      fx <- function(x) dnorm(x, 0, sd(xr))
    } else if (error == "laplacian") {
      fx <- function(x) dlaplace(x, 0, sd(xr) / sqrt(2))
    }

    fz <- densprf(zr, df = dfz)

    fx <- fx(x.x) / integrate(fx, x1, x2)$value
    fz <- fz(z.x) / integrate(fz, z1, z2)$value

    x0 <- (x2 - x1) / length(x.x) * sum(fx)
    z0 <- (z2 - z1) / length(z.x) * sum(fz)

    # Adjust the interval if necessary
    if (abs(FT) > 1.001)
      if (l > 1) l <- r

    if (lambda == 1) {
      q1 <- -l
      q2 <- l
      q <- seq(q1, q2, delta)
    }
    else {
      if (error == "unknown") {
        for (l in FT.grid) {
          ax <- (x2 - x1) / length(x.x) * sum(fx * cos(l * x.x))
          bx <- (x2 - x1) / length(x.x) * sum(fx * sin(l * x.x))
          az <- (z2 - z1) / length(z.x) * sum(fz * cos(l * z.x))
          bz <- (z2 - z1) / length(z.x) * sum(fz * sin(l * z.x))
          FT <- ((az + 1i * bz) / (ax + 1i * bx) * (x0 / z0))^k
          r <- FT.grid[which(FT.grid == l) - 1]
          if (abs(FT) > 1.001) break
          else if (abs(FT) < eps) {
            ax <- (x2 - x1) / length(x.x) * sum(fx * cos((l + delta) * x.x))
            bx <- (x2 - x1) / length(x.x) * sum(fx * sin((l + delta) * x.x))
            az <- (z2 - z1) / length(z.x) * sum(fz * cos((l + delta) * z.x))
            bz <- (z2 - z1) / length(z.x) * sum(fz * sin((l + delta) * z.x))
            FT2 <- ((az + 1i * bz) / (ax + 1i * bx) * (x0 / z0))^k
            if (abs(FT2) < eps) break
          }
        }
      }
      # Handle the normal error model during the interval selection
      else if (error == "normal") {
        for (l in FT.grid) {
          ax <- Re(exp(-1/2 * sigma^2 * (l * a)^2))
          bx <- Im(exp(-1/2 * sigma^2 * (l * a)^2))
          az <- Re(1/length(zr) * sum(exp(1i * l * zr)))
          bz <- Im(1/length(zr) * sum(exp(1i * l * zr)))
          FT <- ((az + 1i * bz) / (ax + 1i * bx))^k
          r <- FT.grid[which(FT.grid == l) - 1]
          if (abs(FT) > 1.001) break
          else if (abs(FT) < eps) {
            ax <- Re(exp(-1/2 * sigma^2 * (l * a)^2))
            bx <- Im(exp(-1/2 * sigma^2 * (l * a)^2))
            az <- Re(1/length(zr) * sum(exp(1i * (l + delta) * zr)))
            bz <- Im(1/length(zr) * sum(exp(1i * (l + delta) * zr)))
            FT2 <- ((az + 1i * bz) / (ax + 1i * bx))^k
            if (abs(FT2) < eps) break
          }
        }
      }

      # Handle the Laplacian error model during the interval selection
      else if (error == "laplacian") {
        for (l in FT.grid) {
          ax <- Re(1 / (1 + sigma^2 * (l * a)^2))
          bx <- Im(1 / (1 + sigma^2 * (l * a)^2))
          az <- Re(1/length(zr) * sum(exp(1i * l * zr)))
          bz <- Im(1/length(zr) * sum(exp(1i * l * zr)))
          FT <- ((az + 1i * bz) / (ax + 1i * bx))^k
          r <- FT.grid[which(FT.grid == l) - 1]
          if (abs(FT) > 1.001) break
          else if (abs(FT) < eps) {
            ax <- Re(1 / (1 + sigma^2 * (l * a)^2))
            bx <- Im(1 / (1 + sigma^2 * (l * a)^2))
            az <- Re(1/length(zr) * sum(exp(1i * (l + delta) * zr)))
            bz <- Im(1/length(zr) * sum(exp(1i * (l + delta) * zr)))
            FT2 <- ((az + 1i * bz) / (ax + 1i * bx))^k
            if (abs(FT2) < eps) break
          }
        }
      }

      # Adjust the interval if necessary
      if (abs(FT) > 1.001)
        if (l > 1) l <- r

      q1 <- -l
      q2 <- l
      q <- seq(q1, q2, delta)
    }

    # Final Fourier transform to obtain the estimated density y
    FTy <- numeric(length(q))
    for (ll in seq(along=q)) {
      l <- q[ll]
      ax <- (x2 - x1) / length(x.x) * sum(fx * cos(l * x.x))
      bx <- (x2 - x1) / length(x.x) * sum(fx * sin(l * x.x))
      az <- (z2 - z1) / length(z.x) * sum(fz * cos(l * z.x))
      bz <- (z2 - z1) / length(z.x) * sum(fz * sin(l * z.x))
      FTy[ll] <- ((az + 1i * bz) / (ax + 1i * bx) * (x0 / z0))^k
    }
    ReFTy <- Re(FTy)
    ImFTy <- Im(FTy)

    # Determine the interval for y based on z and x
    if (error == "unknown") {
      s1 <- min(z) - max(x)
      s2 <- max(z) - min(x)
    } else {
      s1 <- min(z) - max(x.x)
      s2 <- max(z) - min(x.x)
    }

    xy <- seq(s1, s2, length.out = Ly)
    y <- numeric(length(xy))
    for (ll in seq(along = xy)) {
      l <- xy[ll]
      aRe <- (q2 - q1) / length(q) * sum(ReFTy * cos(-l * q))
      bRe <- (q2 - q1) / length(q) * sum(ImFTy * sin(-l * q))
      y[ll] <- 1 / (2 * pi) * (aRe - bRe)
    }

    # Ensure positivity of the estimated density if required
    if (positive == TRUE) y[which(y < 0)] <- 0

    # Optional error calculation
    if (calc.error == FALSE) err <- NA
    else if (calc.error == TRUE) {
      if (error == "unknown") {
        x1C <- min(x)
        x2C <- max(x)
      } else if (error == "normal") {
        x1C <- -6 * sigma
        x2C <- -x1C
      } else if (error == "laplacian") {
        x1C <- -8 * sigma
        x2C <- -x1C
      }

      if (error == "unknown") fxC <- densprf(x, df = dfx)
      else if (error == "normal") fxC <- function(x) dnorm(x, 0, sd(xr))
      else if (error == "laplacian") fxC <- function(x) dlaplace(x, 0, sd(xr) / sqrt(2))

      xC1 <- seq(min(xy), min(x1C, xy), -mean(diff(xy)))[-1]
      xC2 <- seq(max(xy), max(x2C, xy), mean(diff(xy)))[-1]
      xC <- c(xC1, xy, xC2)

      yC <- c(rep(0, length(xC1)), y, rep(0, length(xC2)))

      qC <- seq(min(z), max(z), length.out = Lz)
      zE <- numeric(length(qC))
      for (l in seq(qC)) {
        zE[l] <- (max(xC) - min(xC)) / length(xC) * sum(fxC(qC[l] - xC) * yC)
      }
      densprz <- densprf(z)(qC)
      err <- 10 * (max(qC) - min(qC)) / length(qC) * sum((densprz - zE)^2)

      # Optional plotting of f_x * f_y
      if (plot == TRUE) {
        old_par <- par(no.readonly = TRUE)
        on.exit(par(old_par))

        plot(NULL, xlim = range(qC), ylim = c(0, max(densprz, zE)), xlab = "x", ylab = "Density")
        lines(qC, densprz, col = "blue", lwd = 2)
        par(new = TRUE)
        lines(qC, zE, col = "orange", lwd = 3)
        if (legend == TRUE) legend("topright", legend = c(expression(f[z]), expression(f[x] * " * " * f[y]^{NPFD})),
                                   col = c("blue", "orange"), lwd = c(2, 2))
      }
    }
  }
  else if (mode == "empirical") {
    for (i in seq(N)) {  # Iterate over different values of k in N
      k <- N[i]

      # Scaling and translating z based on the current k
      a <- 1 / sqrt(k)
      bz <- (1/k - 1/sqrt(k)) * mean(z)
      zr <- a * z + bz  # Transformed z

      # Fourier transform approach based on the error model
      if (error == "unknown") {

        # Scaling and translating x based on the current k
        bx <- (1/k - 1/sqrt(k)) * mean(x)
        xr <- a * x + bx  # Transformed x

        for (l in FT.grid) {
          ax <- Re(1/length(xr) * sum(exp(1i * l * xr)))
          bx <- Im(1/length(xr) * sum(exp(1i * l * xr)))
          az <- Re(1/length(zr) * sum(exp(1i * l * zr)))
          bz <- Im(1/length(zr) * sum(exp(1i * l * zr)))
          FT <- ((az + 1i * bz) / (ax + 1i * bx))^k
          r <- FT.grid[which(FT.grid == l) - 1]
          if (abs(FT) > 1.001) break
          else if (abs(FT) < eps) {
            ax <- Re(1/length(xr) * sum(exp(1i * (l + delta) * xr)))
            bx <- Im(1/length(xr) * sum(exp(1i * (l + delta) * xr)))
            az <- Re(1/length(zr) * sum(exp(1i * (l + delta) * zr)))
            bz <- Im(1/length(zr) * sum(exp(1i * (l + delta) * zr)))
            FT2 <- ((az + 1i * bz) / (ax + 1i * bx))^k
            if (abs(FT2) < eps) {
              K <- 1
              break
            }
          }
        }
        if (K == 1) break
      }

      # Fourier transform for normal error model
      else if (error == "normal") {

        if (is.null(sigma)) stop("Argument sigma is missing.")

        for (l in FT.grid) {
          ax <- Re(exp(-1/2 * sigma^2 * (l * a)^2))
          bx <- Im(exp(-1/2 * sigma^2 * (l * a)^2))
          az <- Re(1/length(zr) * sum(exp(1i * l * zr)))
          bz <- Im(1/length(zr) * sum(exp(1i * l * zr)))
          FT <- ((az + 1i * bz) / (ax + 1i * bx))^k
          r <- FT.grid[which(FT.grid == l) - 1]
          if (abs(FT) > 1.001) break
          else if (abs(FT) < eps) {
            ax <- Re(exp(-1/2 * sigma^2 * (l * a)^2))
            bx <- Im(exp(-1/2 * sigma^2 * (l * a)^2))
            az <- Re(1/length(zr) * sum(exp(1i * (l + delta) * zr)))
            bz <- Im(1/length(zr) * sum(exp(1i * (l + delta) * zr)))
            FT2 <- ((az + 1i * bz) / (ax + 1i * bx))^k
            if (abs(FT2) < eps) {
              K <- 1
              break
            }
          }
        }
        if (K == 1) break
      }

      # Fourier transform for Laplacian error model
      else if (error == "laplacian") {

        if (is.null(sigma)) stop("Argument sigma is missing.")

        for (l in FT.grid) {
          ax <- Re(1 / (1 + sigma^2/2 * (l * a)^2))
          bx <- Im(1 / (1 + sigma^2/2 * (l * a)^2))
          az <- Re(1/length(zr) * sum(exp(1i * l * zr)))
          bz <- Im(1/length(zr) * sum(exp(1i * l * zr)))
          FT <- ((az + 1i * bz) / (ax + 1i * bx))^k
          r <- FT.grid[which(FT.grid == l) - 1]
          if (abs(FT) > 1.001) break
          else if (abs(FT) < eps) {
            ax <- Re(1 / (1 + sigma^2/2 * (l * a)^2))
            bx <- Im(1 / (1 + sigma^2/2 * (l * a)^2))
            az <- Re(1/length(zr) * sum(exp(1i * (l + delta) * zr)))
            bz <- Im(1/length(zr) * sum(exp(1i * (l + delta) * zr)))
            FT2 <- ((az + 1i * bz) / (ax + 1i * bx))^k
            if (abs(FT2) < eps) {
              K <- 1
              break
            }
          }
        }
        if (K == 1) break
      }
    }

    # Rescale and recompute the transformed variables with the final k
    k <- k * lambda
    a <- 1 / sqrt(k)

    if(error == "unknown"){
      bx <- (1/k - 1/sqrt(k)) * mean(x)
      xr <- a * x + bx
    }

    bz <- (1/k - 1/sqrt(k)) * mean(z)
    zr <- a * z + bz

    # Adjust the interval l if necessary based on the Fourier transform value
    if (abs(FT) > 1.001)
      if (l > 1) l <- r

    # Define the interval for the Fourier transform and proceed with deconvolution
    if (lambda == 1) {
      q1 <- -l
      q2 <- l
      q <- seq(q1, q2, delta)
    }
    else {
      if (error == "unknown") {
        for (l in FT.grid) {
          ax <- Re(1/length(xr) * sum(exp(1i * l * xr)))
          bx <- Im(1/length(xr) * sum(exp(1i * l * xr)))
          az <- Re(1/length(zr) * sum(exp(1i * l * zr)))
          bz <- Im(1/length(zr) * sum(exp(1i * l * zr)))
          FT <- ((az + 1i * bz) / (ax + 1i * bx))^k
          r <- FT.grid[which(FT.grid == l) - 1]
          if (abs(FT) > 1.001) break
          else if (abs(FT) < eps) {
            ax <- Re(1/length(xr) * sum(exp(1i * (l + delta) * xr)))
            bx <- Im(1/length(xr) * sum(exp(1i * (l + delta) * xr)))
            az <- Re(1/length(zr) * sum(exp(1i * (l + delta) * zr)))
            bz <- Im(1/length(zr) * sum(exp(1i * (l + delta) * zr)))
            FT2 <- ((az + 1i * bz) / (ax + 1i * bx))^k
            if (abs(FT2) < eps) break
          }
        }
      }

      # Handle the normal error model during the interval selection
      else if (error == "normal") {
        for (l in FT.grid) {
          ax <- Re(exp(-1/2 * sigma^2 * (l * a)^2))
          bx <- Im(exp(-1/2 * sigma^2 * (l * a)^2))
          az <- Re(1/length(zr) * sum(exp(1i * l * zr)))
          bz <- Im(1/length(zr) * sum(exp(1i * l * zr)))
          FT <- ((az + 1i * bz) / (ax + 1i * bx))^k
          r <- FT.grid[which(FT.grid == l) - 1]
          if (abs(FT) > 1.001) break
          else if (abs(FT) < eps) {
            ax <- Re(exp(-1/2 * sigma^2 * (l * a)^2))
            bx <- Im(exp(-1/2 * sigma^2 * (l * a)^2))
            az <- Re(1/length(zr) * sum(exp(1i * (l + delta) * zr)))
            bz <- Im(1/length(zr) * sum(exp(1i * (l + delta) * zr)))
            FT2 <- ((az + 1i * bz) / (ax + 1i * bx))^k
            if (abs(FT2) < eps) break
          }
        }
      }

      # Handle the Laplacian error model during the interval selection
      else if (error == "laplacian") {
        for (l in FT.grid) {
          ax <- Re(1 / (1 + sigma^2/2 * (l * a)^2))
          bx <- Im(1 / (1 + sigma^2/2 * (l * a)^2))
          az <- Re(1/length(zr) * sum(exp(1i * l * zr)))
          bz <- Im(1/length(zr) * sum(exp(1i * l * zr)))
          FT <- ((az + 1i * bz) / (ax + 1i * bx))^k
          r <- FT.grid[which(FT.grid == l) - 1]
          if (abs(FT) > 1.001) break
          else if (abs(FT) < eps) {
            ax <- Re(1 / (1 + sigma^2/2 * (l * a)^2))
            bx <- Im(1 / (1 + sigma^2/2 * (l * a)^2))
            az <- Re(1/length(zr) * sum(exp(1i * (l + delta) * zr)))
            bz <- Im(1/length(zr) * sum(exp(1i * (l + delta) * zr)))
            FT2 <- ((az + 1i * bz) / (ax + 1i * bx))^k
            if (abs(FT2) < eps) break
          }
        }
      }

      # Adjust the interval l if necessary based on the Fourier transform value
      if (abs(FT) > 1.001)
        if (l > 1) l <- r

      q1 <- -l
      q2 <- l
      q <- seq(q1, q2, delta)
    }

    # Final Fourier transform to obtain the estimated density y
    FTy <- numeric(length(q))
    if (error == "unknown") {
      for (ll in seq(along=q)) {
        l <- q[ll]
        ax <- Re(1/length(xr) * sum(exp(1i * l * xr)))
        bx <- Im(1/length(xr) * sum(exp(1i * l * xr)))
        az <- Re(1/length(zr) * sum(exp(1i * l * zr)))
        bz <- Im(1/length(zr) * sum(exp(1i * l * zr)))
        FTy[ll] <- ((az + 1i * bz) / (ax + 1i * bx))^k
      }
    }

    # Fourier transform for normal error model
    else if (error == "normal") {
      for (ll in seq(along=q)) {
        l <- q[ll]
        ax <- Re(exp(-1/2 * sigma^2 * (l * a)^2))
        bx <- Im(exp(-1/2 * sigma^2 * (l * a)^2))
        az <- Re(1/length(zr) * sum(exp(1i * l * zr)))
        bz <- Im(1/length(zr) * sum(exp(1i * l * zr)))
        FTy[ll] <- ((az + 1i * bz) / (ax + 1i * bx))^k
      }
    }

    # Fourier transform for Laplacian error model
    else if (error == "laplacian") {
      for (ll in seq(along=q)) {
        l <- q[ll]
        ax <- Re(1 / (1 + sigma^2/2 * (l * a)^2))
        bx <- Im(1 / (1 + sigma^2/2 * (l * a)^2))
        az <- Re(1/length(zr) * sum(exp(1i * l * zr)))
        bz <- Im(1/length(zr) * sum(exp(1i * l * zr)))
        FTy[ll] <- ((az + 1i * bz) / (ax + 1i * bx))^k
      }
    }
    ReFTy <- Re(FTy)
    ImFTy <- Im(FTy)

    # Determine the interval for y based on z and x
    if (error == "unknown") {
      s1 <- min(z) - max(x)
      s2 <- max(z) - min(x)
    } else {
      if (error == "normal") {
        x1 <- -6 * a * sigma
        x2 <- -x1
        x.x <- seq(x1, x2, length.out = Lx)
      } else if (error == "laplacian") {
        x1 <- -8 * a * sigma
        x2 <- -x1
        x.x <- seq(x1, x2, length.out = Lx)
      }
      s1 <- min(z) - max(x.x)
      s2 <- max(z) - min(x.x)
    }

    xy <- seq(s1, s2, length.out = Ly)
    y <- numeric(length(xy))
    for (ll in seq(along = xy)) {
      l <- xy[ll]
      aRe <- (q2 - q1) / length(q) * sum(ReFTy * cos(-l * q))
      bRe <- (q2 - q1) / length(q) * sum(ImFTy * sin(-l * q))
      y[ll] <- 1 / (2 * pi) * (aRe - bRe)
    }

    # Ensure positivity of the estimated density if required
    if (positive == TRUE) y[which(y < 0)] <- 0

    # Optional error calculation
    if (calc.error == FALSE) err <- NA
    else if (calc.error == TRUE) {
      if(error == "unknown"){
        x1C <- min(x)
        x2C <- max(x)
      }
      else if(error == "normal"){
        x1C <- min(x.x)
        x2C <- max(x.x)
      }
      else if(error == "laplacian"){
        x1C <- min(x.x)
        x2C <- max(x.x)
      }
      if(error == "unknown") fxC <- densprf(x, df = dfx)
      else if(error == "normal") fxC <- function(x) dnorm(x, 0, sigma)
      else if(error == "laplacian") fxC <- function(x) dlaplace(x, 0, sigma/sqrt(2))

      xC1 <- seq(min(xy), min(x1C, xy), -mean(diff(xy)))[-1]
      xC2 <- seq(max(xy), max(x2C, xy), mean(diff(xy)))[-1]
      xC <- c(xC1, xy, xC2)

      yC <- c(rep(0, length(xC1)), y, rep(0, length(xC2)))

      qC <- seq(min(z), max(z), length.out = Lz)
      zE <- numeric(length(qC))
      for (l in seq(qC)) {
        zE[l] <- (max(xC) - min(xC)) / length(xC) * sum(fxC(qC[l] - xC) * yC)
      }
      densprz <- densprf(z)(qC)
      err <- 10 * (max(qC) - min(qC)) / length(qC) * sum((densprz - zE)^2)

      # Optional plotting of f_x * f_y
      if (plot == TRUE) {
        old_par <- par(no.readonly = TRUE)
        on.exit(par(old_par))

        plot(NULL, xlim = range(qC), ylim = c(0, max(densprz, zE)), xlab = "x", ylab = "Density")
        lines(qC, densprz, col = "blue", lwd = 2)
        par(new = TRUE)
        lines(qC, zE, col = "orange", lwd = 3)
        if (legend == TRUE) legend("topright", legend = c(expression(f[z]), expression(f[x] * " * " * f[y]^{NPFD})),
                                   col = c("blue", "orange"), lwd = c(2, 2))
      }
    }
  }

  lis <- list(xy, y, k, err)
  names(lis) <- c("x", "y", "N", "error")
  class(lis) <- "deconvolve"
  return(lis)
}
#' @export
plot.deconvolve <- function(x, ...) {
  plot(x$x, x$y, type = "l", xlab = "x", ylab = "Density", lwd = 2, ...)
}
