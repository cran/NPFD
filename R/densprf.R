#' Density Estimation Function
#'
#' @description
#'
#' This function estimates the density using a Poisson GLM with natural splines.
#'
#' @details
#' \code{densprf} is a modification of the \code{denspr} function from the \strong{siggenes} package.
#'
#' For more details, see the documentation in the \strong{siggenes} package.
#'
#' @param x Input data vector.
#' @param n.interval Number of intervals (optional).
#' @param df Degrees of freedom for the splines.
#' @param knots.mode Boolean to determine if quantiles should be used for knots.
#' @param type.nclass Method for determining number of classes.
#' @param addx Add \eqn{x} values (optional).
#'
#' @return The function \code{densprf(x)} returns a function that, for a given input \code{z}, computes the estimated density evaluated at the position values of \code{z} as a result.
#'
#' @importFrom stats predict
#' @importFrom KernSmooth dpih
#' @importFrom siggenes denspr
#'
#' @export
#'
#' @examples
#' # Set seed for reproducibility
#' set.seed(123)
#'
#' # Generate random data
#' z <- rnorm(1000)
#'
#' # Apply densprf function
#' f <- densprf(z)
#'
#' # Define sequences for evaluation
#' x1 <- seq(-4, 4, 0.5)
#' x2 <- seq(-5, 5, 0.1)
#'
#' # Evaluate the density function at specified points
#' f1 <- f(x1)
#' f2 <- f(x2)
densprf <- function(x, n.interval = NULL, df = 5, knots.mode = TRUE,
                    type.nclass = c("wand", "scott", "FD"), addx = FALSE) {

  # Helper function for connecting intervals
  connect <- function(L) {
    for(i in seq(L)) {
      force(L[[i]])
    }
    function(x, y) {
      X <- vector("list", length(L))
      for(i in seq(X)) {
        X[[i]] <- 1
      }
      for(i in 1:length(x)) {
        X[[i]] <- X[[i]] * ifelse(y <= x[i], 1, 0)
        X[[i+1]] <- X[[i+1]] * ifelse(y > x[i], 1, 0)
      }
      V <- list()
      for(i in seq(X)) {
        V[[i]] <- X[[i]] * L[[i]](y)
      }
      S <- 0
      for(i in seq(V)) {
        S <- S + V[[i]]
      }
      S
    }
  }

  # Calculate Quantiles for Knots
  getQuantiles <- function(n.knots, mode) {
    den <- (n.knots + 1) / 2
    tmp1 <- (1:floor(n.knots / 2)) * mode / den
    tmp2 <- 1 - (1 - mode) * (ceiling(n.knots / 2):1) / den
    c(tmp1, tmp2)
  }

  nclass.wand <- function(x) {
    ceiling(diff(range(x)) / KernSmooth::dpih(x))
  }

  # Force evaluation of inputs to ensure they are available within the function
  force(x)
  force(n.interval)
  force(df)
  force(knots.mode)
  force(type.nclass)
  force(addx)

  if (is.null(n.interval)) {
    type <- match.arg(type.nclass)
    if (type == "wand") {
      # Directly reference the function from KernSmooth
      FUN <- function(x, level = 1) {
        ceiling(diff(range(x)) / KernSmooth::dpih(x, level = level))
      }
    } else {
      # For other types, use match.fun as before
      FUN <- match.fun(paste("nclass", type, sep = "."))
    }
    n.interval <- FUN(x)  # Calculate the number of intervals
  }
  else type <- NULL

  # Define the main function that will be returned
  function(z) {
    breaks <- seq(min(x), max(x), length = n.interval + 1)  # Create equally spaced breakpoints
    valHist <- graphics::hist(x, breaks = breaks, plot = FALSE)  # Compute the histogram
    center <- valHist$mids  # Get the midpoints of the histogram bins
    counts <- valHist$counts  # Get the counts for each bin
    ids <- which(counts > 0)  # Find the bins that have non-zero counts
    x.mode <- center[which.max(counts)]  # Determine the mode of x

    if (knots.mode) {  # If knots.mode is TRUE, calculate quantiles for knots
      x.q <- mean(center <= x.mode)  # Calculate the proportion of centers less than or equal to the mode
      q.knots <- getQuantiles(df - 1, x.q)  # Calculate the quantiles based on df and mode proportion
    }
    center <- center[ids]  # Filter out the center points corresponding to non-zero counts

    if (knots.mode) {
      knots <- stats::quantile(center, q.knots)  # Determine the knots based on quantiles
      tmp <- ns.out <- splines::ns(center, knots = knots)  # Create a natural spline using the knots
    } else {
      tmp <- ns.out <- splines::ns(center, df = df)  # Create a natural spline with the specified degrees of freedom
    }

    class(tmp) <- "matrix"
    mat <- data.frame(Number = counts[ids], tmp)
    glm.out <- stats::glm(Number ~ ., data = mat, family = stats::poisson)  # Fit a Poisson GLM to the data
    scale <- sum(diff(breaks) * counts)  # Calculate the scale factor for density normalization
    newx <- predict(ns.out, z)  # Predict spline values for new data z
    class(newx) <- "matrix"  # Set the class of newx to "matrix"
    preds <- predict(glm.out, data.frame(newx), type = "response")  # Predict the response based on the GLM
    out <- (preds / scale)  # Normalize the predictions by the scale factor
    class(out) <- "denspr"  # Set the class of the output to "denspr"

    nullf <- function(x) 0  # Define a null function that returns 0 for any input
    minx <- min(siggenes::denspr(x, addx = TRUE)$x)  # Get the minimum value of x from the density estimate
    maxx <- max(siggenes::denspr(x, addx = TRUE)$x)  # Get the maximum value of x from the density estimate
    connect(list(nullf, function(x) out, nullf))(x = c(minx, maxx), z)  # Connect the intervals and return the final density estimate
  }
}

