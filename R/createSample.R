#' Create a Sample from a Centered Distribution
#'
#' This function creates a sample from a centered distribution based on replicates of mixed data.
#'
#' @param z1 A numeric vector where \eqn{z_1 = x_1 + y}.
#' @param z2 A numeric vector of the same length as \eqn{z_1} where \eqn{z_2 = x_2 + y}.
#' @return A numeric vector representing a sample from the centered distribution.
#' @export
#'
#' @examples
#' # Set seed for reproducibility
#' set.seed(123)
#'
#' # Generate random data
#' x1 <- rnorm(1000)
#' x2 <- rnorm(1000)
#' y <- rgamma(1000, 10, 2)
#' z1 <- x1 + y
#' z2 <- x2 + y
#'
#' # Use createSample to generate a sample
#' x <- createSample(z1, z2)
#'
#' # Perform density estimation
#' f.x <- stats::density(x, adjust = 1.5)
#' x.x <- f.x$x
#' f <- dnorm(x.x)
#'
#' # Plot the results
#' plot(NULL, xlim = range(f.x$x), ylim = c(0, max(f, f.x$y)), xlab = "x", ylab = "Density")
#' lines(x.x, f, col = "blue", lwd = 2)
#' lines(f.x, col = "orange", lwd = 2)
#' legend("topright", legend = c(expression(f), expression(f[x])), col = c("blue", "orange"), lwd = 2)
createSample <- function(z1, z2) {
  # Check if the vectors have the same length
  if (length(z1) != length(z2)) stop("Error: The vectors must have the same length.")

  # Calculate and return the result
  result <- (z1 - z2) / sqrt(2)
  return(result)
}
