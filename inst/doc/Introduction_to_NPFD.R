## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  fig.width = 7,
  fig.height = 5
)

## ----library setup------------------------------------------------------------
library(NPFD)
library(siggenes)
library(KernSmooth)
library(splines)
library(stats)
library(graphics)
library(VGAM)

## ----example setup------------------------------------------------------------
set.seed(123)
x <- rnorm(1000)
y <- rgamma(1000, 10, 2)
z <- x + y

independent.x <- rnorm(100)

fy.NPFD <- deconvolve(independent.x, z, calc.error = T, plot = T)
fy <- denspr(y, addx = T)

fy.NPFD$N
fy.NPFD$error

plot(NULL, xlim = range(y), ylim = c(0, max(fy$y, fy.NPFD$y)), xlab = "x", ylab = "density")
lines(fy, col = "blue", lwd = 2)
lines(fy.NPFD, col = "orange", lwd = 2)
legend("topright", legend = c(expression(f[y]), expression(f[y]^{NPFD})), 
       col = c("blue", "orange"), lwd = c(2, 2))

