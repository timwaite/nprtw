plotting_grid <- function(data,h) {
  n <- length(data$x)
  t <- seq(from=min(data$x),to=max(data$x),length=200)
  t <- sort(c(t, data$x-h, data$x+h, data$x+h-1e-6, data$x-h+1e-6))
  t <- t[ t>=min(data$x) & t<=max(data$x) ]
  # NB the purpose of this elaborate evaluation grid is so that any discontinuities
  # are shown accurately in plots
  # relies on the support of the kernel being [-1,1]
  t
}
