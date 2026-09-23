############################################################
# Threshold delta for
#
# G(z) = exp(-delta*(1-z) - gamma*(1-z)^alpha)
#
# For each admissible alpha, the threshold is approximated
# as the smallest delta for which p_0,...,p_N are all
# nonnegative, up to the specified numerical tolerance.
#
# If 0 < alpha < 1, the search is allowed to extend down to
#
#     delta = -1.1 * alpha * gamma.
#
# If alpha > 1, the search is restricted to delta >= 0.
############################################################


# ----------------------------------------------------------
# Compute the recurrence coefficients
#
# c_k = gamma * (-1)^k * alpha^{falling (k+1)} / k!
#
# They are computed recursively to avoid separately
# evaluating very large factorials.
# ----------------------------------------------------------

compute_coefficients <- function(alpha, gamma, N = 500) {
  
  coeff <- numeric(N)
  
  if (N >= 1) {
    coeff[1] <- -gamma * alpha * (alpha - 1)
  }
  
  if (N >= 2) {
    for (k in 2:N) {
      coeff[k] <- -((alpha - k) / k) * coeff[k - 1]
    }
  }
  
  coeff
}


# ----------------------------------------------------------
# Compute p_0,...,p_N from the recurrence
# ----------------------------------------------------------

compute_p <- function(alpha, gamma, delta, N = 500) {
  
  p <- numeric(N + 1)
  
  # p_0
  p[1] <- exp(-delta - gamma)
  
  coeff <- compute_coefficients(
    alpha = alpha,
    gamma = gamma,
    N = N
  )
  
  for (n in 0:(N - 1)) {
    
    # Contribution from k = 0
    rhs <- (delta + alpha * gamma) * p[n + 1]
    
    # Contributions from k = 1,...,n
    if (n >= 1) {
      rhs <- rhs + sum(
        coeff[1:n] * rev(p[1:n])
      )
    }
    
    p[n + 2] <- rhs / (n + 1)
  }
  
  p
}


# ----------------------------------------------------------
# Test whether p_0,...,p_N are nonnegative
#
# Values between -pmf_tol and 0 are treated as numerical
# round-off errors.
# ----------------------------------------------------------

is_legitimate <- function(alpha,
                          gamma,
                          delta,
                          N = 500,
                          pmf_tol = 1e-12) {
  
  p <- compute_p(
    alpha = alpha,
    gamma = gamma,
    delta = delta,
    N = N
  )
  
  if (any(!is.finite(p))) {
    return(FALSE)
  }
  
  all(p >= -pmf_tol)
}


# ----------------------------------------------------------
# Lower end of the delta search
#
# For alpha < 1:
#
#     delta >= -1.1 * alpha * gamma.
#
# For alpha > 1:
#
#     delta >= 0.
# ----------------------------------------------------------

delta_search_lower_bound <- function(alpha, gamma) {
  
  if (alpha < 1) {
    -1.1 * alpha * gamma
  } else {
    0
  }
}


# ----------------------------------------------------------
# Locate the threshold delta by bisection
#
# This assumes that the admissible values of delta form an
# interval [delta_min, infinity).
#
# If the PGF is already legitimate at the lower end of the
# prescribed search range, that lower endpoint is returned.
# Thus, for alpha < 1, a returned value of
#
#     -1.1 * alpha * gamma
#
# means that the true threshold may lie still further left.
# ----------------------------------------------------------

find_delta_min <- function(alpha,
                           gamma,
                           N = 500,
                           delta_tol = 1e-8,
                           pmf_tol = 1e-12,
                           max_delta = 1e6) {
  
  lo <- delta_search_lower_bound(alpha, gamma)
  
  # If the lower endpoint is already admissible, return it.
  if (is_legitimate(
    alpha = alpha,
    gamma = gamma,
    delta = lo,
    N = N,
    pmf_tol = pmf_tol
  )) {
    return(lo)
  }
  
  # Find an admissible upper endpoint.
  hi <- max(1, lo + 1)
  
  while (!is_legitimate(
    alpha = alpha,
    gamma = gamma,
    delta = hi,
    N = N,
    pmf_tol = pmf_tol
  )) {
    
    hi <- max(2 * hi, hi + 1)
    
    if (hi > max_delta) {
      warning(
        paste(
          "No admissible delta found for alpha =",
          format(alpha, digits = 8)
        )
      )
      
      return(NA_real_)
    }
  }
  
  # At this point:
  #
  #   lo is not admissible,
  #   hi is admissible.
  #
  # Bisect until the desired delta accuracy is reached.
  
  while ((hi - lo) > delta_tol) {
    
    mid <- (lo + hi) / 2
    
    if (is_legitimate(
      alpha = alpha,
      gamma = gamma,
      delta = mid,
      N = N,
      pmf_tol = pmf_tol
    )) {
      hi <- mid
    } else {
      lo <- mid
    }
  }
  
  hi
}


# ----------------------------------------------------------
# Determine the admissible open alpha intervals
#
# gamma > 0:
#
#     (0,1), (2,3), (4,5)
#
# gamma < 0:
#
#     (1,2), (3,4), (5,6)
#
# A small endpoint offset is used to avoid evaluating exactly
# at the integer endpoints.
# ----------------------------------------------------------

alpha_intervals <- function(gamma, endpoint_offset = 0.01) {
  
  if (endpoint_offset <= 0 || endpoint_offset >= 0.5) {
    stop("endpoint_offset must lie strictly between 0 and 0.5.")
  }
  
  if (gamma > 0) {
    
    list(
      c(0 + endpoint_offset, 1 - endpoint_offset),
      c(2 + endpoint_offset, 3 - endpoint_offset),
      c(4 + endpoint_offset, 5 - endpoint_offset)
    )
    
  } else if (gamma < 0) {
    
    list(
      c(1 + endpoint_offset, 2 - endpoint_offset),
      c(3 + endpoint_offset, 4 - endpoint_offset),
      c(5 + endpoint_offset, 6 - endpoint_offset)
    )
    
  } else {
    
    stop(
      paste(
        "gamma = 0 is degenerate:",
        "G(z) is then the PGF of a Poisson distribution",
        "whenever delta >= 0."
      )
    )
  }
}


# ----------------------------------------------------------
# Compute one branch of the threshold curve
# ----------------------------------------------------------

compute_branch <- function(alpha_lo,
                           alpha_hi,
                           gamma,
                           n_alpha = 100,
                           N = 500,
                           delta_tol = 1e-8,
                           pmf_tol = 1e-12) {
  
  alpha_grid <- seq(
    from = alpha_lo,
    to = alpha_hi,
    length.out = n_alpha
  )
  
  delta_grid <- vapply(
    alpha_grid,
    FUN = find_delta_min,
    FUN.VALUE = numeric(1),
    gamma = gamma,
    N = N,
    delta_tol = delta_tol,
    pmf_tol = pmf_tol
  )
  
  data.frame(
    alpha = alpha_grid,
    delta = delta_grid
  )
}


# ----------------------------------------------------------
# Compute and plot all three branches
#
# If add = FALSE, a new plot with axes is created.
#
# If add = TRUE, the curves are added to the current plot
# without creating new axes or changing the plotting region.
#
# All three branches are plotted in the same colour.
#
# The returned value is a list of three data frames, one for
# each branch.
# ----------------------------------------------------------

plot_threshold_curve <- function(gamma,
                                 n_alpha = 100,
                                 N = 500,
                                 delta_tol = 1e-8,
                                 pmf_tol = 1e-12,
                                 endpoint_offset = 0.01,
                                 curve_colour = "blue",
                                 curve_width = 2,
                                 curve_type = 1,
                                 add = FALSE) {
  
  intervals <- alpha_intervals(
    gamma = gamma,
    endpoint_offset = endpoint_offset
  )
  
  branches <- lapply(
    intervals,
    function(interval) {
      compute_branch(
        alpha_lo = interval[1],
        alpha_hi = interval[2],
        gamma = gamma,
        n_alpha = n_alpha,
        N = N,
        delta_tol = delta_tol,
        pmf_tol = pmf_tol
      )
    }
  )
  
  all_deltas <- unlist(
    lapply(
      branches,
      function(branch) branch$delta
    )
  )
  
  finite_deltas <- all_deltas[is.finite(all_deltas)]
  
  if (length(finite_deltas) == 0) {
    stop("No finite threshold values were obtained.")
  }
  
  # Create a new plotting region only when add = FALSE.
  if (!add) {
    
    y_min <- min(0, finite_deltas)
    y_max <- max(0, finite_deltas)
    
    # Avoid a zero-height plotting range.
    if (y_min == y_max) {
      y_min <- y_min - 1
      y_max <- y_max + 1
    }
    
    # Add a small amount of vertical padding.
    y_padding <- 0.04 * (y_max - y_min)
    
    plot(
      NA,
      xlim = c(0, 6),
      ylim = c(
        y_min - y_padding,
        y_max + y_padding
      ),
      xlab = expression(alpha),
      ylab = expression(delta[min]),
      main = paste(
        "Threshold delta for gamma =",
        gamma
      ),
      xaxt = "n"
    )
    
    axis(
      side = 1,
      at = 0:6
    )
    
    # Draw a horizontal reference line at delta = 0.
    abline(
      h = 0,
      col = "grey75",
      lty = 2
    )
    
    grid()
  }
  
  # Draw all three branches in the same colour.
  for (branch in branches) {
    
    finite <- is.finite(branch$delta)
    
    lines(
      x = branch$alpha[finite],
      y = branch$delta[finite],
      col = curve_colour,
      lwd = curve_width,
      lty = curve_type
    )
  }
  
  invisible(branches)
}


############################################################
# USER SETTINGS
############################################################

# Set gamma to any nonzero value between -3 and 3.
#
# If gamma > 0, the alpha intervals are:
#     (0,1), (2,3), (4,5).
#
# If gamma < 0, the alpha intervals are:
#     (1,2), (3,4), (5,6).

gamma <- 1


############################################################
# COMPUTE AND PLOT
############################################################

results <- plot_threshold_curve(
  gamma = -0.05,
  n_alpha = 100,
  N = 300,
  delta_tol = 1e-6,
  pmf_tol = 1e-12,
  endpoint_offset = 0.01,
  curve_colour = "grey",
  curve_width = 3,
  add = TRUE
)


############################################################
# OPTIONAL: inspect or export the computed values
############################################################

# The three branches are stored separately:
#
# results[[1]]
# results[[2]]
# results[[3]]
#
# For example:
#
# print(results[[1]])
#
# To combine them into one data frame:
#
# combined_results <- do.call(
#   rbind,
#   lapply(
#     seq_along(results),
#     function(i) {
#       data.frame(
#         branch = i,
#         results[[i]]
#       )
#     }
#   )
# )
#
# To save the combined values:
#
# write.csv(
#   combined_results,
#   file = "delta_thresholds.csv",
#   row.names = FALSE
# )
#
# Example with negative gamma:
#
# gamma <- -1
#
# results <- plot_threshold_curve(
#   gamma = gamma,
#   n_alpha = 100,
#   N = 500,
#   curve_colour = "blue"
# )
############################################################