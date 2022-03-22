# Computes the continuous ranked probability score (CRPS)
calc_crps <- function(cdf, x, lower = -Inf, upper = Inf) {
  -integrate(function(y, x) (cdf(y) - 1*(y >= x))^2, lower, upper, x = x)$value
}

generate_crps <- function(bs, obs) {
  cdf_list <- purrr::map(bs, cdf_factory) %>% 
    imap_dbl(function(cdf, t) calc_crps(cdf, obs[1, t]))
}

cdf_factory <- function(bs) {
  Xdf_factory(bs, pnorm)
}

pdf_factory <- function(bs) {
  Xdf_factory(bs, dnorm)
}


Xdf_factory <- function(bs, Xdf) {
  force(bs)
  if (ncol(bs$mu_t) == 1 && length(bs$p_t) != 1) {
    f <- function(x) {
      Xdf(x, mean = bs$mu_t[1, 1, 1], sd = sqrt(bs$Sigma_t[1, 1, 1, 1]))
    }
  }
  else {
    list2env(marginal_x(bs), environment())  # Unpack w, mu_list, Sigma_list
    if (nrow(mu_list) > 1) stop("Can only produce cdf for univariate variables")
    piece <- function(i, x) {
      w[i]*Xdf(x, mean = mu_list[1, i], sd = sqrt(Sigma_list[1, 1, i]))
    }
    f <- function(x) {
      purrr::map(seq_along(w), piece, x = x) %>% 
        purrr::reduce(`+`, .init = 0)
    }
  }
  return(f)
}